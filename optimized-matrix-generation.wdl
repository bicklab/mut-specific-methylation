version 1.0

# Main workflow for creating methylation matrix from CX_report files
workflow optimized_matrix_generation {
    input {
        String directory
        Float min_sample_fraction = 0.8
        Int high_coverage_threshold = 25
        Int min_coverage_threshold = 10
    }

    # Step 1: Find all CX_report files
    call find_cx_files {
        input:
            directory = directory
    }

    # Step 2: Process each CX_report file to extract 10x coverage sites
    scatter(idx in range(length(find_cx_files.cx_files))) {
        call extract_10x_sites {
            input:
                cx_file = find_cx_files.cx_files[idx],
                min_coverage = min_coverage_threshold,
                high_coverage = high_coverage_threshold
        }
    }

    # Collect all scatter outputs
    Array[File] site_list_files = select_all(extract_10x_sites.site_list)
    Array[File] high_coverage_site_files = select_all(extract_10x_sites.high_coverage_sites)
    Array[File] filtered_data_files = select_all(extract_10x_sites.filtered_data)
    Array[String] sample_names = select_all(extract_10x_sites.sample_names)

    # Step 3: Combine all site lists to find common sites
    call find_common_sites {
        input:
            high_coverage_files = high_coverage_site_files,
            sample_count = length(find_cx_files.cx_files),
            min_sample_fraction = min_sample_fraction   
    }

    # Step 4: Extract methylation data for selected sites from each sample
    scatter(idx in range(length(filtered_data_files))) {
        call extract_methylation_data {
            input:
                filtered_data = filtered_data_files[idx],
                selected_sites = find_common_sites.selected_sites,
                sample_name = sample_names[idx]
            }
    }

    # Collect methylation data outputs
    Array[File] sorted_methylation_files = select_all(extract_methylation_data.sorted_data)

    call create_methylation_matrix {
        input:
            sorted_methylation_files = sorted_methylation_files
    }

    output {
        File beta_matrix = create_methylation_matrix.beta_matrix
    }
}

# Task 1: Find all CX_report files in directory
task find_cx_files {
    input {
        String directory
        Int cpu = 1
        Float memory = 2
    }

    command <<<
        set -eu -o pipefail
        
        # List all CX_report.txt.gz files
        gsutil ls ~{directory}*.CX_report.txt.gz > cx_files.txt
        
        echo "Found $(wc -l < cx_files.txt) CX_report files"
    >>>
    
    output {
        Array[File] cx_files = read_lines("cx_files.txt")
    }

    runtime {
        docker: "google/cloud-sdk:latest"
        memory: memory + " GB"
        cpu: cpu
    }
}

# Task 2: Extract sites with ≥10x coverage from each CX_report file
task extract_10x_sites {
    input {
        File cx_file
        Int min_coverage = 10
        Int high_coverage = 25
        Float memory = 8
        Int cpu = 2
        Int preemptible = 2
        Int maxRetries = 3
    }

    Int disk_size = ceil(20 + 3 * size(cx_file, "GiB"))
    String sample_name = basename(cx_file, ".CX_report.txt.gz")
    
    command <<<
        set -eu -o pipefail

        echo "Processing sample: ~{sample_name}"
        echo "Minimum coverage: ~{min_coverage}x"
        echo "High coverage threshold: ~{high_coverage}x"

        # Extract CG sites with ≥10x coverage and create outputs
        zcat ~{cx_file} | \
        awk -v min_cov="~{min_coverage}" -v high_cov="~{high_coverage}" -v sample="~{sample_name}" '
            BEGIN {
                OFS="\t"; 
                sites_10x=0; 
                sites_25x=0; 
                total_cov=0
            } 
            $6 == "CG" && $4 + $5 >= min_cov {
                chr=$1; pos=$2; strand=$3; 
                coverage=$4+$5; meth=$4;
                methylation_pct = (coverage > 0) ? (meth/coverage)*100 : 0;
                site_id = chr ":" pos;
                
                # Output full data for matrix creation
                print site_id, chr, pos, strand, coverage, methylation_pct > sample "_filtered_data.txt";
                
                # Output just site IDs for common site detection
                print site_id > sample "_sites_10x.txt";
                
                # Track high coverage sites
                if(coverage >= high_cov) {
                    print site_id > sample "_sites_25x.txt";
                    sites_25x++;
                }
                sites_10x++;
                total_cov += coverage;
            }
            END {
                mean_cov = (sites_10x > 0) ? total_cov/sites_10x : 0;
                print sample "," sites_10x "," sites_25x "," mean_cov > sample "_stats.txt";
            }'

        # Ensure high coverage file exists
        touch ~{sample_name}_sites_25x.txt

        echo "Sample ~{sample_name} processed:"
        echo "  Sites ≥10x: $(wc -l < ~{sample_name}_sites_10x.txt)"
        echo "  Sites ≥25x: $(wc -l < ~{sample_name}_sites_25x.txt)"
    >>>

    output {
        File site_list = "~{sample_name}_sites_10x.txt"
        File high_coverage_sites = "~{sample_name}_sites_25x.txt"
        File filtered_data = "~{sample_name}_filtered_data.txt"
        File sample_stats = "~{sample_name}_stats.txt"
        String sample_names = sample_name
    }

    runtime {
        docker: "gcr.io/nygc-public/genome-utils:v8"
        memory: memory + " GB"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

# Task 3: Find sites common across specified fraction of samples
task find_common_sites {
    input {
        Array[File] high_coverage_files
        Int sample_count
        Float min_sample_fraction
        Int cpu = 8
        Float memory = 30
        Int preemptible = 0
    }

    Int min_samples_required = ceil(min_sample_fraction * sample_count) - 1
    Int disk_size = ceil(100 + 5*size(high_coverage_files, "GiB"))

    command <<<
        set -eu -o pipefail

        echo "=== FINDING COMMON SITES ==="
        echo "Total samples: ~{sample_count}"
        echo "Minimum sample fraction: ~{min_sample_fraction}"
        echo "Minimum samples required: ~{min_samples_required}"

        # Set minimum samples threshold
        min_samples_required=~{min_samples_required}

        # Count site occurrences for 25x coverage
        echo "Counting sites with ≥25x coverage across samples..."

        # Use the WDL file arrays directly
        echo "Processing ~{sample_count} high coverage files"

        # Concatenate all high coverage files, sort, count occurrences, and filter
        cat ~{sep=' ' high_coverage_files} | sort | uniq -c | \
            awk -v min_samples="$min_samples_required" '$1 > min_samples {print $2}' > selected_sites_25x.txt

        sites_25x_count=$(wc -l < selected_sites_25x.txt)
        echo "Sites found in >$min_samples_required samples: $sites_25x_count"
    >>>

    output {
        File selected_sites = "selected_sites_25x.txt"
    }

    runtime {
        docker: "gcr.io/nygc-public/genome-utils:v8"
        memory: memory + " GB"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}

# Task 4: Extract methylation data for selected sites from each sample
task extract_methylation_data {
    input {
        File filtered_data
        File selected_sites
        String sample_name
        Float memory = 4
        Int cpu = 2
        Int preemptible = 0
    }

    Int disk_size = ceil(5 + size(filtered_data, "GiB") + size(selected_sites, "GiB"))

    command <<<
        set -eu -o pipefail

        echo "Extracting methylation data for sample: ~{sample_name}"
        awk 'NR==FNR{sites[$1]=1; next} $1 in sites' ~{selected_sites} ~{filtered_data} > ~{sample_name}_methylation_data.txt
        echo "Extracted $(wc -l < ~{sample_name}_methylation_data.txt) sites for ~{sample_name}"

        echo "Sorting ~{sample_name}..."
        sort -k1,1 ~{sample_name}_methylation_data.txt > ~{sample_name}_sorted.txt
    >>>

    output {
        File methylation_data = "~{sample_name}_methylation_data.txt"
        File sorted_data = "~{sample_name}_sorted.txt"
    }

    runtime {
        docker: "gcr.io/nygc-public/genome-utils:v8"
        memory: memory + " GB"
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}

# Task 5: Create final methylation matrix (OPTIMIZED VERSION)
task create_methylation_matrix {
    input {
        Array[File] sorted_methylation_files
        Float memory = 30
        Int cpu = 8
        Int preemptible = 0
    }

    Int disk_size = ceil(50 + 2 * size(sorted_methylation_files, "GiB"))

    command <<<
        set -eu -o pipefail

        echo "=== CREATING BETA VALUE MATRIX FOR EWAS (OPTIMIZED) ==="
        
        # Extract sample names from file paths
        echo "Extracting sample names from file paths..."
        for file in ~{sep=' ' sorted_methylation_files}; do
            # Extract basename and remove common extensions
            sample_name=$(basename "$file" | sed 's/_sorted\.\(txt\|tsv\|gz\)$//' | sed 's/\.txt$//' | sed 's/\.tsv$//' | sed 's/\.gz$//')
            echo "$sample_name" >> sample_names.txt
        done
        
        echo "Sample names extracted:"
        cat sample_names.txt
                
        # Get site coordinates from first file (assuming all files have same sites in same order)
        first_file=$(echo ~{sep=' ' sorted_methylation_files} | cut -d' ' -f1)
        echo "Using $first_file to establish site order"
        awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' "$first_file" > site_coordinates.txt
        echo "Found $(wc -l < site_coordinates.txt) sites"
        
        # Extract beta values from each file in parallel (much faster than sequential)
        echo "Extracting beta values from all files in parallel..."
        sample_idx=0
        for file in ~{sep=' ' sorted_methylation_files}; do
            sample_idx=$((sample_idx + 1))
            echo "Queuing sample $sample_idx: $(basename $file)"
        done
        
        # Use xargs to process files in parallel
        echo ~{sep=' ' sorted_methylation_files} | tr ' ' '\n' | \
        nl -nln | \
        xargs -P ~{cpu} -I {} sh -c '
            set -- {}
            idx=$1
            file=$2
            echo "Processing sample $idx in parallel: $(basename $file)"
            awk "{print \$6/100}" "$file" > "beta_values_$idx.txt"
        '

        echo "All beta value extraction completed"
        
        # Create the matrix using paste (much faster than iterative joins)
        echo "Creating matrix with paste..."
        
        {
            # Header line
            printf "site_id\tchr\tpos\tstrand"
            while read sample; do
                printf "\t%s" "$sample"
            done < sample_names.txt
            printf "\n"
            
            # Data lines - paste site coordinates with all beta value columns
            paste site_coordinates.txt beta_values_*.txt
            
        } > beta_matrix.tsv
        
        # Cleanup
        rm -f beta_values_*.txt site_coordinates.txt sample_names.txt
        
        echo "Beta matrix created: beta_matrix.tsv"
        echo "Matrix dimensions: $(wc -l < beta_matrix.tsv) rows x $(head -1 beta_matrix.tsv | tr '\t' '\n' | wc -l) columns"
    >>>

    output {
        File beta_matrix = "beta_matrix.tsv"
    }

    runtime {
        docker: "gcr.io/nygc-public/genome-utils:v8"
        memory: memory + " GB"
        disks: "local-disk " + disk_size + " HDD"  # Use SSD for better I/O
        cpu: cpu
        preemptible: preemptible
    }
}