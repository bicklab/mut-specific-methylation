version 1.0

workflow run_logistic_ewas {
  input {
    File beta_matrix
    File meta_data
    String output_file
    String outcome
    String covariates
    File logistic_ewas_R_script
    String sample_column = "sample_id"
    Int chunk_size = 50000
    Int cpu = 16
    Float memory = 104
    Int disk_size = 200
    Int preemptible = 0
    Int maxRetries = 1
  }

  call ewas_logistic {
    input:
      beta_matrix = beta_matrix,
      meta_data = meta_data,
      output_file = output_file,
      outcome = outcome,
      covariates = covariates,
      logistic_ewas_R_script = logistic_ewas_R_script,
      sample_column = sample_column,
      chunk_size = chunk_size,
      cpu = cpu,
      memory = memory,
      disk_size = disk_size,
      preemptible = preemptible,
      maxRetries = maxRetries
  }

  output {
    File ewas_results = ewas_logistic.ewas_summary_stats
  }

  meta {
    author: "Yash Pershad"
    email: "yash.pershad@vanderbilt.edu"
    description: "This WDL runs EWAS as a logistic regression for a binary trait."
  }
}

task ewas_logistic {
  input {
    File beta_matrix
    File meta_data
    String output_file
    String outcome
    String covariates
    File logistic_ewas_R_script
    String sample_column
    Int chunk_size
    Int cpu
    Float memory
    Int disk_size
    Int preemptible
    Int maxRetries
  }

  command <<<
    Rscript ~{logistic_ewas_R_script} \
      ~{beta_matrix} \
      ~{meta_data} \
      ~{output_file} \
      ~{outcome} \
      ~{covariates} \
      ~{sample_column} \
      ~{chunk_size}
  >>>

  output {
    File ewas_summary_stats = "~{output_file}"
  }

  runtime {
    docker: "hpoisner/r_phewas:latest"
    memory: memory + " GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    preemptible: preemptible
    maxRetries: maxRetries
  }
}