version 1.0

workflow run_logistic_dmr_ewas {
  input {
    File beta_matrix
    File meta_data
    String output_file
    String outcome
    String covariates
    File logistic_dmr_ewas_R_script
    String sample_column = "sample_id"
    Int window_size = 2000
    Int step_size = 1000
    Int min_cpgs_per_window = 5
    String summary_method = "mean"
    Int cpu = 16
    Float memory = 128
    Int disk_size = 200
    Int preemptible = 0
    Int maxRetries = 1
  }

  call ewas_logistic_dmr {
    input:
      beta_matrix = beta_matrix,
      meta_data = meta_data,
      output_file = output_file,
      outcome = outcome,
      covariates = covariates,
      logistic_dmr_ewas_R_script = logistic_dmr_ewas_R_script,
      sample_column = sample_column,
      window_size = window_size,
      step_size =  step_size,
      min_cpgs_per_window = min_cpgs_per_window,
      summary_method = summary_method,
      cpu = cpu,
      memory = memory,
      disk_size = disk_size,
      preemptible = preemptible,
      maxRetries = maxRetries
  }

  output {
    File ewas_results = ewas_logistic_dmr.ewas_summary_stats
  }

  meta {
    author: "Yash Pershad"
    email: "yash.pershad@vanderbilt.edu"
    description: "This WDL runs EWAS as a logistic regression for a binary trait."
  }
}

task ewas_logistic_dmr {
  input {
    File beta_matrix
    File meta_data
    String output_file
    String outcome
    String covariates
    File logistic_dmr_ewas_R_script
    String sample_column
    Int window_size = 2000
    Int step_size = 1000
    Int min_cpgs_per_window = 5
    String summary_method = "mean"
    Int cpu
    Float memory
    Int disk_size
    Int preemptible
    Int maxRetries
  }

  command <<<
    Rscript ~{logistic_dmr_ewas_R_script} \
      ~{beta_matrix} \
      ~{meta_data} \
      ~{output_file} \
      ~{outcome} \
      ~{covariates} \
      ~{sample_column} \
      ~{window_size} \
      ~{step_size} \
      ~{min_cpgs_per_window} \
      ~{summary_method}
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