version 1.0

task fastp{
    input {
        String? other_parameters
        Int threads = 2
        File read1 = "sample.R1.fastq.gz"
        File? read2 = "sample.R2.fastq.gz"
        String? adapter_r1 = "auto"
        String? adapter_r2 = "auto"
        String out_read1 = "sample.clean.R1.fastq.gz"
        String? out_read2 = "sample.clean.R2.fastq.gz"
        String html_report = "sample.html"
        String json_report = "sample.json"
        # for runtime
        String docker = "fastp:latest"
        String memory = "1 GiB"
        Int cpu = 1
        String disks = "1 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e 
        ./fastp \
        ~{other_parameters} \
        ~{"--thread " + threads} \
        ~{"-i " + read1} \
        ~{"-I " + read2} \
        ~{"--adapter_sequence " + adapter_r1} \
        ~{"--adapter_sequence_r2 " + adapter_r2} \
        ~{"-o " + out_read1} \
        ~{"-O " + out_read2} \
        ~{"-h " + html_report} \
        ~{"-j " + json_report} 
    >>>

    output {
        File out_read1_file = "~{out_read1}"
        File out_read2_file = "~{out_read2}"
        File html_report_file =  "~{html_report}"
        File json_report_file =  "~{json_report}"
    }

    runtime {
        docker: docker
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "fastp"
        docker: "Dockerfile"
        desc: "A tool designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with multithreading supported to afford high performance. For more detail please refer to https://github.com/OpenGene/fastp"
        logo: "see fastp.png"
        version: "0.21.0"
        basecmd: "./fastp"
    }

    parameter_meta {
        other_parameters: {desc: "其他参数，你可以通过该参数输入一个或多个任何其他当前软件支持的参数，例如'-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        threads: {desc: "Number of threads to use", level: "required", type: "int", range: "", default: "2"}
        read1: {desc: "read1 fastq file", level: "required", type: "infile", range: "", default: "sample.R1.fastq.gz"}
        read2: {desc: "read2 fastq file", level: "optional", type: "infile", range: "", default: "sample.R2.fastq.gz"}
        adapter_r1: {desc: "the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])", level: "optional", type: "str", range: "", default: "auto"}
        adapter_r2: {desc: "the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=auto])", level: "optional", type: "str", range: "", default: "auto"}
        out_read1: {desc: "output read1 fastq file", level: "required", type: "str", range: "", default: "sample.clean.R1.fastq.gz"}
        out_read2: {desc: "output read2 fastq file", level: "optional", type: "str", range: "", default: "sample.clean.R2.fastq.gz"}
        html_report: {desc: "output html report file", level: "required", type: "str", range: "", default: "sample.html"}
        json_report: {desc: "output json report file", level: "required", type: "str", range: "", default: "sample.json"}
    }

}
