version 1.0

task fastp{
    input {
        String? other_parameters
        Int threads = 4
        File read1
        File read2
        String sample_name
        String? adapter_r1
        String? adapter_r2
        # for runtime
        String docker = "gudeqing/fastp:0.21.0"
        String memory = "5 GiB"
        Int cpu = 4
        String disks = "5 GiB"
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
        ~{"-o " + sample_name + ".clean.R1.fastq.gz"} \
        ~{if defined(read2) then "-O " + sample_name + ".clean.R2.fastq.gz" else ""} \
        ~{"-h " + sample_name + ".report.html"} \
        ~{"-j " + sample_name + ".report.json"}
    >>>

    output {
        File out_read1_file = sample_name + ".clean.R1.fastq.gz"
        File? out_read2_file = sample_name + ".clean.R2.fastq.gz"
        File html_report_file = sample_name + ".report.html"
        File json_report_file = sample_name + ".report.json"
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
        sample_name: {desc: "sample name, will be used in output files' name", level: "required", type: "str", range: "", default: "sample_name"}
        adapter_r1: {desc: "the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])", level: "optional", type: "str", range: "", default: "auto"}
        adapter_r2: {desc: "the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=auto])", level: "optional", type: "str", range: "", default: "auto"}
    }

}
