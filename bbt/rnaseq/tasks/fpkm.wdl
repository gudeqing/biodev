version development

task fpkm{
    input {
        String? other_parameters
        File bam
        File gtf
        File info
        String sample_id = "sample_name"
        # for runtime
        String docker = "?"
        String memory = "10 GiB"
        Int cpu = 2
        String disks = "10 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e 
        FPKM-UQ.py \
        ~{other_parameters} \
        ~{"--bam " + bam} \
        ~{"--gtf " + gtf} \
        ~{"--info " + info} \
        ~{"--output " + sample_id} 
    >>>

    output {
        Array[File] outputs = glob(".*")
    }

    runtime {
        docker: docker
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "fpkm"
        image: "dockerfile"
        desc: "Calculate count, FPKM, and FPKM-UQ values defined by TCGA"
        version: "4.0.0"
        source: "http://rseqc.sourceforge.net/"
        basecmd: "FPKM-UQ.py"
    }

    parameter_meta {
        other_parameters: {desc: "other arguments, you could set any other argument with a string such as '-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        bam: {desc: "Alignment file in BAM format. BAM file should be sorted and indexed.", level: "required", type: "infile", range: "", default: ""}
        gtf: {desc: "Gene model in GTF format", level: "required", type: "infile", range: "", default: ""}
        info: {desc: "Gene model information file.", level: "required", type: "infile", range: "", default: ""}
        sample_id: {desc: "prefix for output file name", level: "required", type: "str", range: "", default: "sample_name"}
    }

}
