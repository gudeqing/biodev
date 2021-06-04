version 1.0

task geneBodyCoverage{
    input {
        String? other_parameters
        Array[File] bam
        File refGene
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
        geneBody_coverage.py \
        ~{other_parameters} \
        ~{if defined(bam) then "-i " else ""}~{sep="," bam} \
        ~{"-r " + refGene} \
        ~{"-o " + sample_id} 
    >>>

    output {
        Array[File] outputs = glob("*")
    }

    runtime {
        docker: docker
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "geneBodyCoverage"
        image: "dockerfile"
        desc: "Calculate the RNA-seq reads coverage over gene body."
        version: "4.0.0"
        source: "http://rseqc.sourceforge.net/"
        basecmd: "geneBody_coverage.py"
    }

    parameter_meta {
        other_parameters: {desc: "other arguments, you could set any other argument with a string such as '-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        bam: {desc: "Alignment file in BAM format. BAM file should be sorted and indexed.", level: "required", type: "infile", range: "", default: ""}
        refGene: {desc: "Reference gene model in bed format.", level: "required", type: "infile", range: "", default: ""}
        sample_id: {desc: "prefix for output file name", level: "required", type: "str", range: "", default: "sample_name"}
    }

}
