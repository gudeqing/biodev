version development

task read_distribution{
    input {
        String? other_parameters
        File bam
        File refGene
        String sample_id = "sample_name"
        # for runtime
        String docker = "cdiasgurjao/rseqc:latest"
        String memory = "10 GiB"
        Int cpu = 2
        String disks = "10 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e 
        read_distribution.py \
        ~{other_parameters} \
        ~{"-i " + bam} \
        ~{"-r " + refGene} \
        ~{"> " + sample_id + ".read_distribution.txt"}
    >>>

    output {
        File read_distr = "${sample_id}.read_distribution.txt"
    }

    runtime {
        docker: docker
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "read_distribution"
        image: "dockerfile"
        desc: "Provided a BAM/SAM file and reference gene model, this module will calculate how mapped reads were distributed over genome feature (like CDS exon, 5’UTR exon, 3’ UTR exon, Intron, Intergenic regions)."
        version: "?"
        source: "http://rseqc.sourceforge.net/"
        basecmd: "read_distribution.py"
    }

    parameter_meta {
        other_parameters: {desc: "other arguments, you could set any other argument with a string such as '-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        bam: {desc: "Alignment file in BAM format. BAM file should be sorted and indexed.", level: "required", type: "infile", range: "", default: ""}
        refGene: {desc: "Reference gene model in bed format.", level: "required", type: "infile", range: "", default: ""}
        sample_id: {desc: "prefix for output file name", level: "required", type: "str", range: "", default: "sample_name"}
    }

}
