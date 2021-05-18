version 1.0

task rnaseqc{
    input {
        String? other_parameters
        File gtf
        File bam
        String outdir = "."
        String sample_id
        File? bed
        String? strand
        # for runtime
        String docker = "gcr.io/broad-cga-aarong-gtex/rnaseqc:latest"
        String memory = "6 GiB"
        Int cpu = 2
        String disks = "6 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e 
        rnaseqc \
        ~{other_parameters} \
        ~{gtf} \
        ~{bam} \
        ~{outdir} \
        ~{"--sample " + sample_id} \
        ~{"--bed " + bed} \
        ~{"--stranded " + strand} 
    >>>

    output {
        File gene_tpm = "${sample_id}.gene_tpm.gct"
        File gene_counts =  "${sample_id}.gene_reads.gct"
        File exon_counts =  "${sample_id}.exon_reads.gct"
        File metrics = "${sample_id}.metrics.tsv"
        File insertsize_distr =  "${sample_id}.fragmentSizes.txt"
    }

    runtime {
        docker: docker
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "rnaseqc"
        image: "gcr.io/broad-cga-aarong-gtex/rnaseqc:latest"
        desc: "Fast efficient RNA-Seq metrics for quality control and process optimization"
        version: "2.4.2"
        source: "https://github.com/getzlab/rnaseqc"
        basecmd: "rnaseqc"
    }

    parameter_meta {
        other_parameters: {desc: "other arguments, you could set any other argument with a string such as '-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        gtf: {desc: "The input GTF file containing features to check the bam against", level: "required", type: "infile", range: "", default: ""}
        bam: {desc: "The input SAM/BAM file containing reads", level: "required", type: "infile", range: "", default: ""}
        outdir: {desc: "Output directory", level: "required", type: "str", range: "", default: "."}
        sample_id: {desc: "prefix for output file name", level: "optional", type: "str", range: "", default: ""}
        bed: {desc: "Optional input BED file containing non-overlapping exons used for fragment size calculations", level: "optional", type: "infile", range: "", default: ""}
        strand: {desc: "Use strand-specific metrics. Only features on the same strand of a read will be considered.", level: "optional", type: "str", range: "RF,FR", default: ""}
    }

}
