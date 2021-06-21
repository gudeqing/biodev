version development

task CollectRnaSeqMetrics{
    input {
        String? other_parameters
        File input_bam
        String sample_id
        String strand = "NONE"
        File ref_flat
        File? ribosomal_intervals
        # for runtime
        String docker = "trinityctat/starfusion:1.10.0"
        String memory = "6 GiB"
        Int cpu = 2
        String disks = "6 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e 
        jar -jar /usr/local/src/picard.jar CollectRnaSeqMetrics \
        ~{other_parameters} \
        ~{"-I " + input_bam} \
        ~{"-O " + sample_id + ".RnaSeqMetrics.txt"} \
        ~{"--STRAND_SPECIFICITY " + strand} \
        ~{"--REF_FLAT " + ref_flat} \
        ~{"--RIBOSOMAL_INTERVALS " + ribosomal_intervals} \
        ~{"--CHART_OUTPUT " + sample_id + ".coverage.pdf"}
    >>>

    output {
        File rnaseq_metrics = "~{sample_id}.RnaSeqMetrics.txt"
        File coverage_pdf = "~{sample_id}.coverage.pdf"
    }

    runtime {
        docker: docker
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "CollectRnaSeqMetrics"
        docker: "broadinstitute/picard:latest"
        desc: "Produces RNA alignment metrics for a SAM or BAM file"
        logo: "none"
        version: "2.23.3-1-g4ac48fc-SNAPSHOT"
        basecmd: "java -jar /usr/picard/picard.jar CollectRnaSeqMetrics"
    }

    parameter_meta {
        input_bam: {desc: "One or more input SAM or BAM files to analyze. Must be coordinate sorted", level: "required", type: "infile", range: "", default: ""}
        strand: {desc: "For strand-specific library prep. For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.  Required. Possible values: {NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND}", level: "required", type: "str", range: "NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND", default: "NONE"}
        ref_flat: {desc: "Gene annotations in refFlat form.  Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat  Required.", level: "required", type: "infile", range: "", default: ""}
        ribosomal_intervals: {desc: "Location of rRNA sequences in genome, in interval_list format.  If not specified no bases will be identified as being ribosomal.", level: "optional", type: "infile", range: "", default: ""}
    }

}
