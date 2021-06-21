version development

task markDuplicates{
    input {
        String? other_parameters
        File input_bam
        String sample_id
        String assume_sort_order = "coordinate"
        String optical_dup_pixel_distance = "2500"
        String program_record_id = "null"
        String tagging_policy = "DontTag"
        String create_index = "true"
        # for runtime
        String docker = "trinityctat/starfusion:1.10.0"
        String memory = "6 GiB"
        Int cpu = 2
        String disks = "6 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e 
        java -jar /usr/picard/picard.jar MarkDuplicates \
        ~{other_parameters} \
        ~{"-I " + input_bam} \
        ~{"-O " + sample_id + ".markdup.bam"} \
        ~{"--METRICS_FILE " + sample_id + ".markdup.metrics.txt"} \
        ~{"--ASSUME_SORT_ORDER " + assume_sort_order} \
        ~{"--OPTICAL_DUPLICATE_PIXEL_DISTANCE " + optical_dup_pixel_distance} \
        ~{"--PROGRAM_RECORD_ID " + program_record_id} \
        ~{"--TAGGING_POLICY " + tagging_policy} \
        ~{"--CREATE_INDEX " + create_index} 
    >>>

    output {
        File bam_file = "${sample_id}.markdup.bam"
        File bam_index = "${sample_id}.markdup.bam.bai"
        File metrics = "${sample_id}.markdup.metrics.txt"
    }

    runtime {
        docker: docker
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "markDuplicates"
        docker: "broadinstitute/picard:latest"
        desc: "This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR"
        version: "2.23.3-1-g4ac48fc-SNAPSHOT"
        basecmd: "java -jar /usr/picard/picard.jar MarkDuplicates"
    }

    parameter_meta {
        other_parameters: {desc: "other arguments, you could set any other argument with a string such as '-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        input_bam: {desc: "One or more input SAM or BAM files to analyze. Must be coordinate sorted", level: "required", type: "infile", range: "", default: ""}
        assume_sort_order: {desc: "Assume that the input file has this order even if the header says otherwise.", level: "required", type: "str", range: "unsorted, queryname, coordinate, duplicate,", default: "coordinate"}
        optical_dup_pixel_distance: {desc: "The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is moreappropriate. For other platforms and models, users should experiment to find what works best.  Default value: 100.", level: "required", type: "str", range: "unsorted, queryname, coordinate, duplicate,", default: "2500"}
        program_record_id: {desc: "The program record ID for the @PG record(s) created by this program. Set to null to disable PG record creation.", level: "required", type: "str", range: "", default: "null"}
        tagging_policy: {desc: "Determines how duplicate types are recorded in the DT optional attribute. Possible values: {DontTag, OpticalOnly, All}", level: "required", type: "str", range: "DontTag, OpticalOnly, All", default: "DontTag"}
        create_index: {desc: "Whether to create a BAM index when writing a coordinate-sorted BAM file.", level: "required", type: "str", range: "true,false", default: "true"}
    }

}
