version 1.0

task STAR_Alignment{
    input {
        String? other_parameters
        Int runThreadN = 8
        # https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
        String genomeDir
        Array[File] read1
        Array[File]? read2
        String sample
        String platform = "Illimina"
        String outSAMtype = "BAM SortedByCoordinate"
        String outSAMunmapped = "Within"
        String readFilesCommand = "zcat"
        String twopassMode = "Basic"
        Int outFilterMultimapNmax = 20
        Int alignSJoverhangMin = 8
        Int alignSJDBoverhangMin = 10
        Int chimSegmentMin = 12
        Float outFilterMismatchNoverLmax = 0.3
        String outFilterType = "Normal"
        String outSAMstrandField = "intronMotif"
        String quantMode = "TranscriptomeSAM"
        String outSAMattrRGline = "ID:~{sample} SM:~{sample} PL:~{platform}"
        Int limitBAMsortRAM = 35000000000
        Int limitIObufferSize = 150000000
        Int outSAMattrIHstart = 0
        Int alignMatesGapMax = 500000
        Int alignIntronMax = 500000
        String alignSJstitchMismatchNmax = "5 -1 5 5"
        Int chimJunctionOverhangMin = 8
        Int chimMultimapScoreRange = 3
        Int chimMultimapNmax = 20
        Int chimNonchimScoreDropMin = 10
        Int chimOutJunctionFormat = 1
        Int peOverlapNbasesMin = 12
        Float peOverlapMMp = 0.1
        String chimOutType = "WithinBAM"
        String alignInsertionFlush = "Right"
        Int chimScoreJunctionNonGTAG = -4
        Float alignSplicedMateMapLminOverLmate = 0
        Int alignSplicedMateMapLmin = 30
        # for runtime
        String docker = "STAR2.7.8a:latest"
        String memory = "32 GiB"
        Int cpu = 1
        String disks = "50 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e 
        STAR \
        ~{other_parameters} \
        ~{"--runThreadN " + runThreadN} \
        ~{"--genomeDir " + genomeDir} \
        --readFilesIn ~{sep="," read1}  ~{sep="," read2}  \
        ~{"--outFileNamePrefix " + sample + "."} \
        --outSAMtype ~{outSAMtype} \
        ~{"--outSAMunmapped " + outSAMunmapped} \
        ~{"--readFilesCommand " + readFilesCommand} \
        ~{"--twopassMode " + twopassMode} \
        ~{"--outFilterMultimapNmax " + outFilterMultimapNmax} \
        ~{"--alignSJoverhangMin " + alignSJoverhangMin} \
        ~{"--alignSJDBoverhangMin " + alignSJDBoverhangMin} \
        ~{"--chimSegmentMin " + chimSegmentMin} \
        ~{"--outFilterMismatchNoverLmax " + outFilterMismatchNoverLmax} \
        ~{"--outFilterType " + outFilterType} \
        ~{"--outSAMstrandField " + outSAMstrandField} \
        ~{"--quantMode " + quantMode} \
        ~{"--outSAMattrRGline " + outSAMattrRGline} \
        ~{"--limitBAMsortRAM " + limitBAMsortRAM} \
        ~{"--limitIObufferSize " + limitIObufferSize} \
        ~{"--outSAMattrIHstart " + outSAMattrIHstart} \
        ~{"--alignMatesGapMax " + alignMatesGapMax} \
        ~{"--alignIntronMax " + alignIntronMax} \
        ~{"--alignSJstitchMismatchNmax " + alignSJstitchMismatchNmax} \
        ~{"--chimJunctionOverhangMin " + chimJunctionOverhangMin} \
        ~{"--chimMultimapScoreRange " + chimMultimapScoreRange} \
        ~{"--chimMultimapNmax " + chimMultimapNmax} \
        ~{"--chimNonchimScoreDropMin " + chimNonchimScoreDropMin} \
        ~{"--chimOutJunctionFormat " + chimOutJunctionFormat} \
        ~{"--peOverlapNbasesMin " + peOverlapNbasesMin} \
        ~{"--peOverlapMMp " + peOverlapMMp} \
        ~{"--chimOutType " + chimOutType} \
        ~{"--alignInsertionFlush " + alignInsertionFlush} \
        ~{"--chimScoreJunctionNonGTAG " + chimScoreJunctionNonGTAG} \
        ~{"--alignSplicedMateMapLminOverLmate " + alignSplicedMateMapLminOverLmate} \
        ~{"--alignSplicedMateMapLmin " + alignSplicedMateMapLmin} 
    >>>

    output {
        File bam = glob("*.bam")
        File align_log = glob("*Log.final.out")[0]
        File? chimeric_out = glob("*Chimeric.out.junction")
        File? sj = glob("SJ.out.tab")
    }

    runtime {
        docker: docker
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "STAR_Alignment"
        image: "docker build -t STAR2.7.8a -f ./Dockerfile"
        desc: "STAR (Spliced Transcript Alignment to a Reference) is an ultrafast RNA-seq aligner, is capable of mapping full length RNA sequences and detecting de novo canonical junctions, non-canonical splices, and chimeric (fusion) transcripts. It aligns short and long RNA-seq reads to a reference genome using sequential maximum mappable seed search in uncompressed suffix arrays followed by seed clustering and stitching procedure."
        logo: "star.png"
        version: "2.7.8a"
        source: "https://github.com/alexdobin/STAR"
        basecmd: "STAR"
    }

    parameter_meta {
        other_parameters: {desc: "other arguments, you could set any other argument with a string such as '-i x -j y'", level: "optional", type: "str", value_candidates: ""}
        runThreadN: {desc: "Number of threads to use", level: "required", type: "int", value_candidates: ""}
        genomeDir: {desc: "reference index directory", level: "required", type: "indir", value_candidates: ""}
        readFilesIn: {desc: "fastq files, separate by white space, such as 'R1.fq R2.fq'", level: "required", type: "infile", value_candidates: ""}
        sample: {desc: "prefix for outfile name", level: "required", type: "str", value_candidates: ""}
        outSAMtype: {desc: "format specification for output bam/sam", level: "required", type: "str", value_candidates: ""}
        outSAMunmapped: {desc: "format specification for output bam/sam", level: "required", type: "str", value_candidates: ""}
        readFilesCommand: {desc: "Un-compression Command for each fastq，could be 'zcat' or 'gunzip -c' or '-'", level: "required", type: "str", value_candidates: ""}
        twopassMode: {desc: "set as 'Basic' for basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly", level: "required", type: "str", value_candidates: ""}
        outFilterMultimapNmax: {desc: "max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped", level: "required", type: "int", value_candidates: ""}
        alignSJoverhangMin: {desc: "minimum overhang for unannotated junctions", level: "required", type: "int", value_candidates: ""}
        alignSJDBoverhangMin: {desc: "minimum overhang for annotated junctions", level: "required", type: "int", value_candidates: ""}
        chimSegmentMin: {desc: "parameter controls the minimum mapped length of the two segments that is allowed. For example, if you have 2x75 reads and used --chimSegmentMin 20, a chimeric alignment with 130b on one chromosome and 20b on the other will be output, while 135 + 15 won’t be", level: "required", type: "int", value_candidates: ""}
        outFilterMismatchNoverLmax: {desc: "alignment will be output only if its ratio of mismatches to mapped length is less than or equal to this value.", level: "required", type: "float", value_candidates: ""}
        outFilterType: {desc: "type of filtering", level: "required", type: "str", value_candidates: ""}
        outSAMstrandField: {desc: "include for potential use with StringTie for assembly", level: "required", type: "str", value_candidates: ""}
        quantMode: {desc: "output transcriptome bam for expression quantification", level: "required", type: "str", value_candidates: ""}
        outSAMattrRGline: {desc: "SAM/BAM read group line. The first word contains the read group identifier and must start with 'ID:'", level: "required", type: "str", value_candidates: ""}
        limitBAMsortRAM: {desc: "int>=0: maximum available RAM (bytes) for sorting BAM. If =0, it will be set to the genome index size. 0 value can only be used with –genomeLoad NoSharedMemory option.", level: "required", type: "int", value_candidates: ""}
        limitIObufferSize: {desc: "int>0: max available buffers size (bytes) for input/output, per thread", level: "required", type: "int", value_candidates: ""}
        outSAMattrIHstart: {desc: "start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie.", level: "required", type: "int", value_candidates: ""}
        alignMatesGapMax: {desc: "maximum gap between two mates", level: "required", type: "int", value_candidates: ""}
        alignIntronMax: {desc: "maximum intron size", level: "required", type: "int", value_candidates: ""}
        alignSJstitchMismatchNmax: {desc: "maximum number of mismatches for stitching of the splice junctions.(1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif.", level: "required", type: "str", value_candidates: ""}
        chimJunctionOverhangMin: {desc: "minimum overhang for a chimeric junction", level: "required", type: "int", value_candidates: ""}
        chimMultimapScoreRange: {desc: "the score range for multi-mapping chimeras below the best chimeric score. Only works with –chimMultimapNmax > 1", level: "required", type: "int", value_candidates: ""}
        chimMultimapNmax: {desc: "maximum number of chimeric multi-alignments", level: "required", type: "int", value_candidates: ""}
        chimNonchimScoreDropMin: {desc: "to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be smaller than this value", level: "required", type: "int", value_candidates: ""}
        chimOutJunctionFormat: {desc: "if 1: add comment lines at the end of the file: command line and Nreads: total, unique, multi", level: "required", type: "int", value_candidates: ""}
        peOverlapNbasesMin: {desc: "minimum number of overlap bases to trigger mates merging and realignment", level: "required", type: "int", value_candidates: ""}
        peOverlapMMp: {desc: "maximum proportion of mismatched bases in the overlap area", level: "required", type: "float", value_candidates: ""}
        chimOutType: {desc: "type of chimeric output", level: "required", type: "str", value_candidates: ""}
        alignInsertionFlush: {desc: "how to ush ambiguous insertion positions", level: "required", type: "str", value_candidates: ""}
        chimScoreJunctionNonGTAG: {desc: "penalty for a non-GT/AG chimeric junction", level: "required", type: "int", value_candidates: ""}
        alignSplicedMateMapLminOverLmate: {desc: "alignSplicedMateMapLmin normalized to mate length", level: "required", type: "float", value_candidates: ""}
        alignSplicedMateMapLmin: {desc: "minimum mapped length for a read mate that is spliced", level: "required", type: "int", value_candidates: ""}
    }

}
