version development
# 如果原始输入数据是cleandata，则可跳过fastp
# 如果不关心转录本定量结果，则可跳过rsem_quant步骤
# 如果不需要fusion分析，则可跳过fusion步骤
# 如果不关心circRNA的分析，则可以跳过circRNA步骤

workflow rnaseq_pipeline {
    input {
        # for skip steps
        Boolean skip_fastp = false
        Boolean skip_rsem_quant = false
        Boolean skip_fusion = false
        Boolean skip_circRNA = false
    }

    call getFastqInfo{}

    scatter (each in keys(getFastqInfo.fastq_info)) {
        String sample_id = each
        File read1_file = getFastqInfo.fastq_info[each][0][0]
        File read2_file = getFastqInfo.fastq_info[each][1][0]

        if (! skip_fastp) {
            call fastp {
                input:
                    read1 = read1_file,
                    read2 = read2_file,
                    sample_name = sample_id
            }
        }

        call star_alignment as align {
            input:
                sample = sample_id,
                read1 = [select_first([fastp.out_read1_file, read1_file])],
                read2 = [select_first([fastp.out_read2_file, read2_file])]
            }

        if (!skip_rsem_quant) {
            call rsem_quant {
                input:
                    sample_name = sample_id,
                    input_is_bam = true,
                    bam = align.transcript_bam
            }
        }

        if (! skip_fusion) {
            call star_fusion as fusion {
                input:
                    sample = sample_id,
                    chimeric_junction = align.chimeric_out,
                    left_fq = [select_first([fastp.out_read1_file, read1_file])],
                    right_fq = [select_first([fastp.out_read2_file, read2_file])]
            }
        }

        call markDuplicates as markdup {
            input:
                input_bam = align.bam,
                sample_id = sample_id
        }

        call rnaseqc {
            input:
                sample_id = sample_id,
                bam = markdup.bam_file,
        }

        call CollectRnaSeqMetrics {
            input:
                sample_id = sample_id,
                input_bam = markdup.bam_file
        }

        call read_distribution {
            input:
                sample_id = sample_id,
                bam = align.bam
        }

#        call geneBodyCoverage {
#            input:
#                sample_id = sample_id,
#                bam = [align.bam],
#                bam_bai = [align.bam_bai]
#        }

        if (! skip_circRNA) {
            call CIRCexplorer2 {
                input:
                    sample = sample_id,
                    chimeric_junction = align.chimeric_out
            }
        }
    }

    output {
        Array[File] align_log = align.align_log
        Array[File] chimeric_out = align.chimeric_out
        Array[File?] fastq_html_qc = fastp.html_report_file
        Array[File?] fastq_json_qc = fastp.json_report_file
        Array[File?] fusion_file = fusion.fusion_predictions_abridged
        Array[File?] genome_mdbam = markdup.bam_file
        Array[File?] genome_mdbam_bai = markdup.bam_index
        Array[File?] rnaseqc_gene_tpm = rnaseqc.gene_tpm
        Array[File?] rnaseqc_gene_counts = rnaseqc.gene_counts
        Array[File?] rnaseqc_metrics = rnaseqc.metrics
        Array[File?] rsem_gene_quant = rsem_quant.genes
        Array[File?] rsem_trans_quant = rsem_quant.isoforms
        Array[File?] read_distr = read_distribution.read_distr
        Array[File?] picard_metrics = CollectRnaSeqMetrics.rnaseq_metrics
        Array[File?] circRNA = CIRCexplorer2.circRNA
    }
}

# following are tasks
task getFastqInfo{
    input {
        Array[Directory]? fastq_dirs
        Array[File]? fastq_files
        String r1_name = '(.*).read1.fastq.gz'
        String r2_name = '(.*).read2.fastq.gz'
        String docker = 'gudeqing/getfastqinfo:1.0'
    }

    command <<<
        set -e
        python /get_fastq_info.py \
            ~{if defined(fastq_dirs) then "-fastq_dirs " else ""}~{sep=" " fastq_dirs} \
            ~{if defined(fastq_files) then "-fastq_files " else ""}~{sep=" " fastq_files} \
            -r1_name '~{r1_name}' \
            -r2_name '~{r2_name}' \
            -out fastq.info.json
    >>>

    output {
        Map[String, Array[Array[File]]] fastq_info = read_json("fastq.info.json")
    }

    runtime {
        docker: docker
    }

    parameter_meta {
        fastq_dirs: {desc: "directory list, target fastq files should be in these directories. All target files in 'fastq_files' or 'fastq_dirs' will be used", level: "optional", type: "indir", range: "", default: ""}
        fastq_files: {desc: "target fastq file list. 'fastq_files' or 'fastq_dirs' must be provided.", level: "optional", type: "infile", range: "", default: ""}
        r1_name: {desc: "python regExp that describes the full name of read1 fastq file name. It requires at least one pair of small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'", level: "required", type: "str", range: "", default: ""}
        r2_name: {desc: "python regExp that describes the full name of read2 fastq file name. It requires at least one pair of small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'", level: "required", type: "str", range: "", default: ""}
    }
}

task fastp{
    input {
        String? other_parameters
        Int threads = 4
        File read1
        File? read2
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
        fastp \
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
        docker: "registry-xdp-v3-pre-yifang.basebit.me/basebitai/fastp:0.21.0"
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "fastp"
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

task star_alignment{
    input {
        String? other_parameters
        Int runThreadN = 6
        # https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/
        # Directory indexDir
        Array[File] indexFiles
        Array[File] read1
        # 下面得如果不默认为[],则cromwell会报错
        Array[File] read2 = []
        String sample
        String platform = "Illumina"
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
        String limitBAMsortRAM =  "35000000000"
        String limitIObufferSize = "150000000"
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
        String chimOutType = "Junctions WithinBAM"
        String alignInsertionFlush = "Right"
        Int chimScoreJunctionNonGTAG = -4
        Float alignSplicedMateMapLminOverLmate = 0
        Int alignSplicedMateMapLmin = 30
        String quantTranscriptomeBan = "IndelSoftclipSingleend"
        # for runtime
        String docker = "trinityctat/starfusion:1.10.0"
        String memory = "35 GiB"
        Int cpu = 6
        String disks = "50 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e
        STAR \
        ~{other_parameters} \
        ~{"--runThreadN " + runThreadN} \
        --genomeDir ~{sub(indexFiles[0], basename(indexFiles[0]), "")} \
        --readFilesIn ~{sep="," read1}  ~{sep="," read2} \
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
        ~{"--alignSplicedMateMapLmin " + alignSplicedMateMapLmin} \
        ~{"--quantTranscriptomeBan " + quantTranscriptomeBan}
        samtools index ~{sample}.Aligned.sortedByCoord.out.bam
        samtools flagstat ~{sample}.Aligned.sortedByCoord.out.bam > ~{sample}.align.flagstat.txt
        samtools idxstats ~{sample}.Aligned.sortedByCoord.out.bam > ~{sample}.align.idxstats.txt
    >>>

    output {
        File bam = "~{sample}.Aligned.sortedByCoord.out.bam"
        File bam_bai = "~{sample}.Aligned.sortedByCoord.out.bam.bai"
        File transcript_bam = "~{sample}.Aligned.toTranscriptome.out.bam"
        File align_log = "~{sample}.Log.final.out"
        File chimeric_out = "~{sample}.Chimeric.out.junction"
        File sj = "~{sample}.SJ.out.tab"
    }

    runtime {
        docker: "registry-xdp-v3-pre-yifang.basebit.me/basebitai/starfusion:1.10.0"
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "STAR_Alignment"
        desc: "STAR (Spliced Transcript Alignment to a Reference) is an ultrafast RNA-seq aligner, is capable of mapping full length RNA sequences and detecting de novo canonical junctions, non-canonical splices, and chimeric (fusion) transcripts. It aligns short and long RNA-seq reads to a reference genome using sequential maximum mappable seed search in uncompressed suffix arrays followed by seed clustering and stitching procedure."
        logo: "star.png"
        version: "2.7.8a"
        source: "https://github.com/alexdobin/STAR"
        basecmd: "STAR"
    }

    parameter_meta {
        other_parameters: {desc: "other arguments, you could set any other argument with a string such as '-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        runThreadN: {desc: "Number of threads to use", level: "required", type: "int", range: "", default: "8"}
        indexFiles: {desc: "reference index files", level: "required", type: "indir", range: "", default: ""}
        read1: {desc: "read1 fastq files", level: "required", type: "infile", range: "", default: ""}
        read2: {desc: "read2 fastq files", level: "required", type: "infile", range: "", default: ""}
        sample: {desc: "prefix for outfile name", level: "required", type: "str", range: "", default: ""}
        outSAMtype: {desc: "format specification for output bam/sam", level: "required", type: "str", range: "", default: "BAM SortedByCoordinate"}
        outSAMunmapped: {desc: "format specification for output bam/sam", level: "required", type: "str", range: "", default: "Within"}
        readFilesCommand: {desc: "Un-compression Command for each fastq，could be 'zcat' or 'gunzip -c' or '-'", level: "required", type: "str", range: "", default: "zcat"}
        twopassMode: {desc: "set as 'Basic' for basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly", level: "required", type: "str", range: "", default: "Basic"}
        outFilterMultimapNmax: {desc: "max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped", level: "required", type: "int", range: "", default: "20"}
        alignSJoverhangMin: {desc: "minimum overhang for unannotated junctions", level: "required", type: "int", range: "", default: "8"}
        alignSJDBoverhangMin: {desc: "minimum overhang for annotated junctions", level: "required", type: "int", range: "", default: "10"}
        chimSegmentMin: {desc: "parameter controls the minimum mapped length of the two segments that is allowed. For example, if you have 2x75 reads and used --chimSegmentMin 20, a chimeric alignment with 130b on one chromosome and 20b on the other will be output, while 135 + 15 won’t be", level: "required", type: "int", range: "", default: "12"}
        outFilterMismatchNoverLmax: {desc: "alignment will be output only if its ratio of mismatches to mapped length is less than or equal to this value.", level: "required", type: "float", range: "", default: "0.3"}
        outFilterType: {desc: "type of filtering", level: "required", type: "str", range: "", default: "Normal"}
        outSAMstrandField: {desc: "include for potential use with StringTie for assembly", level: "required", type: "str", range: "", default: "intronMotif"}
        quantMode: {desc: "output transcriptome bam for expression quantification", level: "required", type: "str", range: "", default: "TranscriptomeSAM"}
        outSAMattrRGline: {desc: "SAM/BAM read group line. The first word contains the read group identifier and must start with 'ID:'", level: "required", type: "str", range: "", default: "ID:sample SM:sample PL:Illumina"}
        limitBAMsortRAM: {desc: "int>=0: maximum available RAM (bytes) for sorting BAM. If =0, it will be set to the genome index size. 0 value can only be used with –genomeLoad NoSharedMemory option.", level: "required", type: "int", range: "", default: "25000000000"}
        limitIObufferSize: {desc: "int>0: max available buffers size (bytes) for input/output, per thread", level: "required", type: "int", range: "", default: "150000000"}
        outSAMattrIHstart: {desc: "start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie.", level: "required", type: "int", range: "", default: "0"}
        alignMatesGapMax: {desc: "maximum gap between two mates", level: "required", type: "int", range: "", default: "500000"}
        alignIntronMax: {desc: "maximum intron size", level: "required", type: "int", range: "", default: "500000"}
        alignSJstitchMismatchNmax: {desc: "maximum number of mismatches for stitching of the splice junctions.(1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif.", level: "required", type: "str", range: "", default: "5 -1 5 5"}
        chimJunctionOverhangMin: {desc: "minimum overhang for a chimeric junction", level: "required", type: "int", range: "", default: "8"}
        chimMultimapScoreRange: {desc: "the score range for multi-mapping chimeras below the best chimeric score. Only works with –chimMultimapNmax > 1", level: "required", type: "int", range: "", default: "3"}
        chimMultimapNmax: {desc: "maximum number of chimeric multi-alignments", level: "required", type: "int", range: "", default: "20"}
        chimNonchimScoreDropMin: {desc: "to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be smaller than this value", level: "required", type: "int", range: "", default: "10"}
        chimOutJunctionFormat: {desc: "if 1: add comment lines at the end of the file: command line and Nreads: total, unique, multi", level: "required", type: "int", range: "", default: "1"}
        peOverlapNbasesMin: {desc: "minimum number of overlap bases to trigger mates merging and realignment", level: "required", type: "int", range: "", default: "12"}
        peOverlapMMp: {desc: "maximum proportion of mismatched bases in the overlap area", level: "required", type: "float", range: "", default: "0.1"}
        chimOutType: {desc: "type of chimeric output", level: "required", type: "str", range: "", default: "WithinBAM"}
        alignInsertionFlush: {desc: "how to flush ambiguous insertion positions", level: "required", type: "str", range: "", default: "Right"}
        chimScoreJunctionNonGTAG: {desc: "penalty for a non-GT/AG chimeric junction", level: "required", type: "int", range: "", default: "-4"}
        alignSplicedMateMapLminOverLmate: {desc: "alignSplicedMateMapLmin normalized to mate length", level: "required", type: "float", range: "", default: "0"}
        alignSplicedMateMapLmin: {desc: "minimum mapped length for a read mate that is spliced", level: "required", type: "int", range: "", default: "30"}
        quantTranscriptomeBan: {desc: "prohibit various alignment type", level: "required", type: "str", range: "IndelSoftclipSingleend,Singleend", default: "IndelSoftclipSingleend"}
    }

}

task rnaseqc{
    input {
        String? other_parameters
        File collapsed_gtf
        File bam
        String sample_id
        File? bed
        String? strand
        # for runtime
        String docker = "gudeqing/rnaseqc:2.4.2"
        String memory = "6 GiB"
        Int cpu = 2
        String disks = "6 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e
        rnaseqc \
        ~{other_parameters} \
        ~{collapsed_gtf} \
        ~{bam} \
        ./ \
        ~{"--sample " + sample_id} \
        ~{"--bed " + bed} \
        ~{"--stranded " + strand}
    >>>

    output {
        File gene_tpm = "${sample_id}.gene_tpm.gct"
        File gene_counts =  "${sample_id}.gene_reads.gct"
        File exon_counts =  "${sample_id}.exon_reads.gct"
        File metrics = "${sample_id}.metrics.tsv"
        File? insertsize_distr =  "${sample_id}.fragmentSizes.txt"
    }

    runtime {
        docker: "registry-xdp-v3-pre-yifang.basebit.me/basebitai/rnaseqc:2.4.2"
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "rnaseqc"
        desc: "Fast efficient RNA-Seq metrics for quality control and process optimization"
        version: "2.4.2"
        source: "https://github.com/getzlab/rnaseqc"
        basecmd: "rnaseqc"
    }

    parameter_meta {
        other_parameters: {desc: "other arguments, you could set any other argument with a string such as '-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        collapsed_gtf: {desc: "The input GTF file containing features to check the bam against", level: "required", type: "infile", range: "", default: ""}
        bam: {desc: "The input SAM/BAM file containing reads", level: "required", type: "infile", range: "", default: ""}
        sample_id: {desc: "prefix for output file name", level: "optional", type: "str", range: "", default: ""}
        bed: {desc: "Optional input BED file containing non-overlapping exons used for fragment size calculations", level: "optional", type: "infile", range: "", default: ""}
        strand: {desc: "Use strand-specific metrics. Only features on the same strand of a read will be considered.", level: "optional", type: "str", range: "RF,FR", default: ""}
    }

}

task rsem_quant{
    input {
        String? other_parameters
        Int threads = 6
        String strandness = "none"
        Boolean estimate_rspd = true
        Boolean append_names = true
        String? aligner
        Boolean input_is_bam = false
        Boolean paired_end = true
        File? bam
        Array[File]? read1
        Array[File]? read2
#        Directory indexDir
        Array[File] indexFiles
        String index_prefix = 'rsem'
        String sample_name = "sample_name"
        # for runtime
        String docker = "gudeqing/rsem:1.3.3"
        String memory = "10 GiB"
        Int cpu = 2
        String disks = "10 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e
        rsem-calculate-expression \
        ~{other_parameters} \
        ~{"-p " + threads} \
        --no-bam-output \
        ~{"--strandedness " + strandness} \
        ~{if estimate_rspd then "--estimate-rspd " else ""} \
        ~{if append_names then "--append-names " else ""} \
        ~{aligner} \
        ~{if input_is_bam then "--alignments " else ""} \
        ~{if paired_end then "--paired-end " else ""} \
        ~{bam} \
        ~{sep="," read1} \
        ~{sep="," read2} \
        ~{sub(indexFiles[0], basename(indexFiles[0]), "")}/~{index_prefix} \
        ~{sample_name}
    >>>

    output {
        File genes = "${sample_name}.genes.results"
        File isoforms = "${sample_name}.isoforms.results"
        File? alleles = "${sample_name}.alleles.results"
    }

    runtime {
        docker: "registry-xdp-v3-pre-yifang.basebit.me/basebitai/rsem:1.3.3"
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "rsem_quant"
        desc: "RSEM is a software package for estimating gene and isoform expression levels from RNA-Seq data. The RSEM package provides an user-friendly interface, supports threads for parallel computation of the EM algorithm, single-end and paired-end read data, quality scores, variable-length reads and RSPD estimation. In addition, it provides posterior mean and 95% credibility interval estimates for expression levels."
        logo: "rsem.png"
        version: "1.3.3"
        source: "https://github.com/deweylab/RSEM"
        basecmd: "rsem-calculate-expression"
    }

    parameter_meta {
        other_parameters: {desc: "other arguments, you could set any other argument with a string such as '-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        threads: {desc: "Number of threads to use", level: "required", type: "int", range: "", default: "8"}
        strandness: {desc: "This option defines the strandedness of the RNA-Seq reads. It recognizes three values: 'none', 'forward', and 'reverse'. 'none' refers to non-strand-specific protocols. 'forward' means all (upstream) reads are derived from the forward strand. 'reverse' means all (upstream) reads are derived from the reverse strand. If 'forward'/'reverse' is set, the '--norc'/'--nofw' Bowtie/Bowtie 2 option will also be enabled to avoid aligning reads to the opposite strand. For Illumina TruSeq Stranded protocols, please use 'reverse'. (Default: 'none')", level: "required", type: "str", range: "none,forward,reverse", default: "'none'"}
        estimate_rspd: {desc: "Set this option if you want to estimate the read start position distribution (RSPD) from data. Otherwise, RSEM will use a uniform RSPD. (Default: off)", level: "required", type: "bool", range: "yes,no", default: "yes"}
        append_names: {desc: "If gene_name/transcript_name is available, append it to the end of gene_id/transcript_id (separated by '_') in files 'sample_name.isoforms.results' and 'sample_name.genes.results'. (Default: off)", level: "required", type: "bool", range: "yes,no", default: "yes"}
        aligner: {desc: "aligner name, could be one of [--bowtie, --bowtie2, --star, --hisat2-hca]", level: "optional", type: "str", range: "--bowtie, --bowtie2, --star, --hisat2-hca", default: ""}
        input_is_bam: {desc: "tell if input is transcript alignment file", level: "required", type: "bool", range: "yes,no", default: "no"}
        paired_end: {desc: "tell if input file contains paired-end reads", level: "required", type: "bool", range: "yes,no", default: "yes"}
        bam: {desc: "input alignment file", level: "optional", type: "infile", range: "", default: ""}
        read1: {desc: "input read1 fastq file", level: "optional", type: "infile", range: "", default: ""}
        read2: {desc: "input read2 fastq file", level: "optional", type: "infile", range: "", default: ""}
        indexFiles: {desc: "reference index files", level: "required", type: "infile", range: "", default: ""}
        index_prefix: {desc: "prefix for index files, default to 'rsem'", level: "required", type: "str", range: "", default: "rsem"}
        sample_name: {desc: "prefix for output file name", level: "required", type: "str", range: "", default: "sample_name"}
    }

}

task star_fusion{
    input {
        String? other_parameters
        Int CPU = 5
        Array[File]? left_fq
        Array[File]? right_fq
        File? chimeric_junction
#        Directory genome_lib_dir
        Array[File] indexFiles
        String sample = "fusion"
        # 如果要使用FusionInspector和denovo_reconstruct功能，必须提供read1和read2
        String? FusionInspector = "inspect"
        Boolean examine_coding_effect = true
        Boolean denovo_reconstruct = true
        # for runtime
        String docker = "trinityctat/starfusion:1.10.0"
        String memory = "32 GiB"
        Int cpu = 5
        String disks = "50 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e
        mkdir -p ~{sample}
        STAR-Fusion \
        ~{other_parameters} \
        ~{"--CPU " + CPU} \
        ~{if defined(left_fq) then "--left_fq " else ""}~{sep=" " left_fq} \
        ~{if defined(right_fq) then "--right_fq " else ""}~{sep=" " right_fq} \
        ~{"-J " + chimeric_junction} \
        --genome_lib_dir ~{sub(indexFiles[0], basename(indexFiles[0]), "")} \
        ~{"--output_dir " + sample} \
        ~{"--FusionInspector " + FusionInspector} \
        ~{if examine_coding_effect then "--examine_coding_effect " else ""} \
        ~{if denovo_reconstruct then "--denovo_reconstruct " else ""}
    >>>

    output {
        File fusion_predictions = "${sample}/star-fusion.fusion_predictions.tsv"
        File fusion_predictions_abridged = "${sample}/star-fusion.fusion_predictions.abridged.tsv"
        File? bam = "${sample}/Aligned.out.bam"
        File? star_log_final = "${sample}/Log.final.out"
        File? junction = "${sample}/Chimeric.out.junction"
        File? sj = "${sample}/SJ.out.tab"
        Array[File]? extract_fusion_reads =  ["${sample}/star-fusion.fusion_evidence_reads_1.fq", "${sample}/star-fusion.fusion_evidence_reads_2.fq"]
        File? fusion_inspector_validate_fusions_abridged = "${sample}/FusionInspector-validate/finspector.FusionInspector.fusions.abridged.tsv"
        File? fusion_inspector_validate_web = "${sample}/FusionInspector-validate/finspector.fusion_inspector_web.html"
        File? fusion_inspector_inspect_web = "${sample}/FusionInspector-inspect/finspector.fusion_inspector_web.html"
        File? fusion_inspector_inspect_fusions_abridged = "${sample}/FusionInspector-inspect/finspector.FusionInspector.fusions.abridged.tsv"
    }

    runtime {
        docker: "registry-xdp-v3-pre-yifang.basebit.me/basebitai/starfusion:1.10.0"
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "STAR_Fusion"
        desc: "STAR-Fusion is a component of the Trinity Cancer Transcriptome Analysis Toolkit (CTAT). STAR-Fusion uses the STAR aligner to identify candidate fusion transcripts supported by Illumina reads. STAR-Fusion further processes the output generated by the STAR aligner to map junction reads and spanning reads to a reference annotation set."
        logo: "star-fusion.png"
        version: "v1.10.0"
        source: "https://github.com/STAR-Fusion/STAR-Fusion"
        basecmd: "STAR-Fusion"
    }

    parameter_meta {
        other_parameters: {desc: "other arguments, you could set any other argument with a string such as '-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        CPU: {desc: "Number of threads to use", level: "required", type: "int", range: "", default: "8"}
        left_fq: {desc: "read1 fastq files, separate by white space", level: "optional", type: "infile", range: "", default: ""}
        right_fq: {desc: "read1 fastq files, separate by white space", level: "optional", type: "infile", range: "", default: ""}
        chimeric_junction: {desc: "generated file called 'Chimeric.out.junction' by STAR alignment", level: "optional", type: "infile", range: "", default: ""}
        indexFiles: {desc: "all files in ctat_genome_lib_build_dir", level: "required", type: "indir", range: "", default: ""}
        sample: {desc: "output directory", level: "required", type: "str", range: "", default: "fusion"}
        FusionInspector: {desc: "FusionInspector that provides a more in-depth view of the evidence supporting the predicted fusions.", level: "required", type: "str", range: "inspect, validate", default: "inspect"}
        examine_coding_effect: {desc: "explore impact of fusions on coding sequences", level: "required", type: "bool", range: "yes, no", default: "yes"}
        denovo_reconstruct: {desc: "attempt to reconstruct fusion transcripts using Trinity de novo assembly (requires --FusionInspector)", level: "required", type: "bool", range: "yes, no", default: "yes"}
    }

}

task geneBodyCoverage{
    input {
        String? other_parameters
        Array[File] bam
        Array[File] bam_bai
        File refGene
        String sample_id = "sample_name"
        # for runtime
        String docker = "gudeqing/rseqc:4.0.0"
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
        File cov_txt = "~{sample_id}.geneBodyCoverage.txt"
        File cov_pdf = "~{sample_id}.geneBodyCoverage.curves.pdf"
    }

    runtime {
        docker: "registry-xdp-v3-pre-yifang.basebit.me/basebitai/rseqc:4.0.0"
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "geneBodyCoverage"
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

task markDuplicates{
    input {
        String? other_parameters
        File input_bam
        String sample_id
        String assume_sort_order = "coordinate"
        String optical_dup_pixel_distance = "2500"
        String tagging_policy = "DontTag"
        String create_index = "true"
        # for runtime
        String docker = "broadinstitute/picard:latest"
        String memory = "6 GiB"
        Int cpu = 2
        String disks = "6 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e
        java -jar /usr/picard/picard.jar MarkDuplicates \
        ~{other_parameters} \
        ~{"I=" + input_bam} \
        ~{"O=" + sample_id + ".markdup.bam"} \
        ~{"METRICS_FILE=" + sample_id + ".markdup.metrics.txt"} \
        ~{"ASSUME_SORT_ORDER=" + assume_sort_order} \
        ~{"OPTICAL_DUPLICATE_PIXEL_DISTANCE=" + optical_dup_pixel_distance} \
        ~{"TAGGING_POLICY=" + tagging_policy} \
        ~{"CREATE_INDEX=" + create_index}
    >>>

    output {
        File bam_file = "${sample_id}.markdup.bam"
        File bam_index = "${sample_id}.markdup.bai"
        File metrics = "${sample_id}.markdup.metrics.txt"
    }

    runtime {
        docker: "registry-xdp-v3-pre-yifang.basebit.me/basebitai/picard:latest"
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "markDuplicates"
        desc: "This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR"
        version: "2.23.3-1-g4ac48fc-SNAPSHOT"
        basecmd: "java -jar /usr/local/src/picard.jar MarkDuplicates"
    }

    parameter_meta {
        other_parameters: {desc: "other arguments, you could set any other argument with a string such as '-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        input_bam: {desc: "One or more input SAM or BAM files to analyze. Must be coordinate sorted", level: "required", type: "infile", range: "", default: ""}
        assume_sort_order: {desc: "Assume that the input file has this order even if the header says otherwise.", level: "required", type: "str", range: "unsorted, queryname, coordinate, duplicate,", default: "coordinate"}
        optical_dup_pixel_distance: {desc: "The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is moreappropriate. For other platforms and models, users should experiment to find what works best.  Default value: 100.", level: "required", type: "str", range: "unsorted, queryname, coordinate, duplicate,", default: "2500"}
        tagging_policy: {desc: "Determines how duplicate types are recorded in the DT optional attribute. Possible values: {DontTag, OpticalOnly, All}", level: "required", type: "str", range: "DontTag, OpticalOnly, All", default: "DontTag"}
        create_index: {desc: "Whether to create a BAM index when writing a coordinate-sorted BAM file.", level: "required", type: "str", range: "true,false", default: "true"}
    }

}

task read_distribution{
    input {
        String? other_parameters
        File bam
        File refGene
        String sample_id = "sample_name"
        # for runtime
        String docker = "gudeqing/rseqc:4.0.0"
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
        docker: "registry-xdp-v3-pre-yifang.basebit.me/basebitai/rseqc:4.0.0"
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "read_distribution"
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

task CollectRnaSeqMetrics{
    input {
        String? other_parameters
        File input_bam
        String sample_id
        String strand = "NONE"
        File ref_flat
        File? ribosomal_intervals
        # for runtime
        String docker = "broadinstitute/picard:latest"
        String memory = "6 GiB"
        Int cpu = 2
        String disks = "6 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e
        java -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
        ~{other_parameters} \
        ~{"I=" + input_bam} \
        ~{"O=" + sample_id + ".RnaSeqMetrics.txt"} \
        ~{"STRAND_SPECIFICITY=" + strand} \
        ~{"REF_FLAT=" + ref_flat} \
        ~{"RIBOSOMAL_INTERVALS=" + ribosomal_intervals} \
        ~{"CHART_OUTPUT=" + sample_id + ".coverage.pdf"}
    >>>

    output {
        File rnaseq_metrics = "~{sample_id}.RnaSeqMetrics.txt"
        File coverage_pdf = "~{sample_id}.coverage.pdf"
    }

    runtime {
        docker: "registry-xdp-v3-pre-yifang.basebit.me/basebitai/picard:latest"
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "CollectRnaSeqMetrics"
        desc: "Produces RNA alignment metrics for a SAM or BAM file"
        logo: "none"
        version: "2.23.3-1-g4ac48fc-SNAPSHOT"
        basecmd: "java -jar /usr/local/src/picard.jar CollectRnaSeqMetrics"
    }

    parameter_meta {
        other_parameters: {desc: "other arguments, you could set any other argument with a string such as '-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        input_bam: {desc: "One or more input SAM or BAM files to analyze. Must be coordinate sorted", level: "required", type: "infile", range: "", default: ""}
        strand: {desc: "For strand-specific library prep. For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.  Required. Possible values: {NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND}", level: "required", type: "str", range: "NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND", default: "NONE"}
        ref_flat: {desc: "Gene annotations in refFlat form.  Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat  Required.", level: "required", type: "infile", range: "", default: ""}
        ribosomal_intervals: {desc: "Location of rRNA sequences in genome, in interval_list format.  If not specified no bases will be identified as being ribosomal.", level: "optional", type: "infile", range: "", default: ""}
    }

}

task CIRCexplorer2{
    input {
        String sample
        File genome
        File genome_annot
        File chimeric_junction
         # for runtime
        String docker = "gudeqing/rseqc:4.0.0"
    }

    command <<<
        set -e
        grep -v '#' ~{chimeric_junction} | sed 1d > clean.chimeric_junction.txt
        CIRCexplorer2 parse -t STAR clean.chimeric_junction.txt -b ~{sample}.back_spliced_junction.bed
        CIRCexplorer2 annotate -r ~{genome_annot} \
        -g ~{genome} -b ~{sample}.back_spliced_junction.bed \
        -o ~{sample}.detected.circRNA.txt
    >>>

    output {
        File circRNA = "~{sample}.detected.circRNA.txt"
    }

    runtime {
        docker: "registry-xdp-v3-pre-yifang.basebit.me/basebitai/rseqc:4.0.0"
    }

    parameter_meta {
        sample: {desc: "sample name"}
        genome: {desc: "genome fasta file"}
        genome_annot: {desc: "genome annotation file, specially formated by CIRCexplore2"}
        chimeric_junction: {desc: "chimeric_junction file generated during star alignment"}
    }

}
