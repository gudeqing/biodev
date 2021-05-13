version 1.0

task rsem_quant{
    input {
        String? other_parameters
        Int threads = 8
        String strandness = "none"
        Boolean estimate_rspd = true
        Boolean append_names = true
        String? aligner
        Boolean input_is_bam = false
        Boolean paired_end = true
        File? bam
        Array[File]? read1
        Array[File]? read2
        File index
        String sample_name = "sample_name"
        # for runtime
        String docker = "rsem:1.3.3"
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
        ~{"--strandedness " + strandness} \
        ~{if estimate_rspd then "--estimate-rspd " else ""} \
        ~{if append_names then "--append-names " else ""} \
        ~{aligner} \
        ~{if input_is_bam then "--alignments " else ""} \
        ~{if paired_end then "--paired-end " else ""} \
        ~{bam} \
        ~{sep="," read1} \
        ~{sep="," read2} \
        ~{index} \
        ~{sample_name} 
    >>>

    output {
        File genes = "${sample_name}.genes.results"
        File isoforms = "${sample_name}.isoforms.results"
        File alleles = "${sample_name}.alleles.results"
    }

    runtime {
        docker: docker
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "rsem_quant"
        image: "rsem:1.3.3"
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
        index: {desc: "reference index files", level: "required", type: "infile", range: "", default: ""}
        sample_name: {desc: "prefix for output file name", level: "required", type: "str", range: "", default: "sample_name"}
    }

}
