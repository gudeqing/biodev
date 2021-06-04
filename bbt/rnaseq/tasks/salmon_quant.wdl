version 1.0

task salmon_quant{
    input {
        String? other_parameters
        Int threads = 8
        Array[File]? indexFiles
        File? transcripts
        Array[File]+ transcript_bam
        Array[File]? read1
        Array[File]? read2
        Array[File]? single_end_reads
        File? geneMap
        String outdir = "salmon_quant"
        Boolean seqBias = true
        Boolean gcBias = true
        Boolean posBias = true
        # for runtime
        String docker = "combinelab/salmon:latest"
        String memory = "1 GiB"
        Int cpu = 1
        String disks = "1 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e 
        salmon quant \
        -l a \
        ~{other_parameters} \
        ~{"-p " + threads} \
        ~{if defined(indexFiles) then "-i " + sub(indexFiles[0], basename(indexFiles[0]), "") else ""} \
        ~{"-t " + transcripts} \
        ~{if defined(transcript_bam) then "-a " else ""}~{sep=" " transcript_bam} \
        ~{if defined(read1) then "-1 " else ""}~{sep=" " read1} \
        ~{if defined(read2) then "-2 " else ""}~{sep=" " read2} \
        ~{if defined(single_end_reads) then "-r " else ""}~{sep=" " single_end_reads} \
        ~{"-g " + geneMap} \
        ~{"-o " + outdir} \
        ~{if seqBias then "--seqBias " else ""} \
        ~{if gcBias then "--gcBias " else ""} \
        ~{if posBias then "--posBias " else ""} 
    >>>

    output {
        File outfile = outdir + "/quant.sf"
        File outfile2 = outdir + "/quant.genes.sf"
    }

    runtime {
        docker: docker
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "salmon"
        docker: "docker pull combinelab/salmon"
        source: "https://github.com/COMBINE-lab/salmon"
        desc: "Salmon is a wicked-fast program to produce a highly-accurate, transcript-level quantification estimates from RNA-seq data. Salmon achieves its accuracy and speed via a number of different innovations, including the use of selective-alignment (accurate but fast-to-compute proxies for traditional read alignments), and massively-parallel stochastic collapsed variational inference. The result is a versatile tool that fits nicely into many different pipelines. For example, you can choose to make use of our selective-alignment algorithm by providing Salmon with raw sequencing reads, or, if it is more convenient, you can provide Salmon with regular alignments (e.g. an unsorted BAM file with alignments to the transcriptome produced with your favorite aligner), and it will use the same wicked-fast, state-of-the-art inference algorithm to estimate transcript-level abundances for your experiment. For more detail please refer to https://github.com/COMBINE-lab/salmon"
        logo: "salmon.png"
        version: "1.4.0"
        basecmd: "salmon quant"
    }

    parameter_meta {
        other_parameters: {desc: "其他参数，你可以通过该参数输入一个或多个任何其他当前软件支持的参数，例如'-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        threads: {desc: "Number of threads to use during indexing or quantification", level: "required", type: "int", range: "", default: "8"}
        index_dir: {desc: "Existing directory containing transcripts indexing files for salmon quantification", level: "optional", type: "indir", range: "", default: ""}
        transcripts: {desc: "Transcript fasta file.", level: "optional", type: "infile", range: "", default: ""}
        transcript_bam: {desc: "input transcript based alignment (BAM) file(s)", level: "optional", type: "infile", range: "", default: ""}
        read1: {desc: "read1对应的fastq路径，支持用空格分隔多个文件路径的输入", level: "optional", type: "infile", range: "", default: ""}
        read2: {desc: "read2对应的fastq路径，支持用空格分隔多个文件路径的输入", level: "optional", type: "infile", range: "", default: ""}
        single_end_reads: {desc: "sing-end read对应的fastq路径, 支持输入多个文件, 空格分隔", level: "optional", type: "infile", range: "", default: ""}
        geneMap: {desc: "File containing a mapping of transcripts to genes. If this file is provided salmon will output both quant.sf and quant.genes.sf files, where the latter contains aggregated gene-level abundance estimates. The transcript to gene mapping should be provided as either a GTF file, or a ina simple tab-delimited format where each line contains the name of a transcript and the gene to which it belongs separated by a tab.", level: "optional", type: "infile", range: "", default: ""}
        outdir: {desc: "Output quantification directory.", level: "required", type: "str", range: "", default: "salmon_quant"}
        seqBias: {desc: "Bool argument，Perform sequence-specific bias correction.", level: "optional", type: "bool", range: "no, yes", default: "yes"}
        gcBias: {desc: "Bool argument. Perform fragment GC bias correction.", level: "optional", type: "bool", range: "no, yes", default: "yes"}
        posBias: {desc: "Bool argument. Perform positional bias correction.", level: "optional", type: "bool", range: "no, yes", default: "yes"}
    }

}
