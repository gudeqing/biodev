version development

import "tasks/fastp.wdl" as fastp
import "tasks/star_alignment.wdl" as star_align
import "tasks/salmon_quant.wdl" as salmon_quant
import "tasks/rsem_quant.wdl" as rsem_quant
import "tasks/star_fusion.wdl" as star_fusion
import "tasks/markDuplicates.wdl" as markdup
import "tasks/geneBodyCoverage.wdl" as geneBodyCoverage
import "tasks/read_distribution.wdl" as read_distribution
import "tasks/read_duplication.wdl" as read_duplication
import "tasks/fpkm.wdl" as fpkm
import "tasks/CollectRnaSeqMetrics.wdl" as CollectRnaSeqMetrics
import "tasks/rnaseqc.wdl" as rnaseqc


workflow rnaseq_pipeline {
    input {
        String sample_id
        File read1
        File read2

        # for skip steps
        Boolean skip_fastp = false
        Boolean skip_rsem_quant = false
        Boolean skip_salmon_quant = true
        Boolean skip_fusion = false
    }

    if (! skip_fastp) {
        call fastp.fastp as fastp{
            input:
                read1 = read1,
                read2 = read2,
                sample_name = sample_id
        }
    }

    call star_align.star_alignment as align {
        input:
            sample = sample_id,
            read1 = [select_first([fastp.out_read1_file, read1])],
            read2 = [select_first([fastp.out_read2_file, read2])]
        }

    if (!skip_salmon_quant) {
        call salmon_quant.salmon_quant as salmon_quant {
            input:
                transcript_bam = [align.transcript_bam],
                outdir = sample_id
        }
    }

    if (!skip_rsem_quant) {
        call rsem_quant.rsem_quant as rsem_quant {
            input:
                sample_name = sample_id,
                input_is_bam = true,
                aligner = "--star",
                bam = align.transcript_bam
        }
    }

    if (! skip_fusion) {
        call star_fusion.star_fusion as fusion {
            input:
                sample = sample_id,
                chimeric_junction = align.chimeric_out,
        }
    }

    call markdup.markDuplicates as markdup {
        input:
            input_bam = align.bam,
            sample_id = sample_id
    }

    call rnaseqc.rnaseqc as rnaseqc {
        input:
            sample_id = sample_id,
            bam = markdup.bam_file,
    }

    call CollectRnaSeqMetrics.CollectRnaSeqMetrics as CollectRnaSeqMetrics {
        input:
            sample_id = sample_id,
            input_bam = markdup.bam_file
    }

    call read_distribution.read_distribution as read_distribution {
        input:
            sample_id = sample_id,
            bam = align.bam
    }

    call read_duplication.read_duplication as read_duplication {
        input:
            sample_id = sample_id,
            bam = align.bam
    }

    call geneBodyCoverage.geneBodyCoverage as geneBodyCoverage {
        input:
            sample_id = sample_id,
            bam = [align.bam]
    }

    output {
        File? fastq_html_qc = fastp.html_report_file
        File? fusion_file = fusion.fusion_predictions_abridged
        File? genome_bam = align.bam
    }

}