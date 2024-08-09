#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/reads/*_R{1,2}_001.fastq.gz"
params.transcriptome_file = "$projectDir/data/transcriptome/hsapiens_hg38.fa.gz"
params.multiqc = "$projectDir/multiqc"
params.quant_reads = "$projectDir/results/*/quant.sf"
params.outdir = "results"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    input:
    path transcriptome

    output:
    path 'hg38_index'

    script:
    """
    salmon index --threads 6 -t $transcriptome -i hg38_index
    """
}

process QUANTIFICATION {
    tag "Salmon on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    salmon quant --threads 6 --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}

process FASTQC {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

/* This process exists because the output filename of salmon is exactly the same
 * So I needed to access one file at a time to prevent Nextflow from crashing
 * then I wrote the estimated counts in a properly named file
*/
process EXTRACT_COUNTS{
    publishDir params.outdir, mode:'copy'

    input:
    tuple path(quant_file), val(parent_folder)

    output:
    path "${parent_folder}.tsv"

    script:
    """
    Rscript ${projectDir}/bin/read_counts.R ${quant_file} ${parent_folder}
    """


}

process MERGE_COUNTS{
    publishDir params.outdir, mode:'copy'

    input:
    path list_count_tables

    output:
    path "Merged_SalmonOutput.tsv"

    script:
    """
    Rscript ${projectDir}/bin/merge_counts.R ${list_count_tables}
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    Channel
        .fromPath(params.quant_reads, checkIfExists: true)
        .map { quant_file ->
            def parent_folder = quant_file.parent.name  // Extract the parent folder name
            return tuple(quant_file, parent_folder)
        }
        .set { quant_reads_ch }
    Channel
        .fromPath("results/*.tsv")
        .collect()
        .set { quant_list_ch}    
    

    index_ch = INDEX(params.transcriptome_file)
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(quant_ch.mix(fastqc_ch).collect())
    EXTRACT_COUNTS(quant_reads_ch)
    MERGE_COUNTS(quant_list_ch)


}
