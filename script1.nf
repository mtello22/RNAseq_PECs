params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"

log.info """
	RNAseq - NF Pipeline
	==============================================
	reads:		$params.reads
	transcriptome:	$params.transcriptome_file
	MultiQC:	$params.multiqc
	outdir:		$params.outdir
"""
.stripIndent(true)

