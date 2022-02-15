# ChIP-seq-analysis

## 1. Trimming/QC des reads.

Je te conseille d'utiliser *fastp* pour ça.

La ligne de commande:
`fastp -A -h prefix_qc.html -i prefix_reads1.fastq.gz -I prefix_reads2.fastq.gz`

Link: https://github.com/OpenGene/fastp

## 2. Alignement sur génome de référence
