import os
import glob

ERR,FRR = glob_wildcards("rawReads/{err}_{frr}.fastq")


rule all:
    input:
        expand("reducedReads/{err}_{frr}_reduced.fastq",  err=ERR, frr=FRR),
        #expand("featureCounts/{err}.sorted.bam", err=ERR)

rule reduce:
    input:
        rawReads1="rawReads/{err}_1.fastq",
        rawReads2="rawReads/{err}_2.fastq"
    output:
        "reducedReads/{err}_1_reduced.fastq",
        "reducedReads/{err}_2_reduced.fastq"
    threads:1
    params: 
        reduction_ratio=0.05,
    shell:
        """
        seqtk sample -s11 {input.rawReads1} {params.reduction_ratio} > {output[0]}
        seqtk sample -s11 {input.rawReads2} {params.reduction_ratio} > {output[1]}
        """
    
rule align:
    input:
        reducedRead1="reducedReads/{err}_1_reduced.fastq",
        reducedRead2="reducedReads/{err}_2_reduced.fastq"
    output:
        "alignedReads/{err}.sam"
        "alignedReads/{err}_summary.txt"
    threads: 1
    params:
        ref_gene = "/local/work/biocore/mol8008/RNA/Human_hg20/genome_tran"
    shell:
        """
        hisat2 -p 1 --dta -x {params.ref_gene} -1 {input.reducedRead1} -2 {input.reducedRead2} -S {output[0]} 2>{output[1]}
        """

rule sam2bam:
    input:
        sam="alignedReads/{err}.sam"
    output:
        "alignedReads/{err}_sorted.bam"
    threads: 1
    params:
        ref_gene = "/local/work/biocore/mol8008/RNA/Human_hg20/genome_tran"
    shell:
        """
        samtools view -bS {input.sam} | samtools sort -o {output}
        """
    
rule featureCounts:
    input:
        bam = "alignedReads/{err}_sorted.bam"
    output:
        bam = "alignedReads/{err}_sorted.bam"
    threads:2
    params:
        ref_trans = "/local/work/biocore/mol8008/RNA/Homo_sapiens.GRCh38.84.gtf"
    shell:
        """
		featureCounts -T 1 -t exon -g gene_id -O -a {params.ref_trans} -o count-1smp.txt {err}_sorted.bam
        """