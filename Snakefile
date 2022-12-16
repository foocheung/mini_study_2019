from os import listdir
import glob
import os

configfile: "config.yaml.single"
(config['samples'])
print (config['samples'])



rule all:
    input:
        expand("outdir/{sample}_output/Aligned.sortedByCoord.out.bam", sample = config['samples']),
	expand("outdir/{sample}_output/Aligned.sortedByCoord.out.bam.sort", sample = config['samples']),
        expand("outdir/{sample}_output/Aligned.sortedByCoord.out.bam.sort.count", sample = config['samples'])
             
rule pass1:
        input:
            R1L1 = '{sample}_L001_R1_001.fastq.gz', # may need adjustment if your fastq file name format is different
            R2L1 = '{sample}_L001_R2_001.fastq.gz',
            refdir = directory('Genome')
        params:
            outdir = 'outdir/{sample}_output'
        output:
            "outdir/{sample}_output/Aligned.sortedByCoord.out.bam"
        conda:
            "mapping.yaml"
        threads: 1 # set the maximum number of available cores
        shell:
            'rm -rf {params.outdir} &&' 
            'mkdir {params.outdir} && ' 
            'cd {params.outdir} && '
            'STAR --runThreadN {threads} '
            '--genomeDir {input.refdir} '
            '--readFilesIn {input.R1L1} {input.R2L1} '
            '--readFilesCommand zcat '
            '--outSAMtype BAM SortedByCoordinate && cd ..'

rule sort:
        input:
            sortbam = 'outdir/{sample}_output/Aligned.sortedByCoord.out.bam'
        params:
            outdir = 'outdir/{sample}_output'
        output:
             outdir = 'outdir/{sample}_output/Aligned.sortedByCoord.out.bam.sort'
        threads:1  
        conda:
            "mapping.yaml"
        shell:
            'samtools  sort -@ 6 -n {input.sortbam} -o  {output.outdir}'



rule count:
    input:
       sortbam='outdir/{sample}_output/Aligned.sortedByCoord.out.bam.sort'
    params:
       outdir ='outdir/{sample}_output'
    output:
       outdir ='outdir/{sample}_output/Aligned.sortedByCoord.out.bam.sort.count'
    conda:
       "mapping.yaml"
    shell:
       'featureCounts -p  -T 10 -a Homo_sapiens.GRCh38.95.gtf -s 1 -o {output.outdir} {input.sortbam} ' 
