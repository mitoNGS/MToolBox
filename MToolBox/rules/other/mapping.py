
def all_bams(wildcards):
    samples = all_samples()
    path = "{path}/mapping/bams/".format(path = wildcards.path)
    return [pjoin(path,s + ".bam") for s in samples]




rule bbmap_index:
    input : "{path}/{name}/{assembler}/{name}.fna"
    output : folder = "{path}/{name}/{assembler}/mapping",
             gz = "{path}/{name}/{assembler}/mapping/ref/genome/1/chr1.chrom.gz"
    shell : """
    module load bioinfo-tools
    module load bbmap

    bbmap.sh ref={input} path={output.folder}
    """

rule samnple_wise_bbmap :
    input : index = "{path}/mapping/ref/genome/1/chr1.chrom.gz",
            fwd = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_1P.fastq.gz",
            rev = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_2P.fastq.gz"
    output : bam = "{path}/mapping/bams/{sample}.bam",
#             wdups_stats = "{path}/mapping/bams/{sample}_sorted.stats",
#             stats = "{path}/mapping/bams/{sample}.stats"
    threads : THREADS
    shell : """
        module load bioinfo-tools
        module load bbmap
        module load samtools

        home=`pwd`
        cd `dirname {input.index}`/../../../

        bbmap.sh  in=$home/{input.fwd} in2=$home/{input.rev} threads={threads} bamscript=/scratch/{wildcards.sample}.sh out=/scratch/{wildcards.sample}.sam

        /scratch/{wildcards.sample}.sh
        sambamba flagstat -t {threads} /scratch/{wildcards.sample}_sorted.bam > $home/{wildcards.sample}_sorted.stats
        samtools rmdup  /scratch/{wildcards.sample}_sorted.bam /scratch/{wildcards.sample}.bam
        samtools index /scratch/{wildcards.sample}.bam
        sambamba flagstat  -t {threads} /scratch/{wildcards.sample}.bam > $home/{wildcards.sample}.stats

        rm /scratch/{wildcards.sample}.sam
        rm /scratch/{wildcards.sample}_sorted.bam*
        mv /scratch/{wildcards.sample}.bam* bams/
    """


#rule summrize_bbmap:
#    input

rule bbmap_all_samples:
    input : "{path}/mapping/ref/genome/1/chr1.chrom.gz", all_bams
    output : "{path}/mapping/map_table.tsv",
    threads : 1
    shell : """
    jgi_summarize_bam_contig_depths --outputDepth {input[0]}/map_table.tsv  --pairedContigs {input[0]}/paired_contigs.tsv  {input[0]}/bams/*.bam
    """

rule bbmap_binning_samples:
    params : home = pjoin(HOME,"1000_processed_reads"),
    input : index = "{path}/mapping/ref/genome/1/chr1.chrom.gz",
    	    path = "{path}/mapping/",
            bams = all_bin_samples
    output : "{path}/mapping/binmap_table.tsv",
    threads : 1
    shell : """
    jgi_summarize_bam_contig_depths --outputDepth {input.path}/binmap_table.tsv  --pairedContigs {input.path}/binpaired_contigs.tsv {input.bams}
    """


rule bbmap_diagnostic:
    params : home = "/crex/proj/uppstore2018116/moritz6/1000_processed_reads",
             reads = 100000
    input : "{path}/mapping/ref/genome/1/chr1.chrom.gz"
    output : "{path}/mapping/mapping_rates.txt",
    threads : THREADS
    shell : """
        module load bioinfo-tools
        module load bbmap
        module load samtools

        home=`pwd`
        cd {wildcards.path}/mapping

        echo $'sample\tmated\tfwd\trev' > $home/{output}

        for s in `echo """ + " ".join(all_samples()) + """ | tr ' ' "\n"`
        do
            echo mapping $s to {wildcards.path}
            base={params.home}/$s/reads/trimmomatic/$s
            bbmap.sh  in=${{base}}_1P.fastq.gz in2=${{base}}_2P.fastq.gz threads={threads} out=/dev/null reads={params.reads} 2> tmp
            echo -n $s $'\t' >>  $home/{output}
            cat tmp | grep -E "mated|mapped" | cut -f2  | tr -d ' ' | tr -d % | tr '\n'  '\t' >>  $home/{output}
            echo >> $home/{output}
            rm tmp
        done
    """
