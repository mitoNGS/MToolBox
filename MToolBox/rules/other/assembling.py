rule megahit_single:
    params : temp_folder = TEMP_DIR
    input : fwd = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_1P.fastq.gz",
            rev = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_2P.fastq.gz",
            u_rev = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_2U.fastq.gz",
            u_fwd = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_1U.fastq.gz"
    output : assembly = "1000_processed_reads/{sample}/megahit/{sample}.fna",
             folder = "1000_processed_reads/{sample}/megahit/data",
             foldermain = "1000_processed_reads/{sample}/megahit"
    threads : THREADS
    shell : """
        unpigz -c {input.fwd}  >  {params.temp_folder}/temp_R1.fastq
        unpigz -c {input.rev}  >  {params.temp_folder}/temp_R2.fastq
        unpigz -c {input.u_rev}  >  {params.temp_folder}/temp_U.fastq
        unpigz -c {input.u_fwd}  >>  {params.temp_folder}/temp_U.fastq
        megahit  -1 {params.temp_folder}/temp_R1.fastq -2 {params.temp_folder}/temp_R2.fastq -r {params.temp_folder}/temp_U.fastq -t {threads} -o {params.temp_folder}/temp_data --out-prefix megahit --continue
        rm -r {params.temp_folder}/temp_data/intermediate_contigs/
        mv {params.temp_folder}/temp_data/ {output.folder}
        cp 1000_processed_reads/{wildcards.sample}/megahit/data/megahit.contigs.fa {output.assembly}
    """

rule spades_single:
    params : temp_folder = TEMP_DIR
    input : fwd = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_1P.fastq.gz",
            rev = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_2P.fastq.gz",
            u_rev = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_2U.fastq.gz",
            u_fwd = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_1U.fastq.gz"
    output : assembly = "1000_processed_reads/{sample}/spades/{sample}.fna",
             folder = "1000_processed_reads/{sample}/spades/data",
             foldermain = "1000_processed_reads/{sample}/spades"
    threads : THREADS
    shell : """
        unpigz -c {input.fwd}  >  {params.temp_folder}/temp_R1.fastq
        unpigz -c {input.rev}  >  {params.temp_folder}/temp_R2.fastq
        unpigz -c {input.u_rev}  >  {params.temp_folder}/temp_U.fastq
        unpigz -c {input.u_fwd}  >>  {params.temp_folder}/temp_U.fastq
        spades.py --meta -1 {params.temp_folder}/temp_R1.fastq -2 {params.temp_folder}/temp_R2.fastq -s {params.temp_folder}/temp_U.fastq -t {threads} -o {params.temp_folder}/temp_data
        mv {params.temp_folder}/temp_data/ {output.folder}
        cp 1000_processed_reads/{wildcards.sample}/spades/data/scaffolds.fasta {output.assembly}
    """


rule megahit_coas:
    params : temp_folder = TEMP_DIR
    input : unpack(get_coas_libs)
    output :    folder = "1500_assemblies/{name}/megahit/data/",
                assembly = "1500_assemblies/{name}/megahit/{name}.fna",
                out = "1500_assemblies/{name}/megahit"
    threads : THREADS
    shell : """
         unpigz -c {input.fwd}  >  {params.temp_folder}/temp_R1.fastq
         unpigz -c {input.rev}  >  {params.temp_folder}/temp_R2.fastq
         unpigz -c {input.u_rev}  >  {params.temp_folder}/temp_U.fastq
         unpigz -c {input.u_fwd}  >>  {params.temp_folder}/temp_U.fastq
        megahit  -m 0.8 -1 {params.temp_folder}/temp_R1.fastq -2 {params.temp_folder}/temp_R2.fastq -r {params.temp_folder}/temp_U.fastq -t {threads} -o {params.temp_folder}/temp_data --out-prefix megahit --continue
        rm -r {params.temp_folder}/temp_data/intermediate_contigs/
        mv {params.temp_folder}/temp_data/* {output.folder}
        cp 1500_assemblies/{wildcards.name}/megahit/data/megahit.contigs.fa {output.assembly}
    """

rule spades_coas:
    params : temp_folder = TEMP_DIR
    input : unpack(get_coas_libs)
    output :    folder = "1500_assemblies/{name}/spades/data/",
                assembly = "1500_assemblies/{name}/spades/{name}.fna",
                path = "1500_assemblies/{name}/spades"
    threads : THREADS
    shell : """
         unpigz -c {input.fwd}  >  {params.temp_folder}/temp_R1.fastq
         unpigz -c {input.rev}  >  {params.temp_folder}/temp_R2.fastq
         unpigz -c {input.u_rev}  >  {params.temp_folder}/temp_U.fastq
         unpigz -c {input.u_fwd}  >>  {params.temp_folder}/temp_U.fastq
        spades.py --meta -1 {params.temp_folder}/temp_R1.fastq -2 {params.temp_folder}/temp_R2.fastq -s {params.temp_folder}/temp_U.fastq -t {threads} -o {params.temp_folder}/temp_data

        mv {params.temp_folder}/temp_data/* {output.folder}
        cp 1000_processed_reads/{wildcards.name}/spades/data/scaffolds.fasta {output.assembly}

    """
