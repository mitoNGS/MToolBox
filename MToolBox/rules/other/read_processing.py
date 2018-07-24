import pandas as pd
#from snakemake.utils import validate
import os, glob, subprocess

configfile : "params.json"
samples = pd.read_table(config["general"]["samples"], comment='#').set_index("sample", drop=False)

def getFastqNames(wildcards):
    path="data/raw_reads"
    #return samples.loc[(wildcards.sample), ['R1', 'R2']].dropna()
    L = samples.loc[wildcards.sample, ['R1', 'R2']].dropna()
    M = []
    for i in L:
        M.append(path + "/" + i)
    return M

# def find_sample_fastq_files(sample_list):
#     fastq_files = glob.glob(config['general']['home']+"/data/raw_reads/*fastq.gz")
#     fastq_files.sort()
#     fastq_reads = [fastq_file.split('/')[-1].replace('.fastq.gz', '') for fastq_file in fastq_files]
#     return fastq_reads
# 
# def find_fastq(wildcards):
#     path = '/proj/uppstore2018116/moritz6/0000_raws/0100_reads/genomic/'
#     result = [y for x in os.walk(path) for y in glob.glob(os.path.join(x[0], '*.fastq.gz')) if "/" + wildcards.sample +"_" in y]
#     assert len(result) == 2, print(result)
#     return result

# def find_fastq(wildcards):
#     path="data/raw_reads"
#     #result = [y for x in os.walk(path) for y in glob.glob(os.path.join(x[0], '*.fastq.gz')) if "/" + sample +"_" in y]
#     result = [y for x in os.walk(path) for y in glob.glob(os.path.join(x[0], '*.fastq.gz')) if "/" + {wildcards.sample} +"_" in y]
#     print(result)
#     #result = ' '.join([y for x in os.walk(path) for y in glob.glob(os.path.join(x[0], '*.fastq.gz')) if "/" + {wildcards.sample} +"_" in y])
#     assert len(result) == 2, print(result)
#     return result


# 
# 
# def all_dones(wildcards):
#     path = '0000_raws/0100_reads/genomic/'
#     result = [pjoin("1000_processed_reads/","_".join(s.split("_")[:-4]),"done") for s in os.listdir(path) if s.endswith(".fastq.gz")  ]
#     return list(set(result))


# rule fastqc:
#     #input :  find_fastq
#     input : lambda wildcards: find_fastq(wildcards)
#     output : "1000_processed_reads/{sample}/reads/fastqc"
#     params : threads = config['general']['threads']
#     #threads : THREADS
#     shell:
#         """
#         fastqc -t {params.threads} -o {output} {input}
#         """

#READS=collect_reads()
#print("output of function:\n{}".format(",".join(READS)))
#print("\noutput of expand:")
#print(expand("data/raw_reads/{read}.fastq.gz", read=READS))

ALL_FASTQ_FILES=list(samples['R1']) + list(samples['R2'])
FILES_R1_R2=zip(list(samples['R1']),list(samples['R2']))

rule all:
    input:
        expand("results/mash/{sample}.msh", sample=samples["sample"]),
        expand("data/processed_reads/trimmomatic/{sample}_{fwd_rev}{PE_SE}.fastq.gz", sample=samples["sample"], fwd_rev=[1, 2], PE_SE=["P", "U"]),
        expand("results/fastqc/raw_reads/{sample_seqfile}_fastqc", sample_seqfile=ALL_FASTQ_FILES)

rule fastqc:
    input:
        "data/raw_reads/{sample_seqfile}"
    output:
        "results/fastqc/raw_reads/{sample_seqfile}_fastqc"
    threads:
        config['general']['threads']
    version:
        subprocess.getoutput(
            "fastqc -V"
            )
    message:
        "QC of of read file {input} with {version}"
    shell:
        """

        fastqc -t {threads} -o {output} {input}

        """

rule trimmomatic:
    """ QCing and cleaning reads """
    params:
        java_cmd = config['read_processing']['trimmomatic']['java_cmd'],
        jar_file = config['read_processing']['trimmomatic']['jar_file'],
        mem = config['read_processing']['trimmomatic']['java_vm_mem'],
        options = config['read_processing']['trimmomatic']['options'],
        processing_options = config['read_processing']['trimmomatic']['processing_options'],
        temp_folder = config['general']['temp_dir']
    input:
        getFastqNames
        #in1P = "data/raw_reads" + getFastqNames[0],
        #in2P = "data/raw_reads" + getFastqNames[1]
    output:
        out1P = "data/processed_reads/trimmomatic/{sample}_1P.fastq.gz",
        out2P = "data/processed_reads/trimmomatic/{sample}_2P.fastq.gz",
        out1U = "data/processed_reads/trimmomatic/{sample}_1U.fastq.gz",
        out2U = "data/processed_reads/trimmomatic/{sample}_2U.fastq.gz"
    threads:
        config['general']['threads']
    version:
        subprocess.getoutput(
            "java -jar $TRIMMOMATIC_HOME/trimmomatic.jar PE -version"
            )
    message:
        "Filtering read datasets for sample {wildcards.sample} with Trimmomatic v{version}"
    log: 
        "data/processed_reads/trimmomatic/{sample}.log"
    shell:
        """

        {params.java_cmd} -Xmx{params.mem} -Djava.io.tmpdir={params.temp_folder} -jar {params.jar_file} PE \
            {params.options} \
            -threads {threads} \
            data/raw_reads/{input[0]} data/raw_reads/{input[1]} \
            {output.out1P} {output.out1U} {output.out2P} {output.out2U} \
            {params.processing_options} 2> {log}

        """

rule mash:
    params : kmer = config['read_processing']['mash']['kmer'],
             hashes = config['read_processing']['mash']['hashes'],
    input : 
        #"1000_processed_reads/{sample}/reads/trimmomatic/{sample}_1P.fastq.gz","1000_processed_reads/{sample}/reads/trimmomatic/{sample}_2P.fastq.gz"
        in1P = "data/processed_reads/trimmomatic/{sample}_1P.fastq.gz",
        in2P = "data/processed_reads/trimmomatic/{sample}_2P.fastq.gz"
    output :
        #"1000_processed_reads/{sample}/reads/mash/{sample}.msh"
        "results/mash/{sample}.msh"
    log : "logs/mash/{sample}.log"
    threads : config['general']['threads']
    shell :
        "mash sketch -r -p {threads} -k {params.kmer} -s {params.hashes} -o $(echo '{output}' | sed -e 's/.msh//') {input} > {log}"
# 
# rule kaiju:
#     params : db_path = config['read_processing']['kaiju']['db_path'],
#              db = config['read_processing']['kaiju']['db'],
#              #db = config['kaiju']['uppmax']
#     input : "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_1P.fastq.gz","1000_processed_reads/{sample}/reads/trimmomatic/{sample}_2P.fastq.gz"
#     output : "1000_processed_reads/{sample}/reads/kaiju/{sample}_kaiju.out.summary", "1000_processed_reads/{sample}/reads/kaiju/{sample}_kaiju.html"
#     log : "1000_processed_reads/{sample}/reads/kaiju/{sample}_kaiju.log"
#     threads : config['general']['threads']
#     shell : """
#     module load bioinfo-tools
#     module load Krona
#     kaiju -t {params.db_path}/nodes.dmp -f {params.db_path}/{params.db}   -i 1000_processed_reads/{wildcards.sample}/reads/trimmomatic/{wildcards.sample}_1P.fastq.gz -j 1000_processed_reads/{wildcards.sample}/reads/trimmomatic/{wildcards.sample}_2P.fastq.gz -o 1000_processed_reads/{wildcards.sample}/reads/kaiju/{wildcards.sample}_kaiju.out -z {threads} > {log}
#     kaiju2krona -u -t {params.db_path}/nodes.dmp -n {params.db_path}/names.dmp -i 1000_processed_reads/{wildcards.sample}/reads/kaiju/{wildcards.sample}_kaiju.out -o 1000_processed_reads/{wildcards.sample}/reads/kaiju/{wildcards.sample}_kaiju.out.krona >> {log}
#     ktImportText  -o 1000_processed_reads/{wildcards.sample}/reads/kaiju/{wildcards.sample}_kaiju.html 1000_processed_reads/{wildcards.sample}/reads/kaiju/{wildcards.sample}_kaiju.out.krona >> {log}
#     kaijuReport -p -r genus -t {params.db_path}/nodes.dmp -n {params.db_path}/names.dmp -i 1000_processed_reads/{wildcards.sample}/reads/kaiju/{wildcards.sample}_kaiju.out -r family -o 1000_processed_reads/{wildcards.sample}/reads/kaiju/{wildcards.sample}_kaiju.out.summary >> {log}
#     """
# 
# rule matam:
#     params : db_path = config['read_processing']['matam']['db_path'],
#              temp_dir = config['general']['temp_dir'],
#              max_mem = config['read_processing']['matam']['max_mem']
#     threads : config['general']['threads']
#     input :  fwd = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_1P.fastq.gz",
#              rev = "1000_processed_reads/{sample}/reads/trimmomatic/{sample}_2P.fastq.gz"
#     output : dir =  "1000_processed_reads/{sample}/reads/matam/",
#              fasta =  "1000_processed_reads/{sample}/reads/matam/final_assembly.fa",
#     shell : """
#     paste <(unpigz  -p {threads} -c {input.fwd} | paste  - - - -) <(unpigz {threads} -c {input.rev} | paste - - - -) | tr $"\t" $"\n" >  {params.temp_dir}/matam_{wildcards.sample}.fastq
#     matam_assembly.py -d {params.db_path} -i {params.temp_dir}/matam_{wildcards.sample}.fastq --cpu {threads} --max_memory {params.max_mem} -v --perform_taxonomic_assignment -o {output.dir}
#     """
# 
# 
# rule all_kaiju:
#     input : "1000_processed_reads/done"
#     output : "1000_processed_reads/kaiju_table.tsv"
#     run :
#         from ete3 import NCBITaxa
#         from collections import defaultdict
#         import os
#         from tqdm import tqdm
#         from os.path import join as pjoin
#         from pandas import DataFrame, concat
# 
#         out_file="1000_processed_reads/kaiju_table.tsv"
# 
#         pat = "/crex/proj/uppstore2018116/moritz6/1000_processed_reads/"
#         samples = [s for s in os.listdir(pat) if "." not in s]
#         out_dict = {}
# 
#         for s in tqdm(samples):
#             if not out_dict.get(s):
#                 out_dict[s]=defaultdict(int)
#                 with open(pjoin(pat, s, "reads", "kaiju", s+"_kaiju.out")) as handle:
#                     for l in tqdm(handle):
#                         out_dict[s][l.split()[2]]+=1
# 
#         data = DataFrame.from_dict(out_dict)
#         data = data.fillna(0)
# 
#         taxDb = NCBITaxa()
#         line= {i : taxDb.get_lineage(i) for i in tqdm(list(data.index))}
#         out = {i : taxDb.get_taxid_translator(l) if l else None for i,l in tqdm(line.items())}
#         tt = sum([list(l.keys()) for l in tqdm(out.values()) if l], [])
#         tt = list(set(tt))
#         tt = taxDb.get_rank(tt)
# 
#         out_l = {k : {tt[k] : v for k,v in t.items()} if t else None for k,t in tqdm(out.items())}
#         taxos = DataFrame.from_dict(out_l).transpose()
#         taxos = taxos.fillna("NA")
#         all = concat([data, taxos.loc[data.index]], axis=1)
# 	all.to_csv(out_file)
# 
# 
# 
# 
# 
# 
# 
# rule library:
#     input :           "{path}/{sample}/reads/mash/{sample}.msh",
#             "{path}/{sample}/megahit/mapping/mapping_rates.txt",
#             "{path}/{sample}/reads/kaiju/{sample}_kaiju.out.summary"
#     output : "{path}/{sample}/done"
#     shell:
#         "touch {wildcards.path}/{wildcards.sample}/done"
# 
# 
# 
# rule all_libraries :
#     input : all_dones ,
#     output : "1000_processed_reads/done"
#     shell : "touch 1000_processed_reads/done"
