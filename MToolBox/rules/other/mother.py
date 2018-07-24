import os
from glob import glob
from os.path import join as pjoin
import Bio
from tqdm import tqdm

configfile : "params.json"
workdir : config['home']

include : "read_processing.py"


def all_bams(wildcards):
    samples = all_samples()
    path = "{path}/mapping/bams/".format(path = wildcards.path)
    return [pjoin(path,s + ".bam") for s in samples]



def get_coas_libs(wildcards):
    coas_name = wildcards.name
    with open(pjoin(COAS_FILES_DIR, coas_name + ".txt")) as handle:
        samples = [l.strip() for l in handle]
    fwds = ["1000_processed_reads/{sample}/reads/trimmomatic/{sample}_1P.fastq.gz".format(sample = s) for s in samples]
    revs = ["1000_processed_reads/{sample}/reads/trimmomatic/{sample}_2P.fastq.gz".format(sample = s) for s in samples]
    u_fwds = ["1000_processed_reads/{sample}/reads/trimmomatic/{sample}_1U.fastq.gz".format(sample = s) for s in samples]
    u_revs = ["1000_processed_reads/{sample}/reads/trimmomatic/{sample}_2U.fastq.gz".format(sample = s) for s in samples]
    return {'fwd' : fwds, 'rev' : revs, 'u_rev' : u_revs, 'u_fwd' : u_fwds}

def find_fastq(wildcards):
    path = '0000_raws/0100_reads/genomic/'
    result = [y for x in os.walk(path) for y in glob(os.path.join(x[0], '*.fastq.gz')) if "/" + wildcards.sample +"_" in y]
    assert len(result) == 2, print(result)
    return result

def all_samples():
    path = "1000_processed_reads/"
    samples = [d for d in os.listdir(path) if os.path.isdir(pjoin(path,d)) ]
    return samples

def all_dones(wildcards):
    path = '0000_raws/0100_reads/genomic/'
    result = [pjoin("1000_processed_reads/","_".join(s.split("_")[:-4]),"done") for s in os.listdir(path) if s.endswith(".fastq.gz")  ]
    return list(set(result))

def all_mags(wildcards):
    mag_folder = "{path}/{name}/{assembler}/MAGs/".format(path = wildcards.path, name=wildcards.name, assembler = wildcards.assembler)
    name = wildcards.name
    bins = [f.split("-")[-1][:-3] for f in os.listdir(mag_folder) if f.endswith(".fa")]
    results = [mag_folder + "metabat_{name}_{bin}/{name}_{bin}.faa".format(name=name, bin = f)  for f in bins]

    return results





rule MAG_stats:
    input : bin_file = "{path}/{name}/{assembler}/MAGs/metabat_{name}-unbinned.fa",
            mag_list = all_mags
    output : "{path}/{name}/{assembler}/MAGs/{name}.magstats"
    threads : 1
    run :
        from Bio import SeqIO
        from tqdm import tqdm
        from pandas import DataFrame

        out_dict = {}
        mag_folder = "/".join(input.bin_file.split("/")[:-1])
        bin_files = input.mag_list

        out_file = output[0]

#        with open("phylophlan_tax.txt") as handle:
#            tax = [l.split()[:2] for l in handle.readlines()]
#            tax = {t[0] : ":".join([tt for tt in t[1].split(".") if "?" not in tt ]) for t in tax if "all_samples-" in t[0]}
        def process_bin(binl) :
            bin_head = binl[:-4]
            bin_checkm = bin_head + ".checkm"
            bin_genome= bin_head + ".fna"
            bin_proteom = bin_head + ".faa"
            bin_id = bin_head.split("/")[-1]
            out_dict = {}

            if "unbinned" not in bin_id:
                checkm_fields = ['Bin Id', 'Marker lineage', 'UID', '# genomes', '# markers', '# marker sets', '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity']
                with open(bin_checkm) as handle:
                    dat = handle.readlines()
                dat = {l.split()[0] : {k : v for k,v in zip(checkm_fields, l[:-1].split() ) } for l in dat[3:-1]}
                dat = dat[bin_id]

                out_dict['completeness'] = float(dat['Completeness'])
                out_dict['contamination'] = float(dat['Contamination'])
                out_dict['taxo:checkm'] = dat['Marker lineage']
                out_dict['strain_heterogeneity'] = float(dat['Strain heterogeneity'])
            else :
                out_dict['completeness'] = None
                out_dict['contamination'] = None
                out_dict['taxo:checkm'] = 'Unbinned'
                out_dict['strain_heterogeneity'] = None


            with open( bin_genome) as handle:
                fna = [s for s in SeqIO.parse(handle, "fasta")]
            with open( bin_proteom) as handle:
                faa = [s for s in SeqIO.parse(handle, "fasta")]

            out_dict['length'] = sum([len(s.seq) for s in fna])
            out_dict['nb_contigs'] = len(fna)
            out_dict['nb_proteins'] = len(faa)
            out_dict['coding_density'] = (3.0*sum([len(s.seq) for s in faa]))/sum([len(s.seq) for s in fna])
            out_dict['GC'] = float(sum([str(s.seq).count("G")+str(s.seq).count("C") for s in fna]))/out_dict['length']
            #out_dict['taxo:phylophlan'] = tax[bin_id]
            return (bin_id, out_dict)

        dat = dict([process_bin(p) for p in tqdm(bin_files)])
        DataFrame.from_dict(dat, orient = 'index').to_csv(out_file)



rule all_kaiju:
    input : "1000_processed_reads/done"
    output : "1000_processed_reads/kaiju_table.tsv"
    run :
        from ete3 import NCBITaxa
        from collections import defaultdict
        import os
        from tqdm import tqdm
        from os.path import join as pjoin
        from pandas import DataFrame, concat

        out_file="1000_processed_reads/kaiju_table.tsv"

        pat = "/crex/proj/uppstore2018116/moritz6/1000_processed_reads/"
        samples = [s for s in os.listdir(pat) if "." not in s]
        out_dict = {}

        for s in tqdm(samples):
            if not out_dict.get(s):
                out_dict[s]=defaultdict(int)
                with open(pjoin(pat, s, "reads", "kaiju", s+"_kaiju.out")) as handle:
                    for l in tqdm(handle):
                        out_dict[s][l.split()[2]]+=1

        data = DataFrame.from_dict(out_dict)
        data = data.fillna(0)

        taxDb = NCBITaxa()
        line= {i : taxDb.get_lineage(i) for i in tqdm(list(data.index))}
        out = {i : taxDb.get_taxid_translator(l) if l else None for i,l in tqdm(line.items())}
        tt = sum([list(l.keys()) for l in tqdm(out.values()) if l], [])
        tt = list(set(tt))
        tt = taxDb.get_rank(tt)

        out_l = {k : {tt[k] : v for k,v in t.items()} if t else None for k,t in tqdm(out.items())}
        taxos = DataFrame.from_dict(out_l).transpose()
        taxos = taxos.fillna("NA")
        all = concat([data, taxos.loc[data.index]], axis=1)
	all.to_csv(out_file)







rule library:
    input :           "{path}/{sample}/reads/mash/{sample}.msh",
            "{path}/{sample}/megahit/mapping/mapping_rates.txt",
            "{path}/{sample}/reads/kaiju/{sample}_kaiju.out.summary"
    output : "{path}/{sample}/done"
    shell:
        "touch {wildcards.path}/{wildcards.sample}/done"



rule all_libraries :
    input : all_dones ,
    output : "1000_processed_reads/done"
    shell : "touch 1000_processed_reads/done"
