rule clean_metabat:
    input : "{path}/{name}/{assembler}/mapping/metabat/metabat_{name}"
    output : "{path}/{name}/{assembler}/MAGs/metabat_{name}-unbinned.fa"
    run :
        import os
        from Bio import SeqIO
        from os.path import join as pjoin
        from tqdm import tqdm

        ipath = pjoin(os.path.dirname(input[0]))
        opath = output[0]
        os.makedirs(opath)
        for f in tqdm(os.listdir(ipath)):
            if f[-3:] == ".fa":
                with open(pjoin(ipath, f)) as handle:
                    seqs = [s for s in SeqIO.parse(handle, "fasta")]
                zeros = len(str(len(seqs)))
                bin_name = f[:-3]
                for i,s in enumerate(seqs):
                    s.id = bin_name.replace(".","-") + "-" + str(i+1).zfill(zeros)
                    s.description = ""
                with open(pjoin(opath, f[:-3].replace(".","-")+".fa"), "w") as handle:
                    SeqIO.write(seqs, handle, "fasta")


rule metabat :
    params : min_len = 1500,
             min_bin_size = 10000,
             maxP = 93,
             minS = 50
    input : mapping = "{path}/{name}/{assembler}/mapping/binmap_table.tsv"
    output : file = "{path}/{name}/{assembler}/mapping/metabat/metabat_{name}"
    threads : THREADS
    shell : """
    metabat2 --maxP {params.maxP} --minS {params.minS} -m {params.min_len}  -s {params.min_bin_size} -i  {wildcards.path}/{wildcards.name}/{wildcards.assembler}/{wildcards.name}.fna -o {output.file} -a {input.mapping}  --saveCls  --unbinned -t {threads}
    """

rule concoct :
    params : min_len = 1500,
             min_bin_size = 10000,
             maxP = 93,
             minS = 50
    input : mapping = "{path}/{name}/{assembler}/mapping/binmap_table.tsv"
    output : file = "{path}/{name}/{assembler}/mapping/concoct/concoct_{name}",
             concoct_abundances = "{path}/{name}/{assembler}/mapping/concoct/concoct_{name}.tsv"
    threads : THREADS
    shell : """
    columns=`head -n1 {input.mapping} | sed 's/\t/\n/g' | grep -n bam | grep -v var | cut -f1 -d":" | sed 's/$/,/' | tr -d '\n'`
    cut -f1,${columns%%,} -d$'\t' {input.mapping} > {output.concoct_abundances}

    """

rule maxbin :
    params : max_exec = "/home/moritz/share/MaxBin-2.2.5/run_MaxBin.pl"
    input : mapping = "{path}/{name}/{assembler}/mapping/binmap_table.tsv",
            assembly = "{path}/{name}/{assembler}/{name}.fna"
    output : file = "{path}/{name}/{assembler}/mapping/maxbin/maxbin_{name}",
             maxbin_abunds = "{path}/{name}/{assembler}/mapping/maxbin/maxbin_{name}.tsv"
    threads : THREADS
    shell : """
    mkdir `dirname  {output.maxbin_abunds}`/abunds
    columns=`head -n1 {input.assembly} | tr "\t" "\n" | grep -n bam | grep -v var | cut -f1 -d":" | sed 's/$/,/' | tr -d '\n'`
    cut -f1,${{columns%%,}} -d$'\t' {input.mapping} > {output.maxbin_abunds}

    for f in `head -n1 {output.maxbin_abunds} | tr '\t' '\n'| grep -v contigName`
    do n=`head -n1 {output.maxbin_abunds} | tr "\t" "\n" | grep -n $f | cut -f1 -d":"`
    cut -f1,$n bla | grep -v contig > `basedir {output.maxbin_abunds}`/abunds/$f.tsv
    done

    ls `basedir {output.maxbin_abunds}`/abunds/$f.tsv > {output.maxbin_abunds}.lst


    {params.max_exec} -contig {input.assembly} -out {output.file} -abund_list {output.maxbin_abunds}.lst -thread {threads}
    """
