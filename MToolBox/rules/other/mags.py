rule phylophlan :
    params : phylophlan_path = "/home/moritz/repos/github/phylophlan",
             phylophlan_exe = "phylophlan.py",
             phylophlan_taxa = "data/ppafull.tax.txt",
             default_genomes = "/home/moritz/repos/github/phylophlan/input/default_genomes"
    input : path = "{path}/{name}/bins/{assembler}/{type}MAGs", annotated = "{path}/{name}/bins/{assembler}/annotated"
    output : "{path}/{name}/bins/{assembler}/{type}{name}.tree"
    threads : THREADS
    shell : """
    DD=`pwd`
    mkdir {params.phylophlan_path}/input/{wildcards.type}{wildcards.name}
    cp {input.path}/*/*.faa {params.phylophlan_path}/input/{wildcards.type}{wildcards.name}
    cd {params.phylophlan_path}/input/{wildcards.type}{wildcards.name}
    for f in `ls *.faa`
    do
        sed -i 's/*//g' $f
    done
    cp {params.default_genomes}/*.faa .
    cd ../..
    python2.7 phylophlan.py --nproc {threads}  -i -t  {wildcards.type}{wildcards.name}
    IFS=$"\n"; for r in `cat data/ppafull.tax.txt`; do id=`echo ${{r}} | cut -f1`; tax=`echo ${{r}} | cut -f2`; sed -i "s/${{id}}/${{id}}_${{tax}}/g" output/{wildcards.type}{wildcards.name}/{wildcards.type}{wildcards.name}.tree.int.nwk; done; unset IFS
    cp output/{wildcards.type}{wildcards.name}/{wildcards.type}{wildcards.name}.tree.int.nwk $DD/{output}
    """


rule annotate_all_mags :
    input : "{path}/{name}/{assembler}/MAGs/metabat_{name}-unbinned.fa"
    output : "{path}/{name}/{assembler}/MAGs/metabat_{name}_unbinned/{name}_unbinned.checkm"
    threads : THREADS
    shell : """

    input_dir=`dirname {input}`

    for b in `ls $input_dir/*.fa`
    do
        bin_id=`basename $b | sed "s#metabat_{wildcards.name}-\\(.*\\).fa#\\1#"`
        out_dir=$input_dir/metabat_{wildcards.name}_$bin_id

        mkdir -p $out_dir

        if [ $bin_id = 'unbinned' ]
        then
            prokka --outdir /scratch/{wildcards.name}_$bin_id/  --metagenome --force --prefix {wildcards.name}_$bin_id --locustag {wildcards.name}_$bin_id --cpus {threads} $b
            mv /scratch/{wildcards.name}_$bin_id/* $out_dir/
            touch $out_dir/{wildcards.name}_$bin_id.checkm
        else
            prokka --outdir /scratch/{wildcards.name}_$bin_id/  --force --prefix {wildcards.name}_$bin_id --locustag {wildcards.name}_$bin_id --cpus {threads} $b
            mv /scratch/{wildcards.name}_$bin_id/* $out_dir/

            checkm lineage_wf -t {threads} -x fna $out_dir $out_dir/data > $out_dir/{wildcards.name}_$bin_id.checkm
            rm -r $out_dir/data
        fi
    done
    """



rule filter_good_MAGs :
# cutoff based on https://www.microbe.net/2017/12/13/why-genome-completeness-and-contamination-estimates-are-more-complicated-than-you-think/
    params : contamination = 5, completeness = 40
    input : path = "{path}/{name}/bins/{assembler}/MAGs",
            anots = "{path}/{name}/bins/{assembler}/annotated",
            checkm = "{path}/{name}/bins/{assembler}/MAGs/{name}.checkm"
    output : "{path}/{name}/bins/{assembler}/good_MAGs"
    run :
        from os.path import join as pjoin
        import os
        import shutil

        checkm_fields = ['Bin Id', 'Marker lineage', 'UID', '# genomes', '# markers', '# marker sets', '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity']
        with open(input.checkm) as handle:
            dat = handle.readlines()
        dat = {l.split()[0] : {k : v for k,v in zip(checkm_fields, l[:-1].split() ) } for l in dat[3:-1]}
        good_bact_MAGs = [k for k,v in dat.items() if float(v['Completeness']) > params.completeness and float(v['Contamination']) < params.contamination]

        os.makedirs(output[0])
        for f in good_bact_MAGs:
            shutil.copytree(pjoin(input.path, f), pjoin(output[0],f))



rule MAG_stats:
    input : bin_file = "{path}/{name}/{assembler}/MAGs/metabat_{name}_unbinned/{name}_unbinned.checkm",
    output : "{path}/{name}/{assembler}/MAGs/{name}.magstats"
    threads : 1
    run :
        from Bio import SeqIO
        from tqdm import tqdm
        from pandas import DataFrame

        out_dict = {}
        mag_folder = "/".join(input.bin_file.split("/")[:-2])
        mag_folder = "{path}/{name}/{assembler}/MAGs/".format(path = wildcards.path, name=wildcards.name, assembler = wildcards.assembler)
        bins = [f.split("-")[-1][:-3] for f in os.listdir(mag_folder) if f.endswith(".fa")]
        bin_files = [mag_folder + "metabat_{name}_{bin}/{name}_{bin}.checkm".format(name=wildcards.name, bin = f)  for f in bins]

        out_file = output[0]

#        with open("phylophlan_tax.txt") as handle:
#            tax = [l.split()[:2] for l in handle.readlines()]
#            tax = {t[0] : ":".join([tt for tt in t[1].split(".") if "?" not in tt ]) for t in tax if "all_samples-" in t[0]}
        def process_bin(binl) :
            bin_head = binl[:-7]
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
