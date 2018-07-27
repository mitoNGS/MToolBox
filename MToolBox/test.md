# Test MToolBox snakemake

Samples archive:  

`/lustre/home/attimonelli/Epatocellular_carcinoma/new_bam_files/RA_1/OutputRA1/processed_fastq.tar.gz`

Working directory:  


Extract samples:

```bash
cd /lustrehome/domenico.simone/browser_shared/mtoolbox_snakemake/test/data_ros

qsub<<'EOF'
#!/bin/bash
#PBS -q bigmpi2@sauron.recas.ba.infn.it
#PBS -l nodes=1:ppn=1
# PBS -l walltime=48:00:00
# PBS -N 
# PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} 
# PBS -e ${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -m abe -M dome.simone@gmail.com

cd /lustrehome/domenico.simone/browser_shared/mtoolbox_snakemake/test/data_ros

tar --skip-old-files -xvzf \
/lustre/home/attimonelli/Epatocellular_carcinoma/new_bam_files/RA_1/OutputRA1/processed_fastq.tar.gz

cd processed_fastq

for i in $(ls *.fastq); do
    gzip $i
done

EOF
```

Symlink files to directory `00_raw_reads`

```bash
cd /lustrehome/domenico.simone/browser_shared/mtoolbox_snakemake/test/
ln -s data_ros/processed_fastq/TCGA_CC_A7I*fastq.gz 00_raw_reads/
```

Compile file `samples.tsv` file

```bash
cd /lustrehome/domenico.simone/browser_shared/mtoolbox_snakemake/test/
printf "sample\tR1\tR2\n" > samples.tsv
ls TCGA_CC_A7I* | \
awk 'BEGIN{FS="."}{print $1}' | \
sort | uniq | \
awk 'BEGIN{OFS="\t"}{print $1, $1".R1.fastq.gz", $1".R2.fastq.gz"}' >> samples.tsv
```

Run snakemake. Don't forget to activate the conda env `snakemake3`!

```bash
conda activate snakemake3

nohup snakemake -rp -j 50 --cluster-config cluster.json --cluster 'qsub -V -q bigmpi2@sauron.recas.ba.infn.it \
-l nodes=1:ppn=1' &
```

The process fails for 5 samples:

```
OUT_TCGA_CC_A7II_01A_11D/02_map_exome/:
OUT_TCGA_CC_A7IK_10A_01D/02_map_exome/:
OUT_TCGA_CC_A7IH_01A_11D/02_map_exome/:
OUT_TCGA_CC_A7II_10A_01D/02_map_exome/:
OUT_TCGA_CC_A7IG_01A_11D/02_map_exome/:
```

#### Try to run the workflow again to check if the same files fail

```bash
nohup snakemake -rp -j 50 --cluster-config cluster.json --cluster 'qsub -V -q bigmpi2@sauron.recas.ba.infn.it \
-l nodes=1:ppn=1' &
```

Failed (starred ones didn't fail before):

```
OUT_TCGA_CC_A7IK_10A_01D/02_map_exome/:
OUT_TCGA_CC_A7IE_01A_21D/02_map_exome/: *
OUT_TCGA_CC_A7IH_10A_01D/02_map_exome/: *
OUT_TCGA_CC_A7II_01A_11D/02_map_exome/:
OUT_TCGA_CC_A7IJ_01A_11D/02_map_exome/: *
OUT_TCGA_CC_A7IK_01A_12D/02_map_exome/: *
```

Succeeded (starred ones had failed before):

```
OUT_TCGA_CC_A7IF_01A_11D/02_map_exome/:
OUT_TCGA_CC_A7IH_01A_11D/02_map_exome/: *
OUT_TCGA_CC_A7IL_01A_11D/02_map_exome/: 
```

#### Try to run the pipeline with system gsnap

Generate copy of pipeline folder

```bash
cd /lustre/browser/mtoolbox_snakemake
cp -R MToolBox MToolBox_copy && cd MToolBox_copy
```

Edit `map_exome.rules` and `config.yaml` file in order to have:
- `/lustre/browser/bin/MToolBox/bin/gmap/bin/gsnap` as gsnap exec
- `/lustre/browser/bin/MToolBox/gmapdb` as gsnap_db folder.

Generate copy (with symlinks) of test directory

```bash
cd /lustre/browser/mtoolbox_snakemake

mkdir test_native_gsnap && cd test_native_gsnap

ln -s ../MToolBox_copy/MToolBox/Snakefile .
ln -s ../MToolBox_copy/MToolBox/config.yaml .
ln -s ../MToolBox_copy/MToolBox/cluster.json .

cp ../test/samples.tsv .

mkdir 00_raw_reads && cd 00_raw_reads
ln -s ../../test/data_ros/processed_fastq/TCGA_CC_A7I* .

mkdir 01_processed_reads && cd 01_processed_reads
ln -s ../../test/01_processed_reads/TCGA_CC_A7I* .
cd ..

mkdir 00_raw_reads_qc && cd 00_raw_reads_qc
ln -s ../../test/00_raw_reads_qc/TCGA_CC_A7I* .
cd ..

mkdir logs && cd logs
ln -s ../../test/logs/* .
cd ..
```

Run

```bash
nohup snakemake -rp -j 50 --cluster-config cluster.json --cluster 'qsub -V -q bigmpi2@sauron.recas.ba.infn.it \
-l nodes=1:ppn=1' &
```
