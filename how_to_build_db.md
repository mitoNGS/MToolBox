###How to build GSNAP databases and fasta index used by MToolBox

1) Be sure to have GSNAP and samtools installed on your machine:

```which gsnap```

```which samtools```

The current MToolBox pipeline versions only support samtools versions lower than 1. Please visit [this site](http://sourceforge.net/projects/samtools/files/samtools/), to get samtools versions lower than 1 (i.e. 0.1.xx) 

If GSNAP is not installed, please have a look at the [gsnap website](http://research-pub.gene.com/gmap/), to get the latest or one of the previous GSNAP versions 

2) Create the directories where you want to place your GSNAP database and fasta indexes:

```mkdir gmapdb```

```mkdir genomes```

3) For user's convenience, fasta files of rCRS/RSRS mitochondrial reference sequence and nuclear reference sequence (GRCh37/hg19) can be dowloaded from the [MToolBox sourceforge page](https://sourceforge.net/projects/mtoolbox/). If you have wget already installed:

```
cd genomes

wget http://sourceforge.net/projects/mtoolbox/files/genome_fasta/chrM.fa.tar.gz

wget http://sourceforge.net/projects/mtoolbox/files/genome_fasta/hg19RCRS.fa.tar.gz

wget http://sourceforge.net/projects/mtoolbox/files/genome_fasta/chrRSRS.fa.tar.gz

wget http://sourceforge.net/projects/mtoolbox/files/genome_fasta/hg19RSRS.fa.tar.gz
```
then decompress files...

```for i in `ls *.tar.gz` ; do tar xvzf $i; done```

and remove tar.gz files

```rm *tar.gz```

```cd ..```

4) Generate the GSNAP databases in the gmapdb folder by doing (time consuming step):

```cd gmapdb```

```gmap_build -D . -d hg19RCRS  ../genomes/hg19RCRS.fa -s numeric-alpha ```

```gmap_build -D . -d chrM ../genomes/chrM.fa -s numeric-alpha```

```gmap_build -D . -d chrRSRS ../genomes/chrRSRS.fa -s numeric-alpha ``

```gmap_build -D . -d hg19RSRS ../genomes/hg19RSRS.fa -s numeric-alpha```

```cd ..```

With this command line you are using the default kmer length (15). To choose the best kmer length for your machine memory requirements, please have a look at the GMAP/GSNAP [readme](http://research-pub.gene.com/gmap/src/README)

**NOTE ON GSNAP DATABASE FOR MTOOLBOX**: If you do not want to use the mitochondrial fasta sequences available from the [MToolBox sourceforge page](https://sourceforge.net/projects/mtoolbox/files/genome_fasta/) and prefer to use your mitochondrial fasta reference, please be sure to use the **same** mitochondrial chromosome name of the fasta header (e.g. *>chrMT*) as mitochondrial GSNAP database name (*chrMT* in this example)

5) Create fasta index with samtools

```cd genomes```

```for i in `ls *.fa`; do samtools faidx $i; done```



