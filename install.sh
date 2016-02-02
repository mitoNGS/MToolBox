#!/bin/bash



download()
{
    URL=$1
    OFILE=$2
    wget  --no-check-certificate -c -nc -L "$URL" -O $OFILE
}

muscle_install()
{
	echo "Installing muscle vx.x..."
	muscle_file=muscle3.8.31_i86linux64.tar.gz
	muscle_url=http://www.drive5.com/muscle/downloads3.8.31/$muscle_file
	pushd .
	cd bin
	download $muscle_url $muscle_file
	tar -zxvf $muscle_file
	popd
	echo "Installing muscle vx.x...done."
}

samtools_install()
{
	echo "Installing samtools vx.x..."
	samtools_version=1.3
	samtools_file=samtools-1.3.tar.bz2
	samtools_url=https://github.com/samtools/samtools/releases/download/1.3/$samtools_file
	pushd .
	cd bin
	download $samtools_url $samtools_file
	tar xjvf $samtools_file
	cd samtools-${samtools_version}
	make
	make prefix=samtools/ install
	popd
	echo "Installing samtools vx.x...done."
}

#########################################
DIRECTORY=bin
if [ -d "$DIRECTORY" ]; then
	echo "$DIRECTORY already exists..."
else
	echo "creating $DIRECTORY..."
	mkdir bin
fi

samtools_install
muscle_install

exit 0
