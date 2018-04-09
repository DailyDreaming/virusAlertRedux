#!/bin/bash
#  setup.sh
#
#  Script slightly modified from: https://github.com/NCBI-Hackathons/VirusFriends
#  Definitely run this with 'sudo'.
#
#  Copyright 2017 The University of Sydney
#  Author: Jan P Buchmann <jan.buchmann@sydney.edu.au>
#  Author: Robert Edwards <raedwards@gmail.com>
#  Author: Bhavya Papudeshi <npbhavya13@gmail.com>
#
#  Description:
#   Install some dependencies
#
SCRIPT_DIR="$( cd "$(dirname "$0")" ; pwd -P )"
DATA_DIR="$SCRIPT_DIR/../data"
TOOLS_DIR="$SCRIPT_DIR/../tools"
VF_DIR="$TOOLS_DIR/VirFinder"
export R_USER_LIBS="${TOOLS_DIR}/R"
export PATH=$TOOLS_DIR:$PATH
echo "TOOLS_DIR: $TOOLS_DIR"
echo "DATA_DIR: $DATA_DIR"

INSTALL=1

mkdir -p $TOOLS_DIR
mkdir -p $DATA_DIR
sudo chmod 777 $TOOLS_DIR
sudo chmod 777 $DATA_DIR

NEWPATH=""

function setup_edirect()
{
	## check for edirect
	esearch=$(which esearch)
	efetch=$(which efetch)
	xtract=$(which xtract)
	if [ -z $esearch ] || [ -z $efetch ] || [ -z $xtract ]; then
		## edirect
		if [ $INSTALL == 1 ]; then
			echo "INSTALLING edirect in $endovir_tools/edirect"
			cd $TOOLS_DIR
			perl -MNet::FTP -e \
			    '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
			     $ftp->login; $ftp->binary;
			     $ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
			     gunzip -c edirect.tar.gz | tar xf -
			     rm -f edirect.tar.gz
			     ./edirect/setup.sh
			NEWPATH=$NEWPATH:$PWD
			cd $SCRIPT_DIR
		fi
	fi
}


function setup_blast()
{
	makeprofiledb=$(which makeprofiledb)
	if [ -z $makeprofiledb ]; then
		## blast+
		if [ $INSTALL == 1 ]; then
			echo "Installing NCBI blast+";
			cd $TOOLS_DIR
			wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-*-x64-linux.tar.gz -O blast.tgz
			tar xf blast.tgz
			rm -f blast.tgz
			# bit of a hack to find the path name because blast includes the version number
			P=$(find . -name blastn | sed -e 's/blastn$//; s/^\.\///')
			NEWPATH=$NEWPATH:$PWD/$P
			cd $SCRIPT_DIR
		fi
	fi
}

function setup_magicblast()
{
	magicblast=$(which magicblast)
	if [ -z $magicblast ]; then
		if [ $INSTALL == 1 ]; then
			echo "Installing magicblast"
			cd $TOOLS_DIR
			wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/magicblast/LATEST/ncbi-magicblast-1.3.0-x64-linux.tar.gz -O magicblast.tar.gz
			tar  -xvf magicblast.tar.gz
			rm -f magicblast.tar.gz
			P=$(find . -name magicblast | sed -e 's/magicblast$//; s/^\.\///')
			NEWPATH=$NEWPATH:$PWD/$P
			cd $SCRIPT_DIR
		fi
	fi
}

function setup_samtools()
{
	samtools_v=$(samtools --version)
	version=$( echo $samtools_v | cut -d ' ' -f2 )
	if [ -z $version != 1.7 ]; then
                if [ $INSTALL == 1 ]; then
                        echo "Installing SAMtools"
                        cd $TOOLS_DIR
                        wget  https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 -O samtools.tar.bz2
                        tar  -vxjf samtools.tar.bz2
                        rm -f samtools.tar.bz2
			cd samtools-1.7
			./configure --prefix=$TOOLS_DIR/samtools-1.7/
			make && make install
			cd ..
			P=$(find . -name samtools | sed -e 's/samtools$//; s/^\.\///')
                        NEWPATH=$NEWPATH:$PWD/$P
                        cd $SCRIPT_DIR
                fi
        fi
}

function setup_sratoolkit()
{
        sratools=$(which vdb_dump)
        if [ -z $sratools ]; then
                if [ $INSTALL == 1 ]; then
                        echo "Installing SRAtoolkit"
                        cd $TOOLS_DIR
                        wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.0/sratoolkit.2.9.0-centos_linux64.tar.gz -O sratoolkit.tar.gz
                        tar  -xvzf sratoolkit.tar.gz
                        rm -f sratoolkit.tar.gz
                        P=$(find . -name vdb_dump | sed -e 's/sratools$//; s/^\.\///')
                        NEWPATH=$NEWPATH:$PWD/$P
                        cd $SCRIPT_DIR
                fi
        fi
}

function make_endovir_cdd()
{
  echo "Creating a Viral Profile Database of Conserved Domains."
  if [ $INSTALL == 1 ]; then
  	echo "" > "$DATA_DIR/$endovir_pssms"
  	local qry="txid10239[Organism:exp] NOT (predicted OR putative)"
  	for i in $($esearch -db  cdd -query "$qry"                      | \
             $efetch -format docsum                               | \
             $xtract -pattern DocumentSummary -element Accession  | \
             grep -v cl)
    	do
      	     echo $i".smp" >> "$DATA_DIR/$endovir_pssms"
    	done

  	local cdd_ftp='ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cdd.tar.gz'
  	wget $cdd_ftp -O - | tar -C "$DATA_DIR/" -xzvT "$DATA_DIR/$endovir_pssms"
  	cd $DATA_DIR
  	makeprofiledb -title endovir -in $DATA_DIR/Cdd.pn   \
                 -out $DATA_DIR/endovir_cdd -dbtype rps -threshold 9.82 -scale 100 -index true
  fi
}


setup_blast
setup_magicblast
setup_sratoolkit
setup_samtools
make_endovir_cdd

####################
## R Installation ##
####################

echo "installing debian dependencies..."
sudo apt-get install python3-biopython sra-toolkit r-base libxml2 libxml2-dev libcurl4-openssl-dev libssl-dev

echo "installing R dependencies to ${TOOLS_DIR}/R..."
mkdir -p "$R_USER_LIBS"
echo """
# see https://www.r-bloggers.com/permanently-setting-the-cran-repository/
local({
  r <- getOption('repos')
  r['CRAN'] <- 'http://cran.cnr.berkeley.edu/'
  options(repos = r)
})

source('https://bioconductor.org/biocLite.R')
biocLite()

install.packages('tidyverse', dependencies=TRUE, lib='$R_USER_LIBS')
install.packages('glmnet', dependencies=TRUE, lib='$R_USER_LIBS')
install.packages('Rcpp', dependencies=TRUE, lib='$R_USER_LIBS')
""" | R --slave

echo "downloading VirFinder to ${TOOLS_DIR}/VirFinder..."
pushd "$TOOLS_DIR"
if [[ ! -d VirFinder ]]; then
  wget https://codeload.github.com/jessieren/VirFinder/zip/master
  unzip master
  mv VirFinder-master/linux/VirFinder ./
  rm -r *master*
fi
popd

echo "installing VirFinder R package..."
echo """
install.packages('$VF_DIR', repos=NULL, type='source', lib='$R_USER_LIBS')
""" | R --slave

popd
