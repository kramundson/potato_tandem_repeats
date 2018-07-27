Read this first!

This directory functions as a lab notebook for experiments pertaining to k-mer count
analysis in tetraploid potato. Included are descriptions of datasets, code for
data acquisition and analysis, and interpretations of results.

The experiments described do not include analysis of self-generated datasets. All sequencing data was
retrieved from the NCBI Sequencing Read Archive.

The aim of this project is to explore differences in k-mer content between different tetraploid and
diploid accessions of cultivated potato that may reflect structural, copy number or presence absence
variation between accessions. From previous work, potato is unusual among eukaryotes in that its centromeres
are not solly comprised of tandemly arranged satellite repeats. Rather, some chromosomes have satellite-based
centromeres, while others have centromeres comprised of single- or low-copy sequence. Furthermore, centromere
and flanking pericentromere sequences of homologous chromosomes can differ to such a great extent that they
are either 1) virtually unalignable or 2) presence-absence variations in FISH.

For more reading on potato repeat diversity, see the references located under ```papers```

Instructions for downloading datasets and raw data are in the subdirectory ```download/```
Each experiment has its own subdirectory under ```experiments/```
Analysis of informative or interesting experiments are migrated to ```results/```

TODO: description and role for ```protocols```

First, install necessary software for all experiments using the steps below:

1. Download and install miniconda for Python3

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# install
sh Miniconda3-latest-Linux-x86_64.sh

# review license
# accept license
# accept or change home location
# yes to placing it in your path

source  $HOME/.bashrc
```

2. Clone this repo

```
git clone https://github.com/kramundson/tetraploid_kmers
cd tetraploid_kmers
```

3. Build environment using kmer.yaml

```
conda env create --name kmer -f kmer.yaml
```