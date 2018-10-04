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

For more reading on repeat diversity within cultivated potato, see the references located
under ```papers```

Instructions for downloading datasets and raw data are in the subdirectory ```download/```
Each experiment has its own subdirectory under ```experiments/```
Analysis of informative or interesting experiments are migrated to ```results/```

1. Download and install miniconda for Python3

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# install
bashsh Miniconda3-latest-Linux-x86_64.sh

# review license
# accept license
# If you don't have a home folder, change home location to be your folder in shared disk
# space. Exapmle for user jdoe:
/share/comailab/jdoe/miniconda3

# Do not prepend miniconda3 to path. Instead, do the following:
echo 'export PATH="share/comailab/jdoe/miniconda3/bin:$PATH"' > /share/comailab/jdoe/.bashrc
source /share/comailab/jdoe/.bashrc
```

2. Clone this repository and navigate to it

```
git clone https://github.com/kramundson/potato_tandem_repeats
cd potato_tandem_repeats
```

3. If you're running this as a visting scholar or rotation student, chances are you don't
have your own home folder on the cluster. In this case, activate the conda environment
maintained for this project:

```
source activate /share/comailab/kramundson/miniconda3/envs/kmer/
```

If this ran successfully, your command prompt should look different now:

[alt_text]("ImageURL")

4. If you're not running this from a UCD server, you'll have to build the conda environment
yourself, and then activate it.

This can be done using the included ```kmer.yaml``` file. This will work smoothly, provided
you have a home directory you can write things to.

```
# build
conda env create --name kmer -f kmer.yaml

# activate
source activate kmer
```

```