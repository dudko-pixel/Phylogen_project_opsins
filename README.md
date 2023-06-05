# Phylogen_project_opsins
Here's a repo with some code for the free interpretation of the "The diversity of opsins in Lake Baikal amphipods (Amphipoda: Gammaridae)" article

## Aim: 
The aim of this project was to reproduce results from the article. Basically the target was to find out which phylogenetic groups had opsins and how diverse their quantity was.

## Tasks:
To achieve the mentioned aim it was decided to set some necessary tasks:
* Master the techniques used - scroll through the scripts used, explore the new programs
* Beautifully arrange a repo with all methods used - make a fork of a repo with possibly some custom additions
* Estimate the reproducibility of methods described - repeat the analysis from the article (or a part of it) and compare the results
* Add one more database to use - a supplementary task that appeared while reproducing the analysis

## Results:
All tasks have been completed. After running PIA3 on the downloaded data there were some mismatches in the quantity of opsins between the published and the newly obtained data (however, mostly on the high side). While checking through the pipeline for any commands that went missing in the current launch, we discovered that the file with hits not only contained opsins but also other G protein-coupled receptors. To avoid running the whole pipeline again, we wrote an additional script to filter out excessive sequences based on the script from the article. Also, since some other databases of non-opsins vision-related genes were not included in the PIA3 after evolving from PIA2, they were included by concatenating multiple fasta files with them to multi-fasta files and adding them to the PIA3.
The diversity of opsins in amphipods was re-explored. At first, it seemed that the results were not reproduced but after discovering the mistake of not specifying the  –opsin option it was fixed. The tool was also replenished with additional databases from the previous version to increase its flexibility.
The main result of this work is an Excel table with hits of opsins count, it can be found in a **results** directory

## Literature:
Drozdova, P., Kizenko, A., Saranchina, A. et al. The diversity of opsins in Lake Baikal amphipods (Amphipoda: Gammaridae). BMC Ecol Evo 21, 81 (2021). https://doi.org/10.1186/s12862-021-01806-9

## PIA3

Modified from PIA2 (https://github.com/xibalbanus/PIA2).
Only tested in Ubuntu-based Linux systems but should work anywhere else with minimum adjustments in the installation process.

### Installation

This pipeline requires `Conda`. If it is not installed on your computer, you need to:

**1.** Download [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) (look [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda) to decide which is right for you).

**2.** Run bash installation script like.
```commandline
bash Anaconda3-2020.02-Linux-x86_64.sh
``` 
or 
```commandline
bash Miniconda3-latest-Linux-x86_64.sh
```

**3.** Source bashrc to activate conda.

```commandline
source ~/.bashrc
```

**4.**
`Mamba` and `Snakemake` are also required for pipeline running.

```commandline
conda install -n base -c conda-forge mamba
conda create -n snakemake python=3.10.8
conda install -c conda-forge -c bioconda snakemake=7.25.0
```
Please be patient, if solving environment takes some time.

Activate *snakemake* conda environment.
```commandline
conda activate snakemake
```

**5.**
Download and unpack or clone PIA3 directory on your computer. PIA3 is installed as a separate conda environment during its first run. We recommend to run it on our test data.

```commandline
git clone https://github.com/AlenaKizenko/pia3_amphipod_opsins.git
cd pia3_amphipod_opsins/PIA3
```

Replace PATH/TO/PIA_ENV with preferred **full** path for your PIA3 environment.
```commandline
snakemake --cores 8 --use-conda --conda-prefix PATH/TO/PIA_ENV --conda-frontend conda --configfile config.yaml 
```
**6.** Run unit test.

```commandline
 python3 -m unittest test_PIA3.py
```

### Run the pipeline
For PIA3 run, you need to provide config.yaml file, in which all the requied parameters are stored. Please use only absolute paths to avoid an error.

Explanation of config-file fields:

* `in_dir`: path to folder with input reference `.fasta` file(s) **required**

* `out_dir`: output directory path **required**

* `db`: path to database **required** 
  (Here we offer you to specify one of 13 available databases that we added specifically from the PIA2 version)

* `transcripts`: **cds** - perform BLAST search only on coding sequences (longer than 1/2 of mean database sequence and starting with methionine); **all** use all sequences

* `clean`: delete intermediate files **default True**

* `model`: model for IQ-Tree maximum likelihood tree building (put **TEST** if not known)

* `outgroup`: outgroup for phylogenetic tree building; if not defined by user, first sequence from database FASTA file is taken

* `cd_h`: CH-HIT clustering treshhold (if 1 is chosen, CH-HIT clusters only identical sequences)

* `opsins`: check for retinal lysine in found sequences and output only opsins; NOTE: **RHO_Bos_taurus_AAA30674.1** sequence needs to be in the database file (can be found in `classification_opsins_full_aa.fasta`)

**7.**
In case you forgot to specify the **--opsin** option, there is a script filter_distant_seqs_add.py in this directory to filter your hits from excessive results.
