<img src="http://aca.unl.edu/static/images/help_page_pngs/AOminer-logos.png" width="600" height="400">

<center>(c) <a href='http://bcb.unl.edu'>Yin Lab</a>@<a href='https://www.unl.edu'>UNL</a>2022</center>


## Contents:

<a href='#installation'>I. Installation / Dependencies</a>

<a href='#about'>II. About</a>

<a href='#using_acafinder'>III. Using AOminer</a>


<a href='#examples'>V. Examples</a>

<a href='#workflow'>VI. Workflow</a>

<a href='#faq'>VII. FAQ</a>

****

<div id='installation' />

## I. Installation / Dependencies

Git clone the Github directory to your own, and go to the AOminer directory

### Dependences 

Program expects these versions and using other versions can result in unexpected behavior.

`VIBRANT` - Used to search for potential prophage regions from input genomic sequences

Version used v1.2.0. We recommend installing VIBRANT using Anaconda, but you may also install VIBRANT with other methods from https://github.com/AnantharamanLab/VIBRANT

To install using Anaconda, 
Install dependencies. See Requirements section https://github.com/AnantharamanLab/VIBRANT.
Install directly to $PATH using bioconda. 
```sh
conda install -c bioconda vibrant==1.2.0
```

Download and setup databases. This will take some time due to file sizes, but it only needs to be run once. This step requires ~20GB of temporary storage space and ~11GB of final storage space. To do this, run download-db.sh which should be in your system's $PATH. download-db.sh

`Cctyper` - Used for complete CRISPR-Cas system search

Can be installed either through conda or pip.
It is advised to use conda, since this installs CRISPRCasTyper and all dependencies, and downloads the database in one go.

To use conda:
```sh
conda create -n cctyper -c conda-forge -c bioconda -c russel88 cctyper
```

To use pip:
python -m pip install cctyper

When installing with pip, you need to download the database manually:

```sh
# Download and unpack
svn checkout https://github.com/Russel88/CRISPRCasTyper/trunk/data
tar -xvzf data/Profiles.tar.gz
mv Profiles/ data/
rm data/Profiles.tar.gz

# Tell CRISPRCasTyper where the data is:
# either by setting an environment variable (has to be done for each terminal session, or added to .bashrc):
export CCTYPER_DB="/path/to/data/"
# or by using the --db argument each time you run CRISPRCasTyper:
cctyper input.fa output --db /path/to/data/
```

Detailed information can be found at : https://github.com/Russel88/CRISPRCasTyper#install

`PfamScan/Pfam Database` - For protein annotation 

Install PfamScan with Anaconda:
```sh
conda install -c bioconda pfam_scan
#or
conda install -c bioconda/label/cf201901 pfam_scan
```

Downloading Pfam database
```sh
# a. Go to the all_pFam_hmm/ folder

# b. Download database: 
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

#c. Unpack 
gunzip Pfam-A.hmm.gz

#d. Prepare HMM files 
hmmpress Pfam-A.hmm
```

****


<div id='about' />

## II. About 

### AOminer is the first ever automated tool for reliable screening of anti-CRISPR operons (AOs) using machine learning HMM algorithms 

 Using the built-in two-state HMM trained based on known AOs, AOminer is able to screen for potential AOs from both genomic data and individual operons. AOminer holds an accuracy of 89% from our performance tests, and outperformed all other machine learning based anti-CRISPR prediction tools in terms of AO prediction. In addition, AOminer also searches the input genome/operon for Acr homologs using published Acr proteins, potential Aca proteins containing HTH domains in predicted AOs, CRISPR-Cas systems using (CRISPRCasTyper) and an in-house Self-targeting spacer (STS) search tool, prophage regions using (VIBRANT). All information will be integrated with the predicted AOs, providing users with detailed information vital to prediction assessment.

****

<div id='using_acafinder' />

## **III. <span style='color:RebeccaPurple'>Using AOminer</span>**

#### Input 

AOminer needs **.fna**, **.gff** and **.faa** as input. Only **.fna** file as input is also acceptable; in that case, the **.gff** and **.faa** file will be generated by running <a href='https://github.com/hyattpd/Prodigal'>Prodigal</a>. If you have an operon/gene cluster of interest, and wish to see its probability of being AO, AOminer need only the **.faa** of the operon as input (**NOTE:** The input operon or gene cluster could be any sequence or combinations of proteins, and do not need to follow the parameters of a short-gene-operon). 

#### List of Options

| Option | Alternative | Purpose                     |
| -----  | ----------- | --------------------------- |
| -h     | --help      | Shows all available options |
| -n     | --FNA_file  | <span style="color:red">Required</span> fna file |
| -g     | --GFF_file     | <span style="color:red">Required</span> Path to gff file to use/parse |
| -p     | --FAA_file     | <span style="color:red">Required</span> Path to faa file to use/parse |
| -m     | --mode_prodiagal  | Mode prodigal will be run choices=["single","meta"], default=meta |
| -o     | --outputFolder      | Folder containing all output results, default=AcaFinder_Output |
| -j     | --just_operon      | Provide flag if input data is the faa file of "ONE" operon/gene cluster |
| -l     | --all_protein_length_in_AcrAca_operon    | Max proten length in Acr-Aca operon when length of Acr homolog < 200aa, default=600 |
| -i     | --intergenic_dist_in_AcrAca_operon      | Maximum Intergenic distance in Acr-Aca operon, default=250 |
| -r     | --Acr_protein_database     | The Acr proteins that will be used search for Acas, default are the published Acrs + AcrHub predicted Acrs + 2500 high confident Acr prediction of AcrCatalog, default=AcrDatabase.faa |
| -w     | --Prok     | Provide option -w/--Prok if input data is a Prokaryotic assembled genome/contig |
| -d     | --threads     | Number of cpu cores to run the program with, default=1 |
| -z     | --phamDir     | Directory of all pfam hmm files with .dat files and other binaries, default=AcaFinder/all_pFam_hmm | 
| -a     | --HTH_HMM     | Curated HTH databases, recommended to use the default hmm provided from us, default=AOminer/HTH_hmms/HTH_HMM_strict |
| -e     | --PfamScan_evalue     | Evalue cut-off for Pfam annotation, recommended to use default, default=1e-1 |
| -x     | --execute_level     | Selection levels to predicted AOs, recommeded strict, choices=["strict","medium","relaxed"] |


#### Output files

| Name                 | Meaning     |
| -------------------- | ----------- |
|*<output_dir>*/Short_Gene_Operons   | Folder containing intermediate short gene operon result files |
|*<output_dir>*/CRISPR_Cas_Found   | CCtyper direct output folder |
|*<output_dir>*/VIBRANT_*<input_ID>*_genomic   | VIBRANT direct output folder |
|*<output_dir>*<input_ID>*.prodigalOUT   | Prodigal annotation folder of genomic input |
|*<output_dir>*/ALL_AOs_predicted.csv    | Final results from AOminer's AO prediction, generated by scanning genome/contig input |
|*<output_dir>*/AO_predicted.csv    | Final results from AOminer's AO prediction, generated by scanning ONE operon input |
|*<output_dir>*/Acr_homologs.faa    | Protein seuqnces of Acr homologs found from input genomic sequences | 
|*<output_dir>*/CRISPR-Cas_found.csv   | Summary of Complete CRISPR-Cas systems discovered |
|*<output_dir>*/Homology_Search_Results   | Output folder from Acr homology & HTH domian search |
|*<output_dir>*/prophage_locations.csv    | Summary of prophage regions discovered |
|*<output_dir>*/Short_Gene_Operons/SGO_OperonNumber-*<operon_ID>*.faa   | Protein fasta file of short gene operon | 
|*<output_dir>*/Short_Gene_Operons/Homology_Search_Results/SGO_OperonNumber-*<operon_ID>*.faa.hmmout.| hmmscan output of SGO annotation with dbPFhmm |
|*<output_dir>*/Short_Gene_Operons/Homology_Search_Results/SGO_OperonNumber-*<operon_ID>*.faa.hmmout..Coverage_parsed | Coverage parsed hmmscan output of SGO annotation with dbPFhmm |
*<output_dir>*/Short_Gene_Operons/SGO_OperonNumber-*<operon_ID>*.faa.pfamScanOut | Pfam annoatations of proteins of short-Gene-Operons |

****

<div id='examples' />

## **IV. <span style='color:RebeccaPurple'>Examples</span>**

```sh
python3 AOminer_runner.py --FNA_file sample_organisms/GCA_002194095.1_ASM219409v1_genomic.fna --FAA_file sample_organisms/GCA_002194095.1_ASM219409v1_genomic.gff --GFF_file sample_organisms/GCA_002194095.1_ASM219409v1_protein.faa -o [output_dir] 
```
or you can only use **.fna** file as input.

```sh
python3 AOminer_runner.py --FNA_file sample_organisms/GCA_002194095.1_ASM219409v1_genomic.fna -o [output_dir] 
```
If you wish only to input an operon
```sh
python3 AOminer_runner.py --FNA_file sample_organisms/Operon_exp.faa --just_operon -o [output_dir] 
```

You will see the output result in output_dir/. If you dont specifiy an output_dir, result will be in AcaFinder_Output/


****

<div id="workflow" />

## **V. <span style='color:RebeccaPurple'>Workflow of AOminer</span>**

<img src="http://aca.unl.edu/static/images/help_page_pngs/AOminer_Pipeline.png">

With provided input, AOminer will proceed with the following pipeline:

Step 1) The input FAA will be scanned for short-gen operons (SGOs) with the following criteria: (i) All genes < 200aa (if Acr homolog present in SOG and homolog is >200aa, then all genes < Acr homolog length); (ii) All intergenic distances < 250bp; (iii) All genes on the same strand. If input is a single operon, then this step is bypassed.
Step 2) The extracted SGOs or input operon will be annotated with the built-in dbPFhmm, of which contain HMMs of protein families found in abundance within known AOs.
Step 3) The annotated operons will go through the built-in two-state HMM. 
In addition to the identification of anti-CRISPR operons, AOminer will also scan the input fna file for prophages, CRISPR-Cas systems & self-targeting spacers (STSs), Acr homologs, and potential Acas. If input is a single operon, only Acr homologs and potetnial Aca proteins will be searched.

**All generated information will be associated together, and provided to the users as tables.** 

****
