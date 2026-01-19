# Assembly and binning !

The aim of this tutorial is to go from reads to Mags placed in a phylogenetic tree 

The workflow is quite typical and involve

1. [Assembly](#assembly)

2. [Read mapping](#readmapping)

3. [Contig binning](#binning)

4. [Bin quality ](#checkm)

5. [Phylogenetic tree placement](#gtdb)
 


We are now going to perform a basic assembly based metagenomics analysis of these same samples. 
This will involve a collection of different software programs:

1. [MetaMDBG](https://github.com/GaetanBenoitDev/metaMDBG): A highly efficient metagenomics assembler currently our default for most studies

2. [minimap2](https://github.com/lh3/minimap2): Necessary for mapping reads onto contigs

3. [samtools](http://www.htslib.org/download/): Utilities for processing mapped files

4. [Metabat2](https://github.com/BinPro/CONCOCT): an automatic binning algorithm
5. [checkm](https://ecogenomics.github.io/CheckM/#:~:text=CheckM%20provides%20a%20set%20of,copy%20within%20a%20phylogenetic%20lineage.): Tools to assess bin quality
6. [gtdb-tk](https://github.com/Ecogenomics/GTDBTk): Toolkit to place MAG on reference phylogenetic tree and use placement for taxonomy classification. 



## Conda Envs
We use a [conda](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf) env to install all dependencies, you don't need to install anything and all dependencies are available but only inside that environment.   

Try to check the command line help of megahit

    metaMDBG
<details><summary>not working?</summary>
<p>
Conda environment are created as independant environment to everything else, you need to "activate" an environment to be able to access the sets of tools installed inside.

    conda env list
    conda activate MAGs
    metaMDBG

</p>
</details>


## Assembly
<a name="assembly"/>

Let's create a Projects directory and work inside:

    mkdir -p Projects/MAGs


The dataset for this tutorial can be found at: 

    ~/Documents/databases/zymo_gut_mock

Or equivalently here:

    /home/train/Documents/databases/zymo_gut_mock
    

We are going to use metaMDBG for assembly. It is a fast memory efficient metagenomic assembler and particularly usefull for handling large coassembly.

Bioinformatic is mostly about reading documentation, and looking on the internet how to do things. 
Use the -h flag on metaMDBG and try to craft a command line to launch the assembly on 8 cores.

<details><summary>spoiler (click on me if you're desperate)</summary>
<p>

```bash
cd ~/Projects/MAGs
metaMDBG asm --out-dir Assembly --in-ont /home/train/Documents/databases/zymo_gut_mock/sub_sample_1.fastq.gz --threads 16
```
It should take around 800 seconds
</p>
</details><br />


What is the output of an assembler?<br />
What is a coassembly and why would you want to do one?<br />
How good is the assembly? <br />
What are chimeric assemblies<br />
How would we go for estimating the number of organisms in the assembly?<br />
What do we do with an assembly?<br />


<a name="readmapping"/>

## Read mapping

What informations can be used to bins contigs?

We use minimap2 to map reads to the assembly.

Then I want you to try mapping the reads minimap2 inside a subdirectory Map to produce a sam file 'Map/sample1.sam':

<details><summary>the correct command is:</summary>
<p>

```bash
cd ~/Projects/MAGs
mkdir Map
minimap2 -ax lr:hq --sam-hit-only -t 16 Assembly/contigs.fasta.gz /home/train/Documents/databases/zymo_gut_mock/sub_sample_1.fastq.gz > Map/sub_sample_1.sam
```
</p>
</details>

You can look at the sam:
```
tail Map/sub_sample_1.sam
```

It is quite a complex [format](https://en.wikipedia.org/wiki/SAM_(file_format))

The sam file is a bit bulky so we never store alignments in this format instead we would convert it into bam. Can you convert this file using 
'samtools view':


<details><summary> Convert sam to bam command</summary>
<p>

```bash
cd ~/Projects/MAGs/Map
samtools view -h -b -S sub_sample_1.sam > sub_sample_1.bam
```
</p>
</details>

For downstream analysis we needs the bam file to be sorted:
```
samtools sort sub_sample_1.bam -o sub_sample_1.sorted.bam 
```


#### OPTIONAL 
Using samtool we can filter only those reads which are mapped to the assembly.
```bash
samtools view -b -F 4 sub_sample_1.bam > sub_sample_1.mapped.bam
```


When you have multiple samples and you want to map all of them in turn, here is a loop written in bash:

```bash
cd ~/Projects/MAGs

for file in /home/train/Documents/databases/zymo_gut_mock/prev/sub_sample_*.fastq.gz
do 
   
   stub=${file%.fastq.gz}
   name=${stub##*/}
   
   echo $name

  minimap2 -ax lr:hq --sam-hit-only -t 16 assembly/contigs.fasta.gz $file | samtool view -h -b -S - | samtools sort - > Map/$name.sorted.bam
done

```

<a name="binning"/>

## Contig binning

The first step is to derive coverage from bam files. For this we can use metabat2 script. It takes bam files as inpute produce a table of mean coverage depth and std for each contigs in each sample.

```bash
cd ~/Projects/MAGs/Map
jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
```

Make a new subfolder Binning. Move the Coverage file into this and look into crafting the metabat2 command line. Either use the help flag or a quick google search.

<details><summary> Solution</summary>
<p>

```bash
cd ~/Projects/MAGs
mkdir Binning
mv Map/depth.txt Binning/depth.txt
metabat2 -i Assembly/contigs.fasta.gz -a Binning/depth.txt -t 16 -o Binning/Bins/Bin
```
</p>
</details>

How many contigs were clustered? 
```bash
grep -c ">" Binning/Bins/*.fa | awk -F: '{ s+=$2 } END { print s }'
```
How many nucleotide were clustered?
```bash
grep -v ">" Binning/Bins/*.fa | cut -f2 -d":" |wc -m
```
 
## Which bins are Metagenome assembled genomes (MAGs)?

A bin is a group of contigs put together from looking at coverage/composition. How do you assess bin quality?

Checkm2 is an handy automated pipeline which will ML applid on presence absence of genes database to assess contamination/completion.
```bash
CHECKM2_DB=/home/train/Documents/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd
cd ~/Projects/MAGs/Binning
checkm2 predict --threads 16 --database_path $CHECKM2_DB --input Bins -x .fa --force --output-directory checkm2
```

<details><summary> Let's look at the result which bins can be considered MAGs? </summary>
<p>
"near-complete" if its completeness is ≥90% and contamination is ≤5%<br />
"high quality" if its completeness is ≥70% and contamination is ≤10%<br />
"medium quality" if its completeness is ≥50% and contamination is ≤10%<br />
</p>
</details>

Using awk we can select bins from the results of checkm2:
```bash 
awk -F'\t' '($2 > 50) && ($3 < 10)' */quality_report.tsv
```
The return of the bash loop
```bash 
for f in $(awk -F'\t' '($2 > 50) && ($3 < 10)' */quality_report.tsv | cut -f 1)
do
    cp ${f}.fa mags/
done
```

## A better idea of MAG identity

When doing metagenomic, it happens often that the MAGs you obtain are not in database, or just distantly related to references. Using single copy core genes we can place these MAGs in a phylogenetic tree full with known species representatives.

The gtdb toolkit does that for you:

```bash
GTDBTK_DATA_PATH= <PATH> # change the db path to the correct dir 
cd ~/Projects/MAGs/Binning
mkdir -p gtdb/scratch
gtdbtk classify_wf --cpus 16 --genome_dir mags --out_dir gtdb --extension .fa --scratch_dir gtdb/scratch --skip_ani_screen
```
That will take at least, 30 min min, instead lets use the pre-processed results

```bash
cd ~/Projects/MAGs/Binning
rm -r gtdb
ln -s /home/train/Documents/mag_prerun/MAGs/binning/gtdb .
```

We obtain multiple files what are they?
Look at summary files what information can we obtain.
What is the RED, from gtdb?

![alt tag](https://github.com/Sebastien-Raguideau/Ebame/blob/master/Figures/gtdb.jpg)

