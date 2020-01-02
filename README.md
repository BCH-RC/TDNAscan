# Installing TDNAscan 

## Operating System Requirements

TDNAscan has been tested on the following Linux distributions:

* Ubuntu 14.04 LTS
* Ubuntu 16.04 LTS
* CentOS 7.3
* Debian 7 "Wheezy"
* Debian 8 "Jessie"

## Software Dependencies
The following programs need to be installed and the executable commands should be in $PATH of system.
* BWA (Version ="0.7.12")
* Samtools (Version ="1.3.1")
* Python (Version ="2.7.5")


# Using TDNAscan
* Recommendation: trim your NGS reads using Trimmomatics or other NGS trimmers before using TDNAscan, otherwise, the final results will include some false positve insertions.

## Step 1 - Identify complete and truncated T-DNA insertions
  
### Usage: 

`python tdnascan.py -1 forward.fq -2 reverse.fq -t t-dna.fa -g ref_genome.fa -p tdna`

### Parameters:

* REQUIRED -1 the paired read file 1
* REQUIRED -2 the paired read file 2
* REQUIRED -t the T-DNA sequence file in fasta format
* REQUIRED -g the genome sequence file in fasta format
* REQUIRED -p the name of your project (output files will be placed in a directory with the name you provide)
* OPTIONAL:
	* -@ cpu number for BWA and SAMTOOLS [default 8]
	* -a the window size of clustering soft clipped reads [default:3]
	* -b the length of library fragment in NGS data [default:500]

### Example:

An example data set is provided with this repository.

Running the following example code will create a project directory named 'tdna' relative to where you run the command, and will produce example output:

`python tdnascan.py -1 mt4_chr1_20x_mut_tdna_1.fq -2 mt4_chr1_20x_mut_tdna_2.fq -t t-dna_elison.fa -g mt4_chr1_2Mb.fa -p tdna`

## Step 2 - Annotate complete and truncated T-DNA insertions

### Usage: 

`python tdnaAnnot.py -i tdna_insertion.bed -f ref.gff3 -o tdna_insertion_annot.bed`

### Parameters:

* REQUIRED -i T-DNA BED file
* REQUIRED -f gff3 annotation file
* REQUIRED -o annotated insertion file

### Example:

`python tdnaAnnot.py -i tdna40x/tdna40x_insertionBED.txt -f Athaliana_447_Araport11.gene.gff3 -o ./tdna/tdna_insertion_annot.bed`


### Output:

### BED file
TDNAscan produces a single BED file which contains all unique deletions that were identified.

The output is placed in ./tdna (i.e. in the directory named after your project)

Running the above example code would produce the following BED file:

* ./tdna/**tdna_insertion.bed**

Annotated BED file:

* ./tdna/**tdna_insertion_annot.bed**

### Output file structure

* Chr: Chromosome number;	
* Breakpoint: Start position of insertions;
* SuppRead: CLR represents the clipped reads number; DIR represents discordant reads number;
* SuppRead: tdna_st and tdna_end represent the start and end position of T-DNA sequence truncated when inserted to reference genome.
* Orientation: forward or reverse T-DNA inserted to reference genome;
* Freq: Insertion frequency;
* Genes (optional): This column will only show genes if deletions cover.


# Contact

* Dr. Liang Sun    (lsun@noble.org)
* Yinbing  Ge  (yinge@noble.org)
