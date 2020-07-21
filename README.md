# SCCNV
SCCNV: a software tool for identifying copy number variation from single-cell whole-genome sequencing

Version 1.0.2

Updated date: 2020.07.20

#####
## Author and License

Author: Xiao Dong

Email: biosinodx@gmail.com, xiao.dong@einsteinmed.org

Licensed under the GNU Affero General Public License version 3 or later

#####
## DEPENDENCIES

python 2 or python 3

python modules os, argparse, math, numpy

bedtools, samtools

#####
## OPTIONS

-i, (or --ibam): a text file providing a list of bam files of single-cell whole-genome sequencing data; one row one file; required input file.

-o, (or --odir): a directory to output result; required output directory

-g, (or --genome): genome version; default: hg19

-s, (or --genomesize): genome size; default: 3140000000

-b, (or --binsize): binsize (bp); default: 500000; min: 1000

-w, (or --windowsize): windowsize (number of bins); default: 10

-p, (or --ploidy): ploidy of autosomes; default: 2

-m, (or --minmap): minimum mappability; default: 0.3

-n, (or --maxN): maximum N in the bin; default: 0.1

-k, (or --skipbed): skip calulcating bed with bedtools; default: False

-r, (or --report): report intermediate results; default: False

#####
## USAGE

### I. Prepare bam files of every single cell for analysis

Below is an example for cell_A. I recommend filter out reads with mapq<30.

samtools view -b -q 30 ./bam/cell_A.bam > ./bam_mapq30/cell_A.mapq30.bam

samtools index ./bam_mapq30/cell_A.mapq30.bam

### II. Prepare a list of the single-cell bam files.

Provide a file (e.g. “bamlist.txt”) with the following content. See NOTE below for important information when analyzing multiple single cells of a same clone, e.g., tumor cells.

./bam_mapq30/cell_A.mapq30.bam

./bam_mapq30/cell_B.mapq30.bam

./bam_mapq30/cell_C.mapq30.bam

### III. Perform CNV calling

python sccnv.py -i bamlist.txt -o cellsAtoC -k False -r True

### IV. Tips

One can prepare bed files for calculation using bedtools and samtools (for every cell in paralleles),

mkdir cellsAtoC

bedtools makewindows -g ./resource/hg19.chrlength.txt -w 500000 > ./cellsAtoC/hg19.bin500000.bed

samtools bedcov ./cellsAtoC/hg19.bin500000.bed ../bam_mapq30/cell_A.mapq30.bam > ./cellsAtoC/cell_A.depth.bin500000.bed

And skip the related step in SCCNV,

python sccnv.py -i bamlist.txt -o cellsAtoC -k True -r True

### V. Example data

An example dataset is included in "sccnv_example_v1.0.2.zip". It includes (1) intermediate and final result files generated using SCCNV in txt format; and (2) basic R scripts for visualizing CNVs across the genome.

#####
## NOTE

(1) SCCNV aims to discover difference in CNV between every single cell in the bamlist.txt and the other cells in the bamlist.txt. When analyzing CNV of multiple tumor cells, it is not appropriate to include all tumor cells in the bamlist.txt. Instead, please use one tumor cell with two or more normal diploid cells in the bamlist.txt.

(2) SCCNV only used 0-4 because the final copy number call was after multiple testing correction and I wished to minimize the numbers of hypothesis tested, i.e., 5 for 0-4. However, this will result in an under estimation if the real copy number exceed 4. Please check SCCNV intermediate files to make sure not exceeding 4.

#####
## RELEASE NOTES

v1.0.2 2020.07.20, allowed calling CN>4 and included additional example data.

v1.0.1 2020.04.10, revised readme and included example data.

v1.0.0 2018.07.26, release version.

v0.0.3 2018.07.25, rewrote in Python.

v0.0.2 2017.08, fixed bugs.

v0.0.1 2016.09, drafted script in R.
