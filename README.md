# SCCNV
SCCNV: a software tool for identifying copy number variation from single-cell whole-genome sequencing

Version 1.0

Updated date: 2017.07.27

#####
## Author and License

Author: Xiao Dong

Email: biosinodx@gmail.com, xiao.dong@einstein.yu.edu

Licensed under the GNU Affero General Public License version 3 or later

#####
## DEPENDENCIES

python 2 or python 3

python modules os, argparse, math, numpy

bedtools

#####
## USAGE

### I. Prepare bam files of every single cell for analysis

Below is an example for cell_A. I recommend filter out reads with mapq<30.

samtools view -b -q 30 ./bam/cell_A.bam > ./bam_mapq30/cell_A.mapq30.bam

samtools index ./bam_mapq30/cell_A.mapq30.bam

### II. Prepare a list of the single-cell bam files.

Provide a file (e.g. “bamlist.txt”) with the following content,
./bam_mapq30/cell_A.mapq30.bam
./bam_mapq30/cell_B.mapq30.bam
./bam_mapq30/cell_C.mapq30.bam

### III. Perform CNV calling

python sccnv.py -i bamlist.txt -o cellsAtoC -k False -r True

### IV. Tips

One can prepare bed files for calculation using bedtools and samtools (for every cell in paralleles),

mkdir cellsAtoC

bedtools makewindows -g ./resource/hg19.chrlength.txt -w 500000 ./cellsAtoC/hg19.bin500000.bed

samtools ./cellsAtoC/hg19.bin500000.bed ../bam_mapq30/cell_A.mapq30.bam > ./cellsAtoC/cell_A.depth.bin500000.bed

And skip the related step in SCCNV,

python sccnv.py -i bamlist.txt -o cellsAtoC -k True -r True

#####
## RELEASE NOTES

v1.0 2018.07.26, release version.

v0.0.3 2018.07.25, rewrote in Python.

v0.0.2 2017.08, fixed bugs.

v0.0.1 2016.09, drafted script in R.
