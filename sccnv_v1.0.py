### SCCNV: a software tool for identifying copy number variation from single-cell whole-genome sequencing
# Copyright (C) 2018  Dong, et. al.
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.

### Updates
# v1.0 2018.07.26, release version
# v0.0.3 2018.07.25, rewrote in Python.
# v0.0.2 2017.08, fixed bugs
# v0.0.1 2016.09, drafted script in R.

# require: linux or unix; samtools; bedtools
import argparse, os, numpy as np
from math import *

parser=argparse.ArgumentParser(description="sccnv.py, v1.0, Xiao Dong, biosinodx@gmail.com, xiao.dong@einstein.yu.edu")
parser.add_argument("-i", "--ibam", type=str, required=True, help="INPUT: a text file providing a list of bam files of single-cell whole-genome sequencing data; one row one file; required")
parser.add_argument("-o", "--odir", type=str, required=True, help="OUTPUT: a directory to output result; required")

parser.add_argument("-g", "--genome", type=str, default="hg19", help="INPUT: genome version; default: hg19")
parser.add_argument("-s", "--genomesize", type=int, default=3140000000, help="INPUT: genome size; default: 3140000000")
parser.add_argument("-b", "--binsize", type=int, default=500000, help="INPUT: binsize (bp); default: 500000; min: 1000")
parser.add_argument("-w", "--windowsize", type=int, default=10, help="INPUT: windowsize (number of bins); default: 10")
parser.add_argument("-p", "--ploidy", type=int, default=2, help="INPUT: ploidy of autosomes; default: 2")
parser.add_argument("-m", "--minmap", type=float, default=0.3, help="INPUT: minimum mapability; default: 0.3")
parser.add_argument("-n", "--maxN", type=float, default=0.1, help="INPUT: maximum N in the bin; default: 0.1")

parser.add_argument("-k", "--skipbed", type=str, default=False, help="INPUT: skip calulcating bed with bedtools; default: False")
parser.add_argument("-r", "--report", type=str, default=False, help="INPUT: report intermediate results; default: False")
args=parser.parse_args()
print(args)

os.system("mkdir " + args.odir)

### 1. samtools / bedtools to calculate genome fragments ###
cell=[]
infile=open(args.ibam)
# infile=open('bamlist_BCELL.txt')

lines=infile.readlines()
infile.close()
for line in lines:
	if len(line)==0:
		break
	line=line.split('/')
	cell.append((line[len(line)-1]).split('.bam\n')[0])
print(cell)

if not(args.skipbed):
	line="bedtools makewindows -g " + "./resource/"+ args.genome + ".chrlength.txt" + " -w " + str(args.binsize) + " > " + args.odir + "/" + args.genome + ".bin" + str(args.binsize) + ".bed"
	print(line)
	os.system(line)
	for i in range(0, len(cell)):
		x="samtools bedcov " + args.odir + "/" + args.genome+ ".bin" + str(args.binsize) + ".bed " + lines[i].split('\n')[0] + " > " + args.odir + "/" + cell[i] + ".depth.bin"+ str(args.binsize) +".bed"
		print(x); os.system(x)

### Input all data ###
# def readbed(cellname, odir='./tmp', binsize=500000):
def readbed(cellname, odir=args.odir, binsize=args.binsize):
	infile=odir + '/' + cellname + '.depth.bin' + str(binsize) + '.bed'
	infile=open(infile)
	lines=infile.readlines()
	infile.close()
	cc=[];p1=[];p2=[];dd=[]
	for line in lines:
		if len(line)==0:
			break
		line = line.split('\n')[0]
		line = line.split('\t')
		cc.append(line[0]); p1.append(line[1]); p2.append(line[2]); dd.append(line[3])
	return([cc,p1,p2,dd])

dat0=[]
for i in cell:
	tmp = readbed(i)
	if len(dat0)==0:
		dat0.append(tmp[0])
		dat0.append(tmp[1])
		dat0.append(tmp[2])
	dat0.append(tmp[3])

# print(dat0)

##### For future, automatically prepare mapability scores and GC content file #####
# temp = './resource/' + 'hg19' + '.mapability_gc.bin' + str(500000) + '.bed'
temp = './resource/' + args.genome + '.mapability_gc.bin' + str(args.binsize) + '.bed'
normal_file=open(temp)

nn = normal_file.read()
normal_file.close()
nn = nn.split('\n')
t0=[];t1=[];t2=[];t3=[];t4=[];t5=[];t6=[]
for i in nn[1:]:
	if len(i)>0:
		temp = i.split('\t')
		t0.append(temp[0])
		t1.append(temp[1])
		t2.append(temp[2])
		t3.append(temp[3])
		t4.append(temp[4])
		t5.append(temp[5])
		t6.append(temp[6])

nn = [t0,t1,t2,t3,t4,t5,t6]

# define output function
def output(datX, ofile_root):
	ofile=open(ofile_root, 'w')
	ll = 'chr\tpos1\tpos2\tMapability\tAT\tGC\tN'
	for i in cell:
		ll = ll + '\t' + i
	ll = ll + '\n'
	ofile.write(ll)
	for i in range(0, len(datX)):
		ll=datX[i][0]
		for j in datX[i][1:]:
			ll = ll + '\t' + str(j)
		ll = ll + '\n'
		ofile.write(ll)
	ofile.close()

# Merge data file and normalization file #
dat1=[]
for i in range(0, len(dat0[0])):
	for j in range(0, len(nn[0])):
		if dat0[0][i]==nn[0][j] and int(int(dat0[1][i])/10)==int(int(nn[1][j])/10) and int(int(dat0[2][i])/10)==int(int(nn[2][j])/10):
			temp = [dat0[0][i], int(dat0[1][i]), int(dat0[2][i]), float(nn[3][j]), float(nn[4][j]), float(nn[5][j]), float(nn[6][j])]
			for k in range(3, len(dat0)):
				temp.append(float(dat0[k][i]))
			dat1.append(temp)
			break

if args.report:
	output(dat1, args.odir+'/result.dat1_raw.txt')

### Normalize by mapability: dat2 ###
dat2=[]
for i in range(0, len(dat1)):
	# if dat1[i][3]>=0.3 and dat1[i][6]<=0.1:
	if dat1[i][3]>=args.minmap and dat1[i][6]<=args.maxN:
		temp = dat1[i][0:7]
		for j in range(7, len(dat1[0])):
			temp.append(dat1[i][j]/dat1[i][3])
		dat2.append(temp)

if args.report:
	output(dat2, args.odir+'/result.dat2_mapability.txt')

### Normalize by GC content: dat3 ###
gc=[]
for i in range(0, len(dat2)):
	gc.append(dat2[i][5])

gc_percentile=np.percentile(gc, range(0,101))

gc_pp=[]
for i in range(0, len(dat2)):
	for j in range(0, 101):
		if dat2[i][5]>=gc_percentile[j] and dat2[i][5]<=gc_percentile[j+1]:
			gc_pp.append(j)
			break

# Normalize gc for each cell. nr: list of num. reads of each bin; gc_pp: corresponding gc_pp of the bin #
def gcNorm(nr, gc_pp=gc_pp):
	nr_a=[]
	for i in range(0,100):
		temp = []
		for j in range(0, len(nr)):
			if gc_pp[j]==i:
				temp.append(nr[j])
		nr_a.append(np.average(temp))
	# return(nr_a)
	nr_n=[]
	for i in range(0, len(nr)):
		nr_n.append(nr[i] * np.average(nr) / nr_a[gc_pp[i]])
	return(nr_n)

temp = []
for k in range(7,len(dat2[0])):
	cell_nr=[]
	for i in range(0, len(dat2)):
		cell_nr.append(dat2[i][k])
	temp.append(gcNorm(nr=cell_nr))

dat3 = []
for i in range(0, len(dat2)):
	temp1=dat2[i][0:7]
	for j in range(0, len(temp)):
		temp1.append(temp[j][i])
	dat3.append(temp1)

if args.report:
	output(dat3, args.odir+'/result.dat3_gc.txt')

### Cross batch normalize ###
# convert to raw CN estimate: dat4 #
cell_nr=[]
for i in range(7, len(dat3[0])):
	temp=[]
	for j in range(0, len(dat3)):
		if not('X' in dat3[j][0]) and not('Y' in dat3[j][0]):	# considering the sex chr. #
			temp.append(dat3[j][i])
	cell_nr.append(np.average(temp))

dat4=[]
for i in range(0, len(dat3)):
	temp=dat3[i][0:7]
	for j in range(0, len(cell_nr)):
		# temp.append(dat3[i][7 + j] / cell_nr[j] * 2)
		temp.append(dat3[i][7 + j] / cell_nr[j] * args.ploidy)
	dat4.append(temp)

if args.report:
	output(dat4, args.odir+'/result.dat4_cnvraw.txt')

# cross batch normalize: dat5 #
# estimat sex chr. ploidy#
x_ploidy=[]; y_ploidy=[]
for j in range(7, len(dat4[0])):
	tX=[]; tY=[]
	for i in range(0, len(dat4)):
		if 'X' in dat4[i][0]:
			tX.append(dat4[i][j])
		if 'Y' in dat4[i][0]:
			tY.append(dat4[i][j])
	x_ploidy.append(round(np.percentile(tX,50)))
	y_ploidy.append(round(np.percentile(tY,50)))

## the above is correct !! ##
dat5=[]
for i in range(0, len(dat4)):
	temp=dat4[i][0:7]
	for j in range(0, len(cell)):
		temp2=(np.sum(dat4[i][7:]) - dat4[i][7 + j]) / (len(cell)-1)
		if temp2 == 0:
			temp.append(0.0)
		elif not('X' in dat4[i][0]) and not('Y' in dat4[i][0]):	# considering the sex chr. #
			# temp.append(dat4[i][7 + j] / temp2 * 2);
			temp.append(dat4[i][7 + j] / temp2 * args.ploidy);
		elif 'X' in dat4[i][0]:
			temp3=0.0
			for k in range(0, len(cell)):
				if k!=j:
					temp3=temp3 + (3-x_ploidy[k]) * dat4[i][7 + k]
			temp3=temp3/(len(cell)-1)
			temp.append(dat4[i][7 + j]*(3-x_ploidy[j]) / temp3 * x_ploidy[j]);
		elif 'Y' in dat4[i][0]:
			temp3=0.0; temp4=0
			for k in range(0, len(cell)):
				if k!=j and y_ploidy[k]!=0:
					temp3=temp3 + dat4[i][7 + k]; temp4=temp4+1
			if temp4!=0:
				temp3=temp3/temp4
				temp.append(dat4[i][7 + j] / temp3 * y_ploidy[j])
			else:
				temp.append(0);
	dat5.append(temp)

## For each sample, normalize to the median cnv equals to 2 ##
nfactor=[]
for i in range(0, len(cell)):
	tmp=[]
	for j in range(0, len(dat5)):
		tmp.append(float(dat5[j][7+i]))
	nfactor.append(np.median(tmp))
# print(nfactor)

dat51=[]
for i in range(0, len(dat5)):
	temp=dat5[i][0:7]
	for j in range(0, len(cell)):
		# temp.append(2/nfactor[j] * dat5[i][7+j])
		# print([i,j,nfactor[j],dat5[i][7+j]])
		temp.append(args.ploidy/(nfactor[j]) * dat5[i][7+j])
	dat51.append(temp)

dat5=dat51
if args.report:
	output(dat5, args.odir+'/result.dat5_cnvcross.txt')

### Futher smoothing: dat6 ###
# step=floor(10/2)
step=floor(args.windowsize/2)
dat6=[]
for i in range(0, len(dat5)):
	tmp=dat5[i][0:7]
	for j in range(7, len(dat5[0])):
		p1=i-step; p2=i+step+1;
		if p1<0:	p1=0
		if p2>=len(dat5):	p2=len(dat5)-1
		temp=[]
		for k in range(int(p1),int(p2)):
			if dat5[k][0]==dat5[i][0]:
				temp.append(dat5[k][j])
		tmp.append(np.average(temp))
	dat6.append(tmp)

if args.report:
	output(dat6, args.odir+'/result.dat6_cnvsmooth.txt')

### CNV calling: dat7 ###
# Estimate sigma #
nr=[]
for j in range(7, len(dat6[0])):
	temp=[]
	for i in range(0, len(dat6)):
		temp.append(dat6[i][j])
	nr.append(temp)

sigma=[]
for j in range(0, len(dat6[0][7:])):
	sigma.append(abs(np.percentile(nr[j], 30.9)-2) + abs(np.percentile(nr[j], 69.1)-2))

# Define pdf #
def norm_pdf(x, mm, sg):
	return(exp(-(x-mm)**2/2/(sg**2))/sg/sqrt(2*pi))

# Define criteria #
cri=1 - 1.0 / args.genomesize * args.windowsize * args.binsize
# cri=1 - 1.0 / 3140000000 * 10 * 500000

dat7=[]
for i in range(0, len(dat6)):
	tmp=dat6[i][0:7]
	for j in range(0, len(dat6[0][7:])):
		pdf=[norm_pdf(x=dat6[i][7+j],mm=k,sg=sigma[j]) for k in range(0, 5)]
		lh = [k/np.sum(pdf) for k in pdf]
		a=[]
		for k in range(0,len(lh)):
			if lh[k] >= cri:
				a.append(k)
		if len(a)==0:
			a=['nan']
		tmp.append(a[0])
	dat7.append(tmp)

### output ###
output(dat7, args.odir+'/result.sccnv.txt')
