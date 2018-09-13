#!/usr/bin/env python


#encoding:utf-8
"""
tdnaAnnot.py
created by Liang Sun on 2018-08-13
Copyright (c) 2018 Liang Sun. All rights reserved.
Updated and maintained by Liang Sun since Aug 2018

tdnaAnnot annotate identified T-DNA insertions using gff3 annotation file.

usage:
	tdnaAnnot [options] -t type -i <input file> -o <output file>
	
example:
	python tdnaAnnot.py -i result/AF_unique.txt -t s -f Mtruncatula_285_Mt4.0v1.gene.gff3 -o result/AF_unique_annot.txt

"""

import sys
import os
import argparse
import time


        
#read file
def processGFF(file):
    geneInfo = {} # chr = [start, end, geneID]
    FileIN = open(file,"rU")
    for line in FileIN:
        if "##" in line:
            continue
        else:
            data = line.strip().split("\t")
            if data[2] == "gene":
                geneID = data[8].strip().split(";")[1].split("=")[1]
                if geneInfo.has_key(data[0]):
                    geneInfo[data[0]].append((int(data[3]),int(data[4]),geneID))
                else:
                    geneInfo[data[0]] = [(int(data[3]),int(data[4]),geneID)]
    
    FileIN.close()
    return geneInfo
        
def processVar(file):
    varInfo = {}
    FileIN = open(file,"rU")
    FileIN.readline()
    for line in FileIN:
        data = line.strip().split("\t")
        data2 = int(data[1]) + int(1)
        varInfo[(data[0],int(data[1]),data2)] = line.strip()
    
    FileIN.close()
    return varInfo

def annotateVar(varInfo,geneInfo):  
    geneAnnot = {}
    for chr,start,end in varInfo:
        #print chr,start,end
        flag = 0
        if not geneInfo.has_key(chr):
            geneAnnot[(chr,start,end)] =[]
            continue
        for s,e,id  in geneInfo[chr]:
            if (start>= s and start <= e) or (end>=s and end <= e):
                if geneAnnot.has_key((chr,start,end)):
                    geneAnnot[(chr,start,end)].append(id)
                else:
                    flag = 1
                    geneAnnot[(chr,start,end)] = [id]
        if flag == 0:
            geneAnnot[(chr,start,end)] =[]
    
    return geneAnnot
                
            
            
        
    
############################# declare variables #############################
t0 = time.time()

ifile = ""
ffile = ""
ofile = ""
geneInfo = {} # chr = geneID ,start, end
varInfo = {} # (chr, pos) = info
geneAnnot = {} # (chr, pos) = geneID

############################## process command line arguments #############################
parser = argparse.ArgumentParser(description="tdnaAnnot annotate identified T-DNA insertions using gff3 annotation file")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-i', action='store', dest='ifile', help="T-DNA BED file")
parser.add_argument('-f', action='store', dest='ffile', help="gff3 annotation file")
parser.add_argument('-o', action='store', dest='ofile', help="the output file name")
#parser.add_argument('-r', action='store', default='8', type=float, dest='cutoff', help="reads covrage cutoff")
args = parser.parse_args()

ifile = args.ifile
ffile = args.ffile
ofile = args.ofile
#cutoff = args.cutoff (how many )


if len(sys.argv) == 1:
	print parser.print_help()
	sys.exit("error: give me your input and output files!")

############################# read AF files #############################
else :
    geneInfo = processGFF(ffile)
    varInfo = processVar(ifile)
    geneAnnot = annotateVar(varInfo,geneInfo)
    #print geneAnnot

############################# write unique variance to output file #############################
FileOUT = open(ofile,"w")
    #FileOUT.write("Del#\tchr\tstart_position\tend_position\tdeletionLength\tbreakpoint_pos\tsupportRead\tdel_mutant\tdel_control\tHomo_Unique\tGenes\n")
for chr, start,end in sorted(varInfo):
    FileOUT.write(varInfo[(chr,start,end)]+"\t"+str(geneAnnot[(chr,start,end)])+"\n")
    
FileOUT.close()

t1 = time.time()
totalTime = t1-t0
print "Run time: " + str(totalTime) + "sec"