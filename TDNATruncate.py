#!/usr/bin/env python
import os.path
import sys
import getopt
import re
from random import randint
from Bio.Seq import Seq

def main(argv):
    ifile = "C:\\temp\\referenceGonome.fa" 
    ofile = ""
    minLength = 50
    number = 500

    try:
        opts, args = getopt.getopt(argv,"hi:o:n:",["ifile=","ofile=","number="])
    except getopt.GetoptError:
        print("TDNAscan -i <inputfile> -o <outputfile>")
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print("TDNATruncate -i <inputfile> -o <outputfile>")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-o","--ofile"):
            ofile = arg
        elif opt in ("-n", "--number"):
            number = int(arg)

    if ofile == "":
        outputpath = os.path.dirname(ifile)
        ofile = outputpath + "\\truncatedTDNA.fa"
    outputfile = open(ofile,"w+")
    
    def formatFasta(inputStr):
        size = 60
        parts = [inputStr[i:i + size] for i in range(0,len(inputStr),size)]
        return '\n'.join(parts)

    def TruncateTDNA(tdna, name):
        newtdna = tdna.replace("\n","")
        strLength = len(newtdna)
        for i in range(number):
            truncatedLength = randint(minLength, strLength)
            if (truncatedLength % 2 == 0):
                truncatedLength = strLength
            randomLength = strLength - truncatedLength
            randomPos = randint(0, randomLength)
            truncatedTDNA = newtdna[randomPos:truncatedLength + randomPos]
            print("TDNA length: " + str(strLength) + ", start position: " + str(randomPos) + ", truncatedTDNA length: " + str(truncatedLength))
            outputfile.write(">ID: " + str(i))
            outputfile.write(", TDNA length: " + str(strLength) + ", start position: " + str(randomPos) + ", truncatedTDNA length: " + str(truncatedLength))

            if(randint(0,1) == 1):
                seq = Seq(truncatedTDNA)
                truncatedTDNA = str(seq.reverse_complement())
                outputfile.write(",-\n")
            else:
                outputfile.write(",+\n")
            
            formatedtruncatedTDNA = formatFasta(truncatedTDNA)

            outputfile.write(formatedtruncatedTDNA)
            outputfile.write("\n\n")    
        #outputfile.write(">original TDNA\n")
        #outputfile.write(name)
        #outputfile.write(tdna)
        #outputfile.write("\n")


    if not os.path.isfile(ifile):
        print("File does not exist.")
    else:
        with open(ifile,"r") as f:
            content = f.readlines()   

        
    TDNAString = ""
    TDNAName = ""
    for line in content:
        if(line.startswith(">")):
            TDNAName = line
            if(len(TDNAString) > 0):
                TruncateTDNA(TDNAString, TDNAName)

            TDNAString = ""
            continue 
        if (len(line) > 0):
            TDNAString = TDNAString + line.rstrip()
   
    if(len(TDNAString) > 0):
       TruncateTDNA(TDNAString, TDNAName)



    outputfile.close()

if __name__ == "__main__":
    main(sys.argv[1:])
