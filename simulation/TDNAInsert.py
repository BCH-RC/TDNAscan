#!/usr/bin/env python
import os.path
import sys
import getopt
import random
import re

lastLine = ""

def main(argv):
    sfile = "C:\\temp\\truncatedTDNA.fa"
    tfile = "C:\\temp\\referenceGonome.fa"
    ofile = ""

    try:
        opts, args = getopt.getopt(argv,"hs:t:o:",["sfile=","tfile=","ofile="])
    except getopt.GetoptError:
        print("TDNAscan -s <sourcefile> -t <targetfile> -o <outputfile>")
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print("TDNATruncate -s <sourcefile> -t <targetfile> -o <outputfile>")
            sys.exit()
        elif opt in ("-s", "--sfile"):
            sfile = arg
        elif opt in ("-t","--tfile"):
            tfile = arg
        elif opt in ("-o", "--ofile"):
            ofile = arg
   
    if ofile == "": 
        outputpath = os.path.dirname(sfile)
        ofile = outputpath + "\\outputTDNAInsertion.fa"
    parafile = os.path.dirname(ofile) + "\\" + os.path.basename(ofile).split('.')[0] + "_paras.txt"

    outputfile = open(ofile,"w+")
    parafile = open(parafile,"w+")
    parafile.write("ID\tChromosome\tInsert_Position\tTruncated_TDNA_Length\tTrucated_Start_Position\n")
    
    def formatFasta(inputStr, isLast):
        global lastLine
        size = 60
        inputStr = lastLine + inputStr
        parts = [inputStr[i:i + size] for i in range(0,len(inputStr),size)]
        lastLine = parts[len(parts) - 1]
        if isLast == False:
            parts = parts[:-1]
        return '\n'.join(parts) + '\n'

    def InsertTDNA(tdna, isLast):
        tdna = formatFasta(tdna, isLast)
        outputfile.writelines(tdna)

    if not os.path.isfile(sfile) or not os.path.isfile(tfile):
        print("File does not exist.")
    else:
        #with open(sourcefile,"r") as f:
        #    sourcecontent = f.readlines()

        #with open(targetfile,"r") as f:
        #    targetcontent = f.readlines()
        sourcefile = open(sfile,"r").read()
        targetfile = open(tfile,"r").read()
        dnacount = sourcefile.count(">")
        print(dnacount)

        randomstart = 1
        #firstLine = ""
        #if (targetfile.startswith(">")):
            #firstLine = targetfile.split('\n',1)[0]
            #outputfile.write(firstLine + "*****After Insertion*****" + "\n")
            #targetfile = re.sub('\n','',targetfile)
			#targetfile=re.sub(firstLine,"",targetfile)
			#targetfile = re.sub(firstLine,'',targetfile)
        tString = ""
        tLength = 0
        dnaID = ""
        t_dict = dict()
        with open(tfile,"r") as tf:
            t_content = tf.readlines()
        for tline in t_content:
            if(tline.startswith(">")):
                if(len(tString) > 0):						
                    t_dict[dnaID] = tString
                    tLength+=len(tString)
                    tString = ""
                dnaID = tline
                continue
            if(len(tline.strip()) > 0):
                tString = tString + tline.rstrip()
        if(len(tString) > 0):
             t_dict[dnaID] = tString
             tLength+=len(tString)			

        randomarray = random.sample(range(randomstart, tLength),dnacount)
        randomarray.sort()
        randomarray.insert(0, 0)

        with open(sfile,"r") as f:
            s_content = f.readlines()

        TDNAString = ""
        InsertedStringLength = 0
        index = 0
        dnaIDString = ""
        currentDNALength = 0
        trunctedTDNALength = ""
        reverseSymbol = ""
        for line in s_content:
            if (line.startswith(">")):                
                dnaIDString = line.strip()
                if(len(TDNAString) > 0):
                    currentTotalLength = 0
                    previousTotalLength = 0
                    previousString = ""
                    for k,v in t_dict.items():
                        previousTotalLength = currentTotalLength
                        currentTotalLength += len(v) 
                        previousString = v
                        if (currentTotalLength > randomarray[index]):
                            parafile.write(k.split(':')[0].replace(">","") + "\t")
                            currentDNALength = randomarray[index + 1] - previousTotalLength
                            parafile.write(str(currentDNALength) + "\t")
                            parafile.write(str(len(TDNAString)) + "\t")
                            parafile.write(trunctedTDNALength + "\t")
                            parafile.write(reverseSymbol + "\n")
                            if (index == 0):
                                outputfile.writelines(k)
                                targetString = v[randomarray[index] - previousTotalLength:randomarray[index + 1] - previousTotalLength]
                            else:
                                targetString = v[randomarray[index] + 1 - previousTotalLength:randomarray[index + 1] - previousTotalLength]
                            if(index > 0 and randomarray[index - 1] < previousTotalLength):
                                laveString = previousString[len(previousString) - previousTotalLength + randomarray[index - 1]:len(previousString)]
                                outputfile.writelines(laveString + "\n")
                                outputfile.writelines(k)
                            break
                        else:
                            currentDNALength = 0

                    #parafile.write(str(randomarray[index + 1]) + "\t")
                    #parafile.write(str(randomarray[index + 1] -
                    #randomarray[index] + InsertedStringLength) + "\t")
                    #randomarray[index] + InsertedStringLength) + "\t")
                    #parafile.write(str(len(TDNAString)) + "\n")
                    #parafile.write(TDNAString + "\n")

                    #if (index == 0):
                    #    targetString =
                    #    targetfile[randomarray[index]:randomarray[index + 1]]
                    #else:
                    #    targetString = targetfile[randomarray[index] +
                    #    1:randomarray[index + 1]]
                        
                    TDNAString = targetString + TDNAString
                    InsertTDNA(TDNAString, False)
                    
                    InsertedStringLength = InsertedStringLength + len(TDNAString)
                    
                    TDNAString = ""
                    print(index)
                    index += 1
                parafile.write(dnaIDString.split(',',1)[0].replace(">ID:","") + "\t")

                trunctedTDNALength = dnaIDString.split(',')[2].split(':')[1].strip()
                reverseSymbol = dnaIDString.split(',')[4]
                continue

            if(len(line.strip()) > 0):                
                TDNAString = TDNAString + line.rstrip()

        if(len(TDNAString) > 0):
            print(index)
            currentTotalLength = 0
            previousTotalLength = 0
            previousString = ""
            for k,v in t_dict.items():
                previousTotalLength = currentTotalLength
                currentTotalLength += len(v) 
                previousString = v
                if (currentTotalLength > randomarray[index]):
                    parafile.write(k.split(':')[0].replace(">","").split(" ")[0] + "\t")
                    currentDNALength = randomarray[index + 1] - previousTotalLength
                    parafile.write(str(currentDNALength) + "\t")
                    parafile.write(str(len(TDNAString)) + "\t")
                    parafile.write(trunctedTDNALength + "\t")
                    parafile.write(reverseSymbol + "\n")
                    targetString = v[randomarray[index] + 1 - previousTotalLength:randomarray[index + 1] - previousTotalLength]
                    TDNAString = targetString + TDNAString + v[randomarray[index + 1] - previousTotalLength:len(v)]
                    #targetString = previousString[len(previousString) -
                    #previousTotalLength + randomarray[index -
                    #1]:len(previousString)]
                    break
                else:
                    currentDNALength = 0

            #targetString = targetfile[randomarray[index]:randomarray[index +
            #1]]
            #TDNAString = targetString + TDNAString +
            #targetfile[randomarray[index + 1]:len(targetfile)]
            InsertTDNA(TDNAString, True)

            #parafile.write(str(randomarray[index + 1] - randomstart) + "\t")
            #parafile.write(str(randomarray[index + 1] - randomstart +
            #InsertedStringLength) + "\t")
            #parafile.write(str(len(TDNAString)) + "\n")
            #parafile.write(TDNAString + "\n")

    outputfile.close()
    parafile.close()

if __name__ == "__main__":
    main(sys.argv[1:])
