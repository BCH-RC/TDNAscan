#!/usr/bin/env python
import os.path
import sys
import getopt
import re

def main(argv):
    parafile = "C:\\temp\\at-tair10_mut500_paras.txt" 
    ofile = ""
    resultfiles = "C:\\temp\\10x_insertionBED.txt"

    try:
        opts, args = getopt.getopt(argv,"p:r:o:",["parafile=","resultfile=","ofile="])
    except getopt.GetoptError:
        print("CalAccuracy -p <parameterfile> -r <resultfile>")
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print("CalAccuracy -p <parameterfile> -r <resultfile>")
            sys.exit()
        elif opt in ("-p", "--parafile"):
            parafile = arg
        elif opt in ("-r","--resultfile"):
            resultfiles = arg
        elif opt in ("-o", "--ofile"):
            ofile = arg
            print(arg)
        else:
            assert False, "unhandled option"

 
    if ofile == "": 
        outputpath = os.path.dirname(parafile)
        ofile = outputpath + "\\output.txt"
    outputfile = open(ofile,"w+")

    with open(parafile,"r") as f:
        pcontent = f.readlines()

    for resultfile in resultfiles.split(" "):
        print(resultfile)
        outputfile.writelines(resultfile)
        outputfile.write("\n")

        with open(resultfile,"r") as f:
            rcontent = f.readlines()

        matchNumber = 0
        truePositives = []
        falseNegatives = []
        ID_Dictory = dict()
        MatchedChr = dict()
        for rline in rcontent:
            rinfo = rline.split("\t")
            chrID = rinfo[0]
            chrPos = rinfo[1]
            ismatch = False
            for pline in pcontent:
                if(pline.startswith("ID")):
                    continue
                else:
                    pinfo = pline.split("\t")
                    tdnaID = pinfo[0].strip()
                    pchrID = pinfo[1].split(" ")[0].strip(">")
                    pchrPos = pinfo[2].strip()
                    diff = abs(int(chrPos) - int(pchrPos))
                    if (diff < 50 and chrID == pchrID):
                        #print ( tdnaID + ":" + chrID +" "+ chrPos + " vs " +
                        #pchrPos + ", diff: "+ str(diff))
                        matchNumber += 1
                        if (ismatch == True):
                            print("closer matched one: " + pline + " vs " + MatchedChr[tdnaID] + "for the chr: " + rline)
                            outputfile.write("closer matched one: " + pline + " vs " + MatchedChr[tdnaID] + "for the chr: " + rline)
                            exChrPos = MatchedChr[tdnaID].split("\t")[2].strip()
                            extdnaID = MatchedChr[tdnaID].split("\t")[0].split(":")[1].strip()
                            exdiff = abs(int(chrPos) - int(exChrPos))
                            if(diff < exdiff):
                                #remove the previous from dict
                                MatchedChr.pop(extdnaID)
                            else:                                
                                continue
                        else:
                            ismatch = True

                        if(tdnaID in truePositives):
                            exChrPos = ID_Dictory[tdnaID].split("\t")[1]
                            exdiff = abs(int(exChrPos) - int(pchrPos))
                            if(diff < exdiff):
                                print("multiple matched TDNA: " + pline)
                                outputfile.write("mutiple matched TDNA: " + pline)
                                print("previous matched one: " + ID_Dictory[tdnaID])
                                outputfile.write("previous matched one: " + ID_Dictory[tdnaID])
                                ID_Dictory[tdnaID] = rline
                                print("replaced by better matched one: " + rline)
                                outputfile.write("replaced by better matched one: " + rline)
                            else:
                                print("find a closer one: " + str(tdnaID) + ", diff: " + str(diff) + " vs " + str(exdiff) + "\ncurrent chromosome: " + rline + "matched chromosome: " + ID_Dictory[tdnaID] +"\n" +"TDNA: " + pline)
                                outputfile.write("find a closer one: " + str(tdnaID) + ", diff: " + str(diff) + " vs " + str(exdiff) + "\ncurrent chromosome: " + rline)
                                outputfile.write("matched chromosome: " + ID_Dictory[tdnaID] +"\n")
                                outputfile.write("TDNA: " + pline)
                        else:
                            truePositives.append(tdnaID)
                            ID_Dictory[tdnaID] = rline
                            MatchedChr[tdnaID] = pline

            if (ismatch == False):
                #print ("doesn't match: " + rline)
                falseNegatives.append(rline)

        #print ('\n'.join(falseNegatives))
        outputfile.write("\nFalse Negatives: ")
        outputfile.write('\n'.join(falseNegatives))

        outputfile.write("\n\n missing (inserted but didn't find)\n\n")
        for pline in pcontent:
            if(pline.startswith("ID")):
                    continue
            else:
                pinfo = pline.split("\t")
                tdnaID = pinfo[0].strip()
                if(tdnaID not in truePositives):
                    outputfile.write(pline)
                    outputfile.write("\n")

        matchResult = "\ntotal: " + str(len(rcontent)) + ", match:" + str(matchNumber) + "\n"
        print(matchResult)
        outputfile.writelines(matchResult)

        precision = float(len(truePositives)) / float(len(pcontent))
        recall = float(len(truePositives)) / float(len(falseNegatives) + len(truePositives))
        fscore = 2 * precision * recall / (precision + recall)  
    
        outputresult = "fscore: " + str(fscore) + ", precision: " + str(precision) + ", recall: " + str(recall) + "\n"
        print(outputresult)
        outputfile.writelines(outputresult)
        outputfile.write("\n")
    outputfile.close()

if __name__ == "__main__":
    main(sys.argv[1:])