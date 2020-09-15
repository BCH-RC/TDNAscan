#!/usr/bin/env python
import os,sys,argparse,re,errno,time
import glob
import multiprocessing as mp

#def align2tdna(r1,r2,tdna_seq):
    #do something
start_time = time.time()    
IRbag = {}
PRTDNA = [] #paired end reads completely mapped on TDNA
TDNAname = ""
class IR:
    def _initi_(self):
        self.id = ""
        self.chr = ""
        self.pos = int
        self.soft = int
        self.match = int
        self.type = "" #whether soft-clipped, or discordant reads
        # self.genome_chr = ""
        # self.genome_pos = int
        # self.genome_soft = int
        # self.genome_match = int
        # self.genome_type = "" #whether soft-clipped, or discordant reads
    def calulateCigar(self,cigar,Rtype):
        pattern1 = "^(\d+)M(\d+)[SH]$"
        pattern2 = "^(\d+)[SH](\d+)M$"
        pattern3 = "^(\d+)M$"
        # m1_1 = re.search(pattern1,cigar)
        # m1_2 = re.search(pattern2,cigar)
        # m1_3 = re.search(pattern3,cigar)
        # m1_4 = re.search(pattern4,cigar)
        
        if "*" in cigar:
            self.soft = int(0)
            self.match = int(0)
            self.type = "*"
        elif re.search(pattern1,cigar):
            m1_1 = re.search(pattern1,cigar)
            self.soft = int(m1_1.group(2))
            self.match = int(m1_1.group(1))
            if Rtype == "R1":
                self.type = "1_1"
            else:
                self.type = "2_1"
        elif re.search(pattern2,cigar):
            m1_2 = re.search(pattern2,cigar)
            self.soft = int(m1_2.group(1))
            self.match = int(m1_2.group(2))
            if Rtype == "R1":
                self.type = "1_2"
            else:
                self.type = "2_2"
        elif re.search(pattern3,cigar):
            m1_3 = re.search(pattern3,cigar)
            self.soft = int(0)
            self.match = int(m1_3.group(1))
            if Rtype == "R1":
                self.type = "1_3"
            else:
                self.type = "2_3"
        else:
            self.soft = int(0)
            self.match = int(0)
            self.type = "*"

def captureTDNAname(tdna_seq):
    #?????????????????????????????????Check if fasta type or not??????????????????????????????????????????????????????
    FileIN = open(tdna_seq,"r")
    headline = FileIN.readline()
    if headline[0]==">" and len(headline.split(">")[1]) > 0:
        TDNAname = headline[1:].strip().split(" ")[0]
        return TDNAname
    else:
        sys.exit("TDNA is not in fasta format.")

def align2genome(R1,R2,genome,outfile,thread,directory):
    #align2genome("informativeRead_1.fq","informativeRead_2.fq","mt4_chr1_2Mb.fa")
    #bwa mem -T 20 -t 8 mt4_chr1_2Mb.fa informativeRead_1.fq informativeRead_2.fq  >informativeMt4.sam
    # samtools view -@ 8 -buS -q 30 informativeMt4.sam |samtools sort -@ 8 - -O bam -o informativeMt4_sort.bam
    #samtools index informativeMt4_sort.bam
    # bwa index genome data?????????????????????????????????
    cmd0 = "cp "+genome+" "+directory
    genome = os.path.basename(genome)
    genome = directory + "/" + genome
    cmd1 = "bwa index "+ genome
    cmd2 = "bwa mem -T 20 -t "+ str(thread) +" "+ genome + " " + R1 + " "+ R2 + " >"+outfile+".sam"
    #cmd2 = "bwa mem -t "+ str(thread) +" "+ genome + " " + R1 + " "+ R2 + " >"+outfile+".sam"    
    #cmd3 = "samtools view -@  "+ str(thread) + " -buS -q 30 " + outfile+".sam |samtools sort -@ "+ str(thread) + " - -O bam -o "+outfile+"_sort.bam"
    cmd3 = "samtools view -@  "+ str(thread) + " -buS " + outfile+".sam |samtools sort -@ "+ str(thread) + " - -O bam -o "+outfile+"_sort.bam"
    cmd4 = "samtools index " + outfile+"_sort.bam"
    os.system(cmd0)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)
    
def fileSplit(tdnaSAM,tmp_dir):
    cmd1 = "cp "+tdnaSAM+" "+tmp_dir
    #print cmd1
    tmp_sam = tmp_dir + "/" + os.path.basename(tdnaSAM)
    cmd2 = "python FileSplitter.py "+ tmp_sam +" 10000"
    #print cmd2
    
    cmd3 = "rm "+tmp_sam
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)

def filterSam(f):
    global TDNAname
    FileIN = open(f,"r")
    with open(f) as FileIN:
        bag = []
        readbag = []
        read1 = ""
        flag = False
        # directory = "./tmp"
        # try:
        #     os.makedirs(directory)
        # except OSError as e:
        #     if e.errno != errno.EEXIST:
        #         raise
        for line in FileIN:
            if "@" in line:
                continue
            data = line.strip().split("\t")
            if data[0] in readbag:
                if flag:
                    bag.append(read1)
                    bag.append(line.strip())
                    flag = False
                else:
                    if data[0].find(TDNAname)>-1:
                        bag.append(read1)
                        bag.append(line.strip())
                        read1 = ""
                        readbag = []
                        flag = False
                    else:
                        read1 = ""
                        readbag = []
                        flag = False
            else:
                readbag.append(data[0])
                read1 = line.strip() #save the first reads
                if line.strip().find(TDNAname)>-1:
                    flag = True
        # fname = directory+"/"+f
        # fo = open(fname,"w")
        # if fo:
        #     for b in bag:
        #         fo.write(b+"\n")
        # fo.close()
        return bag

def add2IRbag(rd1,data):
    global PRTDNA
    global IRbag
    #print rd1
    IRbag[rd1[0]] = IR()
    IRbag[rd1[0]].id = rd1[0]
    IRbag[rd1[0]].chr = rd1[2]
    IRbag[rd1[0]].pos = int(rd1[3]) #need to modify the position based on the clipped case
    #print "The data type for this R1"
    IRbag[rd1[0]].calulateCigar(rd1[5],"R1")
    
    #if read1 is informative, then continue, else, keep information of Read2
    if IRbag[rd1[0]].type == "1_3":
        # IRbag[data[0]] = IR()
        # IRbag[data[0]].id = data[0]
        # 
        # 
        # #print "The data type for this R1"
        # IRbag[data[0]].calulateCigar(data[5],"R2")
        # if IRbag[data[0]].type == "1_3": # Read1 and Read2 mapped to TDNA completely throw these reads out.
        if re.search("^(\d+)M$",data[5]):
            PRTDNA.append(data[0])
            IRbag.pop(data[0],None)
        else: #the other reads is "*"
            IRbag[data[0]] = IR()
            IRbag[data[0]].chr = data[2]
            IRbag[data[0]].pos = int(data[3]) #need to modify the position based on the clipped case
            IRbag[data[0]].calulateCigar(data[5],"R2")
            
    elif  IRbag[rd1[0]].type == "*":
        IRbag[data[0]].calulateCigar(data[5],"R2") 

def captureIR_tdna(sfile_tdna,IR1,IR2,TDNAname): #?????????????????????????parallel computing?????????????????????,tdna_name to replace t-dna_jiang
    FileIN = open(sfile_tdna,"r") # open fixed/filtered sam file
    
    #IR1 = "genome_informative_r1.fq"
    #IR2 = "genome_informative_r2.fq"
    FileIR1 = open(IR1,"w")
    FileIR2 = open(IR2,"w")
    #TEMP = IR1 + "_informativeReads.sam"
    #FileOUT = open(TEMP,"w")
    readbag = [] #store all read IDs
    read1 ="" #save the first reads
    mapq1 = int
    flag = False
    for line in FileIN:
        if "@" in line:
            continue
        data = line.strip().split("\t")
        # if int(data[4]) < 20:
        #     continue
        if data[0] in readbag:
            if flag:
                rd1 = read1.strip().split("\t")
                if int(data[4]) >= 20 or mapq1 >=20:
                    #FileOUT.write(line)
                    
                    # record the informative reads?????????????????????????????????????
                    add2IRbag(rd1,data) 
                    ##########write informative reads to files.
                    if data[0] in PRTDNA:
                        continue  #skip the paired end reads which are completely mapped to TDNA
                    FileIR1.write("@"+rd1[0]+"/1\n")
                    FileIR1.write(rd1[9] + "\n+\n")
                    FileIR1.write(rd1[10]+"\n")
                    read1 = ""
                    FileIR2.write("@"+data[0]+"/2\n")
                    FileIR2.write(data[9] + "\n+\n")
                    FileIR2.write(data[10]+"\n")
                    flag = False # renew the flag variable for next paired reads
                    del readbag[:]
                else:
                    continue
            else:
                if line.strip().find(TDNAname)>-1:
                    #FileOUT.write(line)
                    rd1 = read1.strip().split("\t")
                    if int(data[4]) >= 20 or mapq1 >=20:
                        add2IRbag(rd1,data) 
                        ##########write informative reads to files.
                        if data[0] in PRTDNA:
                            continue #skip the paired end reads which are completely mapped to TDNA
                        FileIR1.write("@"+rd1[0]+"/1\n")
                        FileIR1.write(rd1[9] + "\n+\n")
                        FileIR1.write(rd1[10]+"\n")
                        read1 = ""
                        FileIR2.write("@"+data[0]+"/2\n")
                        FileIR2.write(data[9] + "\n+\n")
                        FileIR2.write(data[10]+"\n")
                        del readbag[:]
                    else:
                        continue
                else:
                    reads1 = "" #renew the first read variable for next paired reads
                    mapq1 = int
                    del readbag[:]

        else:
            readbag.append(data[0])
            read1 = line #save the forward reads
            mapq1 = int(data[4])
            if line.strip().find(TDNAname)>-1:
                flag = True
                
    FileIN.close()
    FileIR1.close()
    FileIR2.close()
    #FileOUT.close()
    
# def split_list(a_list):
#     a_list.sort()
#     half = len(a_list)/2
#     return a_list[:half], a_list[half:]

def analyzeClust(clustbag):
    #(chr,start,readType,tdna_start,orientation)
    insertionTDNA = ()
    startPos_bag = []
    startPos_bag_estimate = [] #used to estimate the insertion position but not the exact one if only 'DIR' reads identified
    tdna_pos_st_bag = []
    tdna_pos_end_bag = []
    pos_mode = int
    tdna_pos_st_mode = ""
    tdna_pos_end_mode = ""
    orientation_mode = ""
    
    chromosome = ""
    orientation_bag = []
    clr_n = 0
    dir_n = 0
    #print "-----------------------------"
    for (chr,start,readType,tdna_direction,tdna_pos,orientation) in clustbag:
        chromosome = chr
        #print (chr,start,readType,tdna_pos,orientation)
        #remove DIR type for position prediction

        if readType == "CLR":
            clr_n = clr_n + 1
            startPos_bag.append(start)
            orientation_bag.append(orientation)
            if tdna_direction == "TSP":
                tdna_pos_st_bag.append(tdna_pos)
            elif tdna_direction == "TEP":
                tdna_pos_end_bag.append(tdna_pos)
        elif readType == "DIR":
            dir_n = dir_n + 1
            startPos_bag_estimate.append(start)
    
    if not startPos_bag and startPos_bag_estimate: #only 'DIR' reads identified
        pos_mode = max(set(startPos_bag_estimate),key=startPos_bag_estimate.count)
        tdna_pos_st_mode = '-'
        tdna_pos_end_mode = '-'
        orientation_mode = '*'
    elif startPos_bag:
        pos_mode = max(set(startPos_bag),key=startPos_bag.count)
        orientation_mode = max(set(orientation_bag),key=orientation_bag.count)
        if tdna_pos_st_bag:
            tdna_pos_st_mode = str(max(set(tdna_pos_st_bag),key=tdna_pos_st_bag.count))
        else:
            tdna_pos_st_mode = "-"
        if tdna_pos_end_bag:
            tdna_pos_end_mode = str(max(set(tdna_pos_end_bag),key=tdna_pos_end_bag.count))
        else:
            tdna_pos_end_mode = "-"
            
        
    # elif len(tdna_pos_bag) > 1:
    #     pos_mode = max(set(startPos_bag),key=startPos_bag.count)
    #     #print tdna_pos_bag
    #     tdna_pos_st, tdna_pos_end = split_list(tdna_pos_bag)
    #     tdna_pos_st_mode = max(set(tdna_pos_st),key=tdna_pos_st.count)
    #     tdna_pos_end_mode = max(set(tdna_pos_end),key=tdna_pos_end.count)
    #     orientation_mode = max(set(orientation_bag),key=orientation_bag.count)
    # else:
    #     pos_mode = max(set(startPos_bag),key=startPos_bag.count)
    #     tdna_pos_st_mode = tdna_pos_bag[0]
    #     tdna_pos_end_mode = tdna_pos_bag[0]
    #     orientation_mode = max(set(orientation_bag),key=orientation_bag.count)

    
    suppRead = "CLR:"+str(clr_n)+",DIR:"+str(dir_n)
    tdna_info = "tdna_st:"+tdna_pos_st_mode+",tdna_end:"+tdna_pos_end_mode
    insertionTDNA = (chromosome,pos_mode,suppRead,tdna_info,orientation_mode)
    return insertionTDNA

def clusterIR(informativeGenome,insertionRead,minRD,winCLR,winDIR):
    informativeGenomeSAM = informativeGenome + ".sam"
    # FileIN = open(informativeGenomeSAM,"r") # open sam file
    # FileOUT = open(insertionRead,"w")
    # FileOUT1 = open(insertionBED,"w")
    insertion = [] # [(chr,start,readType(CLR,DIR),tdna_start,tdna_end,orientation)],[(chr,start,supportReads,tdna_start,tdna_end,orientation)]
    insertionbag = []   #????????????????????????????????????????????????????????????????????????
    # INSbag = []
    # read1 = ""
    # Rtype = ""
    tmp_reads = open(informativeGenome+"_IR.txt","w") #??????????????????????????????????????????????????????????????????????????????????
    with open(informativeGenomeSAM) as FileIN,open(insertionRead,"w") as FileOUT:
        for line in FileIN:
            if "@" in line:
                continue
            data = line.strip().split("\t")
            if not IRbag.has_key(data[0]):
                print "Warning: this read is not in the IR bag " + data[0]
            #here we only care about the read (one pair of the paired end reads) whether it contains information in pattern 1 or 2 or 3 below.
            if int(data[4]) < 20:   
                continue
            #if data[0] in INSbag:
            # 6 Cases
            # case a: ref 45M105S tdna 45S105M
            # case b: ref 45S105M tdna 45M105S
            # case c: ref 150M     tdna 150M
            # case d: ref 45M105S tdna 105M45S
            # case e: ref 45S105M tdna 105S45M
            # case f: ref 150M     tdna 150M
            pattern1 = "^(\d+)M(\d+)[SH]$"
            pattern2 = "^(\d+)[SH](\d+)M$" #1_2
            pattern3 = "^(\d+)M$"
            m1_1 = re.search(pattern1,data[5])
            m1_2 = re.search(pattern2,data[5])
            m1_3 = re.search(pattern3,data[5])
            #if IRbag[data[0]].type == "1_1" or IRbag[data[0]].type == "1_2": #reverse read
            if m1_1:
                soft = int(m1_1.group(2))
                match = int(m1_1.group(1))
                if IRbag[data[0]].type == "1_2" or IRbag[data[0]].type == "2_2": ##case a
                    if abs(IRbag[data[0]].match -soft)<= 5: # the difference of mapped reads to plant genome and soft clipped reads on TDNA
                        pos = int(data[3]) + match
                        #insertion.append((data[2],pos,"CLR",IRbag[data[0]].pos,"+"))
                        insertion.append((data[2],pos,"CLR","TSP",IRbag[data[0]].pos,"+")) #TSP:T-DNA Start Postion
                        tmp_reads.write(line)
                if IRbag[data[0]].type == "1_1" or IRbag[data[0]].type == "2_1": ##case d
                    if abs(IRbag[data[0]].match -soft)<= 5: # the difference of mapped reads to plant genome and soft clipped reads on TDNA
                        pos = int(data[3]) + match
                        insertion.append((data[2],pos,"CLR","TEP",(IRbag[data[0]].pos + IRbag[data[0]].match),"-")) #TEP:T-DNA End Postion
                        tmp_reads.write(line)
            elif m1_2:
                soft = int(m1_2.group(1))
                match = int(m1_2.group(2))
                if IRbag[data[0]].type == "1_1" or IRbag[data[0]].type == "2_1": ##case b
                    if abs(IRbag[data[0]].match - soft)<= 5:
                        pos = int(data[3])
                        insertion.append((data[2],pos,"CLR","TEP",(IRbag[data[0]].pos + IRbag[data[0]].match),"+"))
                        tmp_reads.write(line)
                if IRbag[data[0]].type == "1_2" or IRbag[data[0]].type == "2_2": ##case e
                    if abs(IRbag[data[0]].match - soft)<= 5:
                        pos = int(data[3])
                        insertion.append((data[2],pos,"CLR","TSP",IRbag[data[0]].pos,"-"))
                        tmp_reads.write(line)
            elif m1_3:
                #discordant reads, we need to guess the insertion location
                #print data[0] + "-----------------------DIR"
                if IRbag[data[0]].type == "1_3" or IRbag[data[0]].type == "2_3":
                    pos = int(data[3]) + 250 #??????????????????????????????????assume the length of sequence fragment is 500bp
                    insertion.append((data[2],pos,"DIR","-","-","*"))
                    tmp_reads.write(line)
            else: #read the reverse read
                #print "-----------------"+str(data)
                continue
    
        win = int
        winflag = False #the previous read is not DIR type
        clustbag = []
        bag1 = ()
        #insertionTDNA = ()
        for (chr,start,readType,tdna_direction,tdna_pos,orientation) in sorted(insertion,key=lambda element:(element[0],element[1])):
            #write the insertion to
            #print (chr,start,readType,tdna_pos,orientation)
            FileOUT.write(str((chr,start,readType,tdna_direction,tdna_pos,orientation))+"\n")
            if winflag:
                #win = 500
                win = winDIR
            elif readType == "CLR":
                #win = 5
                win = winCLR
            elif readType == "DIR":
                #win = 500
                win = winDIR
            #print "-------------windown size used:"+str(win)
            if not clustbag:
                bag1 = (chr,start,readType,tdna_direction,tdna_pos,orientation)
                clustbag.append((chr,start,readType,tdna_direction,tdna_pos,orientation))
            else:
                #print (chr,start,readType,tdna_direction,tdna_pos,orientation)
                if chr == bag1[0] and abs(start - bag1[1])<= win:
                    clustbag.append((chr,start,readType,tdna_direction,tdna_pos,orientation))
                    if readType == "DIR":
                        winflag = True #make the windown size 500 
                        
                    bag1 = (chr,start,readType,tdna_direction,tdna_pos,orientation)
                else: # we should clust the bag and write the result to insertion file
                    insertionTDNA = analyzeClust(clustbag)
                    #write the insertion [(chr,start,supportReads,tdna_start,tdna_end,orientation)]
                    clr_n = int(insertionTDNA[2].split(",")[0].split(":")[1])
                    dir_n = int(insertionTDNA[2].split(",")[1].split(":")[1])
                    t_start = insertionTDNA[3].split(",")[0].split(":")[1]
                    t_end = insertionTDNA[3].split(",")[1].split(":")[1]
                    #if t_start != "-" and t_end != "-":
                    if (clr_n + dir_n) >=minRD:
                        insertionbag.append(insertionTDNA)
                    # if t_start == "-" or t_end == "-":
                    #     if (clr_n + dir_n) >=3:
                    #         insertionbag.append(insertionTDNA)
                    # else:
                    #     t_start = int(t_start)
                    #     t_end = int(t_end)
                    #     if (clr_n + dir_n) >=3 and (t_end-t_start)>=25:
                    #         insertionbag.append(insertionTDNA)
                        #FileOUT1.write(insertionTDNA[0]+"\t"+str(insertionTDNA[1])+"\t"+"\t".join(insertionTDNA[2:5])+"\n") #(chr,pos_mode,suppRead,tdna_info,orientation_mode)
                    clustbag = []
                    bag1 =  (chr,start,readType,tdna_direction,tdna_pos,orientation)
                    clustbag.append((chr,start,readType,tdna_direction,tdna_pos,orientation))
    
            
        # the end of file
        insertionTDNA = analyzeClust(clustbag)
        clr_n = int(insertionTDNA[2].split(",")[0].split(":")[1])
        dir_n = int(insertionTDNA[2].split(",")[1].split(":")[1])
        t_start = insertionTDNA[3].split(",")[0].split(":")[1]
        t_end = insertionTDNA[3].split(",")[1].split(":")[1]

        #if t_start != "-" and t_end != "-":
        if (clr_n + dir_n) >=minRD : #and (t_end-t_start)>=25
            insertionbag.append(insertionTDNA)
        # #     FileOUT1.write(insertionTDNA[0]+"\t"+str(insertionTDNA[1])+"\t"+"\t".join(insertionTDNA[2:5])+"\n")
    tmp_reads.close()
    #os.system("sort -k 3,3 -k 4,4n "+informativeGenome+"_IR.txt >"+informativeGenome+"_IR_sort.txt")
    return insertionbag        
        # #extract sequences of insertion sites and put it into one file and align all reads to it and get the average read depth
        # for (chr,pos_mode,suppRead,tdna_info,orientation_mode) in insertionTDNA:
            

def detectZygosity(insertionbag,fq1,fq2,reference,thread,directory,project):
    insertionBED = directory+ "/5."+project+"_insertion.bed"
    insertionSeq = directory+ "/5.insertionSeq.fa"
    outfile = directory+ "/5.insertion"
    seq_win = 500
    with open(insertionBED,"w") as file_bed:
        file_bed.write("Chr\tBreakpoint\tSuppRead\tTDNA_info\tOrientation\tFreq\n")
        seq_region = ""
        for (chr,pos_mode,suppRead,tdna_info,orientation_mode) in insertionbag:
            seq_region = seq_region + ' "'+chr+":" + str((pos_mode-seq_win))+"-" + str((pos_mode+seq_win))+ '"'
        cmd1 = 'samtools faidx '+reference+' '+seq_region+' >'+insertionSeq
        #print cmd1
        os.system(cmd1)
        cmd2 = "bwa index "+insertionSeq
        cmd3 = "bwa mem -T 20 -t "+ str(thread) +" "+ insertionSeq + " " + fq1 + " "+ fq2 + " >"+outfile+".sam"
        #cmd4 = "samtools view -@  "+ str(thread) + " -buS -q 30 " + outfile+".sam |samtools sort -@ "+ str(thread) + " - -O bam -o "+outfile+"_sort.bam"
        cmd4 = "samtools view -@  "+ str(thread) + " -buS " + outfile+".sam |samtools sort -@ "+ str(thread) + " - -O bam -o "+outfile+"_sort.bam"
        cmd5 = "samtools index "+outfile+"_sort.bam"
        os.system(cmd2)
        os.system(cmd3)
        os.system(cmd4)
        os.system(cmd5)
        for (chr,pos_mode,suppRead,tdna_info,orientation_mode) in insertionbag:
            
            clr_n = int(suppRead.split(",")[0].split(":")[1])
            #5bp around the insertion site
            insertion_region = ' "'+chr+":" + str((pos_mode-seq_win))+"-" + str((pos_mode+seq_win))+ '"'
            region = str((seq_win-5))+"-" + str((seq_win+5))
            tmp_file = directory+ "/tmp.txt"
            #cmd6 = "samtools depth -r "+insertion_region+":"+ region + " " + outfile+"_sort.bam >"+tmp_file
            cmd6 = "samtools view " + outfile+"_sort.bam " +insertion_region+":"+ region + " >"+tmp_file
            os.system(cmd6)
            spanread = 0
            freq = float
            with open(tmp_file) as tf:
                for line in tf:
                    data = line.strip().split("\t")
                    pattern = "^(\d+)M"
                    m = re.search(pattern,data[5])
                    if m and int(data[4])>=30:
                        match = int(m.group(1))
                        pos = int(data[3])+match
                        if int(data[3]) <= (seq_win-5) and pos >= (seq_win+5):
                            spanread = spanread + 1
                            #print line
            if clr_n+spanread == 0:
                freq = float(0)
            else:
                freq = clr_n/float(clr_n+spanread)
            #print freq
            if clr_n == 0:
                file_bed.write(chr+"\t~"+ str(pos_mode)+"\t"+suppRead+"\t"+tdna_info+"\t"+orientation_mode+"\t"+str(freq)+"\n")
            else:
                file_bed.write(chr+"\t"+ str(pos_mode)+"\t"+suppRead+"\t"+tdna_info+"\t"+orientation_mode+"\t"+str(freq)+"\n")
            cmd7 = "rm "+tmp_file
        #os.system(cmd7) #??????????????????????????????????????????
                
                    
        


#define parameter
#---------------------------------------------------------------------------------############
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="detect informative reads")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-1', action='store', dest='fq1', help="read 1 of paired end reads",required=True) ##????????????????????????only 1  parameter??????????????????????
    parser.add_argument('-2', action='store', dest='fq2', help="read 2 of paired end reads",required=True)
    parser.add_argument('-g', action='store', dest='reference', help="plant genome file",required=True)
    parser.add_argument('-t', action='store', dest='tdna_seq', help="t-DNA sequence file",required=True)
    parser.add_argument('-p', action='store', dest='project', help="the project name",required=True)
    parser.add_argument('-n', action='store', default='3', type=int, dest='minRD', help="the minimal total number of informative reads (clipped and discordant reads [default:3]")
    parser.add_argument('-a', action='store', default='5', type=int, dest='winCLR', help="the window size of clustering soft clipped reads [default:3]")
    parser.add_argument('-b', action='store', default='500', type=int, dest='winDIR', help="the length of DNA fragment in NGS data [default:500]")
    parser.add_argument('-@', action='store', default='8', type=int, dest='thread', help="default thread: 8")

    
    args = parser.parse_args()    
    fq1 = args.fq1
    fq2 = args.fq2
    reference = args.reference
    tdna_seq = args.tdna_seq
    #tfile = args.tfile
    project = args.project
    minRD = args.minRD
    winCLR = args.winCLR
    winDIR = args.winDIR
    thread = args.thread
    insertionbag = []
    lines = []

    #python tdnascan.py -1 mt4_chr1_20x_mut_tdna_1.fq -2 mt4_chr1_20x_mut_tdna_2.fq -t t-dna_elison.fa -g mt4_chr1_2Mb.fa -p tdna
    tnd_len = 10000 #10kb ??????????????????????????????
    ##########################Step 1: map all reads to T-DNA sequence to generate a SAM file
    print "Running the step1: map all reads to T-DNA sequence to generate a SAM file "
    #def align2genome(R1,R2,genome,outfile):
    #align2genome("mt4_chr1_20x_mut_tdna_1.fq","mt4_chr1_20x_mut_tdna_2.fq","t-dna_elison.fa","informativeTDNA")
    directory = "./" + project
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
        
    informativeTDNA = directory+ "/1.TDNA"
    #change the first line of T-DNA sequence
    with open(tdna_seq) as tdna_seq_fh:
        lines = tdna_seq_fh.readlines()
    lines[0] = ">"+project+"\n"
    with open (tdna_seq, "w") as tdna_seq_fh:
        tdna_seq_fh.writelines(lines)
    # align all reads to T-DNA sequence
    align2genome(fq1,fq2,tdna_seq,informativeTDNA,thread,directory)
    ##########################Step 2: capture all informative reads from SAM file from Step1
    print "Running the step2: capture all informative reads from SAM file from Step1 "
    tdnaSAM = directory+ "/1.TDNA.sam"
    TDNAname = captureTDNAname(tdna_seq)
    #captureIR_tdna("tdna_sample.sam")
    tmp_dir = "./" + project+"/tmp"
    try:
        os.makedirs(tmp_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
        
    files = glob.glob(tmp_dir +"/*")
    for f in files:
        os.remove(f)
    fileSplit(tdnaSAM,tmp_dir)
    #filter out splitted sam file via parallel computing
    
    print  "parallel Processing"
    start_time1 = time.time()
    
    #thread = 8
    
    #init objects
    pool = mp.Pool(thread)  ##big memory problem from https://stackoverflow.com/questions/14677287/multiprocessing-with-large-data
    
    
    filenames = glob.glob(tmp_dir +"/*.sam")
    #print filenames
    tdnaSAM_fix = directory+ "/2.informativeTDNA_fix.sam"
    with open(tdnaSAM_fix,"w") as f:
        #print pool.imap(filterSam,filenames)
        for result in pool.imap(filterSam,filenames):
            if result:
                for r in result:
                    f.write(r+"\n")
    
        
    pool.close()
    pool.join()
    
    
    end_time1 = time.time()        
     
    print("Time for multiple processes: %ssecs" % (end_time1 - start_time1))
    
    #capture all informative reads
    informative_genome_r1 = directory+ "/2.informative_genome_r1.fq"
    informative_genome_r2 = directory+ "/2.informative_genome_r2.fq"
    
    captureIR_tdna(tdnaSAM_fix,informative_genome_r1,informative_genome_r2,TDNAname)

    # ##########################Step 3: map all informative reads to plant genome
    print "Running the step3: map all informative reads to plant genome "
    informativeGenome = directory+ "/3."+project+"_informativeGenome"
    align2genome(informative_genome_r1, informative_genome_r2, reference,informativeGenome,thread,directory)
    ##########################Step 4: output all candidate insertion location and report the detail information of truncated T-DNA including how many base pairs truncated at both side of TDNA and the insertion direction
    print "Running the step4:extract informative reads and cluster all informative reads based on location "
    #clusterIR("informativeMt4.sam")

    insertionRead = directory+ "/4.insertionRead.txt"
    
    insertionbag = clusterIR(informativeGenome,insertionRead,minRD,winCLR,winDIR)
    
    print "Running the step5: extract sequences blast all reads to them, and identify zygosity "
    #extract reads using samtools
    
    detectZygosity(insertionbag,fq1,fq2,reference,thread,directory,project)
    
    #blast reads to above sequences
    
    #extract the average reads near the 5bp up/downstream of insertion sites and calculate the zygosity = clipped reads/average reads 
    
    end_time = time.time()        
 
    print("Time for TDNA detection: %ssecs" % (end_time - start_time))


# def main():
#     try:
#         TDNAscan();
#     except IOError as err:
#         sys.stderr.write("IOError " + str(err) + '\n');
#         return;
#     except StandardError as err:
#         sys.stderr.write("StandardError\n");
#         handleException(err);
#         sys.exit(1);
# 
# if __name__ == "__main__":
#     sys.exit(main())
#     (END)
