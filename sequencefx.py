#Code by Zach Frye for the experimentally customized functions
#Processes the sample sequences and compiles the mutations

import itertools
import operator
import time
import os
import sys
import datetime

#Sets up the parameters for the set of alignment keys 
KeySize = 20 #Determines the hash size of the WT gene
KeyOverlap = 0 #Overlap between alignment hashes
KeyCoverage = 4 #Required number of keys aligned to sequence
WTkey = sequencefx.KeyMaker(WTseq) #Keys for alignment

FailC = 0   #Counts the number of Failed sequences
MutLimit = 10   #Controls the number of mutations tolerated in a given sequence
NT = '0ACGTN'   #Establishes the columns for the data table
seqOverlap = 20
#Set up the Qscore scale
QscoreKey = {}
QscoreKey = proteinfx.QsKey()

#Makes a set of sequence from the WT gene sequence with specified size and overlap for indexing to sample sequence.
def KeyMaker(seq):
    seqLength = len(seq)
    WTkey = {}

    print('WTlength: %s, Key with overlap: %s.' % 
          (seqLength, int(seqLength / (sequencefx.KeySize - sequencefx.KeyOverlap))+1))

    for i in range(0, seqLength, (sequencefx.KeySize - sequencefx.KeyOverlap)):
        WTkey[i] = seq[i:(i + sequencefx.KeySize)-1]

    #Additional considerations for the first and last WTseq keys
    WTkey[1] = seq[1:sequencefx.KeySize]
    WTkey[seqLength-sequencefx.KeySize] = seq[-(sequencefx.KeySize):]

    sorted_WTkey = sorted(WTkey.items(), key=operator.itemgetter(0))
    WTkey = {}
    for i in sorted_WTkey[:-1]:
        WTkey[i[0]] = i[1]
    return WTkey

#Used as a line by line analysis to save memory
def QuickProcess(seq, WTMat, MutCor, FailSeq, AATable):    
    if len(seq) > 0:
        skip3, AATable, MutCor = sequencefx.KeyAlign(seq, WTMat, MutCor, AATable)
        if  skip3 == 0:
            skip3, AATable, MutCor = sequencefx.KeyAlign(proteinfx.RevComp(seq), WTMat, MutCor, AATable)
            if skip3 == 0:
                #FailC += 1
                FailSeq.append(seq) #In case the failed sequences are wanted
                print(seq)
    return skip3, FailSeq, AATable, MutCor

    #Matches the overlaping sequences (under construction)
def Overlap(Fhit, Rhit, seqF, seqR):
    Dual = ''
    if Fhit < Rhit:
        OL = len(seqF) + Rhit - Fhit
        seq1 = seqF[(Rhit - Fhit - 1):] 
        seq2 = seqR[:(len(seqF) + Rhit - Fhit)]
    else:
        TempSeq = seqF
        seqF = seqR
        seqR = TempSeq

        OL = len(seqR) + Rhit - Fhit
        seq1 = seqF[(Rhit - Fhit):] 
        seq2 = seqR[:(len(seqR) + Rhit - Fhit)]
    return Dual

#Used to find an InDel (under construction)
def FindIndel(seq, StartNT):
    Mcount = 0
    InDelIndex = 0
    MutList = []
    if StartNT == 15:
        Aligned = WTseq[(StartNT):(StartNT + len(seq))]
    else:
        Aligned = WTseq[StartNT:(StartNT + len(seq))]

    if len(Aligned) < len(seq):
        seq = seq[0:len(Aligned)]

    for i in range(0,len(seq)):
        if seq[i] != Aligned[i]:
            MutList.append(i)
    return

#Searches the sample sequence with the Key dictionary
def KeySearch(seq):
    place, StartNT, indel, skip, diff = 0, -999, 0, 0, 0
    places = []

    #Iterate through the dictionary of WTsequence keys to determine the start position of the sequence
    for i, value in WTkey.items():
        try:
            #If a WT gene fragment is found within the sample sequence save the NT place
            place = seq.index(value)
            places.append(i-place)
            #This checks for indels if the keys are found at more than one start position
            if StartNT != (-999) and StartNT != (i - place):
                diff = 1
                if StartNT - (i - place) < 3:
                    indel = 1

            StartNT = i - place  #Finds the start NT first position is 1
            skip += 1
        except ValueError:
            None   

    #If indels detected, then use the most common start position for further analysis
    if diff == 1:
        StartNT = mypyfunctions.most_common(places)

    if StartNT < 0 and StartNT != -999: #Corrects for a negative alignent
        seq = seq[(0-StartNT):]
        StartNT = 0

    if indel == 1:
        seq = sequencefx.FindIndel(seq, StartNT)

    return StartNT

#Process the datafile
def LineReader (path, file, file2, WTMat, MutCor, FailSeq, FailC, AATable):
    count, TemplateLength, TemplateCount, Fhit, Rhit = 0, 0, 0, 0, 0
    EOFile = 0 #End of file

    os.chdir(path)
    Fdatafile = open(file)
    Rdatafile = open(file2)

    print('Opening data file')
    print('%s and %s files opened' % (file, file2))

    FseqID = Fdatafile.readline().rstrip()
    RseqID = Rdatafile.readline().rstrip()

    while EOFile == 0 and count < limit:
        Fseq =      Fdatafile.readline().upper().rstrip()
        FseqQ =     Fdatafile.readline().rstrip()
        FQscore =   Fdatafile.readline().rstrip()

        #Check the score to make sure the p < 0.05
        if len(Fseq) == len(FQscore):
            for i in range(0,len(Fseq)):
                if sequencefx.QscoreKey[FQscore[i]] <= 13:
                    if i == 0:
                        Fseq = 'N' + Fseq[1]
                    elif i == (len(Fseq) - 1):
                        Fseq = Fseq[:-1] + 'N'
                    else:
                        Fseq = Fseq[:(i-1)] + 'N' + Fseq[(i+1):]

        Rseq =      Rdatafile.readline().upper().rstrip()
        RseqQ =     Rdatafile.readline().rstrip()
        RQscore =   Rdatafile.readline().rstrip()

        #Check the score to make sure the p < 0.05
        if len(Rseq) == len(RQscore):
            for i in range(0,len(Rseq)):
                if sequencefx.QscoreKey[RQscore[i]] <= 13:
                    if i == 0:
                        Rseq = 'N' + Rseq[1]
                    elif i == (len(Rseq) - 1):
                        Rseq = Rseq[:-1] + 'N'
                    else:
                        Rseq = Rseq[:(i-1)] + 'N' + Rseq[(i+1):]

        #Finds the start location
        if FseqID[:43] == RseqID[:43]:
            Fhit = sequencefx.KeySearch(Fseq)
            if Fhit == -999:
                Fhit = sequencefx.KeySearch(proteinfx.RevComp(Fseq))
        
            Rhit = sequencefx.KeySearch(Rseq)
            if Rhit == -999:
                Rhit = sequencefx.KeySearch(proteinfx.RevComp(Rseq))

            if abs(Fhit - Rhit) < (len(Fseq) - sequencefx.seqOverlap):
                Oseq = sequencefx.Overlap(Fhit, Rhit, Fseq, Rseq)

            if Fhit != -999 and Rhit != -999:
                if Fhit < Rhit:
                    TemplateLength += len(Rseq) + Rhit - Fhit
                else:
                    TemplateLength += len(Rseq) + Fhit - Rhit
                TemplateCount += 1

            try:
                Fhit, Failseq, AATable, MutCor = sequencefx.QuickProcess(Fseq, WTMat, MutCor, FailSeq, AATable)
            except TypeError:
                None

            try:
                Rhit, Failseq, AATable, MutCor = sequencefx.QuickProcess(Rseq, WTMat, MutCor, FailSeq, AATable)
            except TypeError:
                None
    
        else:
            print('Sequences not the same!')
            Fhit, Failseq, AATable, MutCor = sequencefx.QuickProcess(Fseq, WTMat, MutCor, FailSeq, AATable)
            Rhit, Failseq, AATable, MutCor = sequencefx.QuickProcess(Rseq, WTMat, MutCor, FailSeq, AATable)

        if Fhit != 0:
            count += 1
        else:
            FailC += 1
    
        FseqID = Fdatafile.readline().rstrip()
        RseqID = Rdatafile.readline().rstrip()

        if FseqID == '' or RseqID == '':
            EOFile = 1              
   
        if (count + 1) % 10**6 == 0:
            print(count)

        count += 1
    Fdatafile.close()
    Rdatafile.close()
   
    return count, FailSeq, AATable, MutCor

def AAHandler(AATable):
    AATable[0][0] = AATable[0][0] + '-%s' % (SuccessNum)

    for j in range(1,len(AATable[0])):
        AATable[23][j] = max(mypyfunctions.column(AATable[1:23],j))
        AAind = mypyfunctions.column(AATable[:23],j).index(AATable[23][j])
        AATable[24][j] = Atable[AAind]
        AATable[25][j] = sum(mypyfunctions.column(AATable[1:23], j )) / 22 #Average mutations
        AATable[27][j] = sum(mypyfunctions.column(AATable[1:23], j )) # Total mutations
        if AATable[26][j] != 0:
            AATable[28][j] = AATable[27][j] / AATable[26][j]
            AATable[29][j] = AATable[23][j] / AATable[26][j] #Frequency of the largest mutation
        else:
            AATable[28][j] = 0
            AATable[29][j] = 0
    mypyfunctions.saver('AA%s.csv' % (i), AATable)
    return AATable

#Alternative function that uses a hash Key to find sequences in the WT gene
def KeyAlign(seq, BaseMut, MutCor, AATable):
    skip = 0  #Skips if seq is found, and the Key index
    MutCount = 0
    MutList = []
    StartNT = -999
    NPs = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4, '-':5}
    diff, indel = 0,0

    for i, value in WTkey.items():
        try:
            place = seq.index(value)
            #This checks for indels if the keys are found at more than one start position
            if StartNT != (-999) and StartNT != (i - place):
                diff = 1
                if StartNT - (i - place) < 3:
                    indel = 1

            StartNT = i - place  #Finds the start NT first position is 1
            skip += 1
        except ValueError:
            None         

    if StartNT < 0: #Corrects for a negative alignent
        seq = seq[(0-StartNT):]
        StartNT = 0

    if indel == 1: #Used if indel detected by varied start positions
        seq = sequencefx.FindIndel(seq, StartNT)

    if skip >= sequencefx.KeyCoverage:       
        holder = sequencefx.FindMut(seq, StartNT, BaseMut)
        MutCount = holder[0]
        MutList = holder[1]        

        if MutCount < sequencefx.MutLimit and len(MutList) > 0:
            for i in range(2, len(MutList),2):    
                Mpos = ((MutList[i] - MutList[0])*len(NPs)) - 3 + (NPs[MutList[i+1]])
                MutCor[MutList[0]][Mpos] += 1
           
            AATable = proteinfx.Trans(seq, StartNT, AATable)  
    
    return skip, AATable, MutCor
        
def FindMut(seq, StartNT, BaseMut):
    
    MutCount = 0  #Counts the number of mutations in a sequence
    MutList = [] #Make of list of the mutations contained in a single sequence.

    #Find any mutated bases
    Aligned = WTseq[StartNT:(StartNT + len(seq))]

    if len(Aligned) < len(seq):
        seq = seq[0:len(Aligned)]

    for i in range(0,len(seq)):        
        if (Aligned[i] != seq[i]) and seq[i] != '-':    
            try:
                BaseMut[sequencefx.NT.index(seq[i])][(StartNT + i)+1] += 1
            except ValueError:
                print(seq[i], (StartNT + i)+1)
            MutList.append((StartNT + i)+1) #Saves the position of the Mutation
            MutList.append(seq[i])  #Saves the type of mutation
            MutCount += 1

        #Used to get the coverage of the reads
        BaseMut[7][(StartNT + i)+1] += 1

    if MutCount > sequencefx.MutLimit:
        #indelC += 1
    
        MutList = []
        MutCount = 0

        #Undoes the coverage count if the sequence is found to be above the mutation tolerance
        for i in range(0,len(seq)):  
            BaseMut[7][(StartNT + i)+1] -= 1
    
    return MutCount, MutList

#Makes a report of the locations of mutations and how they correlate to one another.
def MutReport(WTMut, MutCor, filename):
    Tpose = []
    for i in range(0,len(WTMut[0])):
        Tpose.append(mypyfunctions.column(WTMut, i))

    mypyfunctions.saver('Report %s.csv' % (filename[:-9]), Tpose)
    mypyfunctions.saver('Correlation %s.csv' % (filename[:-9]), MutCor)
    return
