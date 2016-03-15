#Code by Zach Frye to align and compile residue mutation rates for HTS data.
#Portions of the output analysis are designed based on enrichment/depletion 
#calculations from (Whitehead et. al, PLOS one, 2015)

import itertools
import operator
import time
import os
import sys
import datetime

#Contains commonly used functions for general python scripting
class mypyfunctions(object):
       
    #Returns column i from a table
    def column(matrix, i):
        return [row[i] for row in matrix]

    #Returns the most common list entry. Used to find major start position (Github)
    def most_common(L):
        # get an iterable of (item, iterable) pairs
        SL = sorted((x, i) for i, x in enumerate(L))
        # print 'SL:', SL
        groups = itertools.groupby(SL, key=operator.itemgetter(0))
        # auxiliary function to get "quality" for an item
        def _auxfun(g):
            item, iterable = g
            count = 0
            min_index = len(L)
            for _, where in iterable:
                count += 1
                min_index = min(min_index, where)
            return count, -min_index
        # pick the highest-count/earliest item
        return max(groups, key=_auxfun)[0]

    #Saves the Matrix file as csv
    def saver(label, m):
        os.chdir(OUTPUTTARGET)
        file = open(label, 'w')
    
        for i in m:
            buffer = ''
            for j in i:
                buffer = buffer + ',' + str(j) 
            file.write('%s\n' % buffer)
        file.close()
        return
    #Opens a file and returns the sequence
    def opener(path, label):
        #Remove this when ready
        os.chdir(path)
        file = open(label)
        #Discard the first Unless necessary for further data management
        line = file.readline()
        line = file.readline().upper()
        file.close()
        return line

#Class containing more general functions for handling sequences 
class proteinfx(object):
    
    # nucleotide to amino acid conversion
    amino_dict = {
        'AAA' : 'K', 'AAC' : 'N', 'AAG' : 'K', 'AAT' : 'N', 'ACA' : 'T', 'ACC' : 'T', 'ACG' : 'T', 'ACT' : 'T',
        'AGA' : 'R', 'AGC' : 'S', 'AGG' : 'R', 'AGT' : 'S', 'ATA' : 'I', 'ATC' : 'I', 'ATG' : 'M', 'ATT' : 'I', 

        'CAA' : 'Q', 'CAC' : 'H', 'CAG' : 'Q', 'CAT' : 'H', 'CCA' : 'P', 'CCC' : 'P', 'CCG' : 'P', 'CCT' : 'P',
        'CGA' : 'R', 'CGC' : 'R', 'CGG' : 'R', 'CGT' : 'R', 'CTA' : 'L', 'CTC' : 'L', 'CTG' : 'L', 'CTT' : 'L', 

        'GAA' : 'E', 'GAC' : 'D', 'GAG' : 'E', 'GAT' : 'D', 'GCA' : 'A', 'GCC' : 'A', 'GCG' : 'A', 'GCT' : 'A', 
        'GGA' : 'G', 'GGC' : 'G', 'GGG' : 'G', 'GGT' : 'G', 'GTA' : 'V', 'GTC' : 'V', 'GTG' : 'V', 'GTT' : 'V',

        'TAA' : '*', 'TAC' : 'Y', 'TAG' : '*', 'TAT' : 'Y', 'TCA' : 'S', 'TCC' : 'S', 'TCT' : 'S', 'TCG' : 'S',
        'TGA' : '*', 'TGC' : 'C', 'TGG' : 'W', 'TGT' : 'C', 'TTA' : 'L', 'TTC' : 'F', 'TTG' : 'L', 'TTT' : 'F',

        'GCN' : 'A', 'GGN' : 'G', 'GTN' : 'V', 'TCN' : 'S', 'CTN' : 'L', 'CCN' : 'P', 'CGN' : 'R', 'ACN' : 'T'}
    
    #Makes a list of translations in the first frames. 
    def trans_amino( dna_seq ):
        # declare empty list
        amino_acid_seq = ''
  
        for j in range( 0, int(len( dna_seq )/3) ):
            dna_trip = dna_seq[j*3:j*3+3]
            try:
                amino_acid_seq += proteinfx.amino_dict[dna_trip]
            except KeyError:
                amino_acid_seq += 'X'

        return amino_acid_seq

    #Saves a FASTA Matrix file for selected sequences
    def FastaSaver(label, m):
        os.chdir(OUTPUTTARGET)
        file = open(label, 'w')
    
        for i in m:
            file.write('>%s\n' % i[0])
            file.write('%s\n'  % i[1])
        file.close()
        return

    #Makes the Reverse Complement sequence (Pulled from GitHub)
    def RevComp(seq):
        reverse = list(seq[::-1])
        reverse_complement = []
    
        for i in range(0, len(reverse)):
            if reverse[i] == "A":
                reverse_complement.append("T")
            elif reverse[i] == "T":
                reverse_complement.append("A")
            elif reverse[i] == "C":
                reverse_complement.append("G")
            elif reverse[i] == "G":
                reverse_complement.append("C")
            elif reverse[i] == "N":
                reverse_complement.append("N")
            elif reverse[i] == "-":
                reverse_complement.append("-")

        return "".join(reverse_complement)

    #Makes the Qscore read Key to check the quality of the sequences at each Nt.
    def QsKey():
        key = {}
        Phred = '''!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI'''
        number = []
        [number.append(i) for i in range(1,42)]
        
        for i in range(0,41): 
            key[Phred[i]] = number[i]
        return key

    #Translator to determine the most prominant mutated amino acid
    def Trans(seq, NT1, AATable):
        Tseq = ''
        
        AA0 = int((NT1 + 2) / 3)   #Finds the first Amino acid in the sequence

        SeqFrame = (AA0 * 3) - NT1 + 1
    
        Tseq = proteinfx.trans_amino(seq[SeqFrame:])
        WTAAalign = WTAAseq[(AA0):(AA0+len(Tseq))]

        for i in range(0, len(Tseq)):
            if i < len(WTAAalign):
                AATable[26][(AA0 + i + 1)] += 1
                if Tseq[i] != WTAAalign[i]:   
                    #If the next AA is also mutated
                    if (i + 1) < len(Tseq) and (i + 1) < len(WTAAalign):
                        if Tseq[i + 1] != WTAAalign[i + 1]:
                            AATable[InvAtable[Tseq[(i + 1)]]][(AA0 + i + 1)] += 1
        return AATable

    
    #Load WT sequence and make matrix for mutations
    def WTReader (WTline, WTMat, MutCor, AATable):
        WTL = len(WTline)
        WTMat = [['WTseq'],['A'],['C'],['G'],['T'],['N'],['-'],['Count']]

        for i in range(0, WTL):
                WTMat[0].append(WTline[i])
        for i in range(1,8):
            for j in range(0, len(WTline)):
                WTMat[i].append(0)
    
        NTs = {0: 'A', 1:'C', 2:'G', 3:'T', 4:'N', 5:'-'}
        LNTs = len(NTs)

        #Make header for the matrix to track the mutation correlations
        MutCor.append(['Base','NT'])
        for i in range(0, LNTs * span):
            MutCor[0].append(NTs[i % LNTs])
        #Initiate the values in the MutCor matrix
        for i in range(1, WTL):
            MutCor.append([i, WTline[i-1]])
            for j in range(0, LNTs * span):
                MutCor[i].append(0)
    
        AATable = [['AAseq']]
        for i in WTAAseq:
            AATable[0].append(i)
        for i, value in Atable.items():
            AATable.append([])
            AATable[i].append(Atable[i])
            for j in range(0,len(WTAAseq)):
                AATable[i].append(0)
    
        return WTMat, MutCor, AATable

#Used to make the Alignment keys, read/process the data-stream, and compile the mutations 
class sequencefx(object):
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


#START PROGRAM


#Compiled Amino acid frequency table with the code for the individual rows
Atable = {1: 'F', 2: 'W', 3: 'Y', 4:'P', 5: 'M', 6: 'I', 7: 'L', 8: 'V', 9:'A', 10: 'G', 
            11: 'C', 12: 'S', 13: 'T', 14:'N', 15: 'Q', 16: 'D', 17: 'E', 18: 'H', 19:'K', 20: 'R', 
            21: 'X', 22: '*', 23: 'Max', 24: 'Mut', 25: 'Ave', 26: 'Read', 27: 'MutSum', 28: 'MutRate',
            29: 'MaxFreq' }

InvAtable = {v: k for k, v in Atable.items()}

#Data files from MiSeq125 sequencing run
FileList = ['File1-F.fastq.gz',       'File1-R.fastq.gz']

FILETARGET =   'Filepath /Data/'
OUTPUTTARGET =   'Filepath /Output/'      
os.chdir(FILETARGET)

span = 150 #Set this based on the expected length of the sequences
indelC = 0
AATable = []
limit = 10**4 #Limits the number of sequences analyzed in the database
SuccessSampling = 10**7

#Open the WT gene to extract the sequence then make the report matrices and hash key(s)
WTseq, WTMat, MutCor = '', [], []
WTseq = mypyfunctions.opener(FILETARGET, 'WTseq.txt')
WTAAseq = proteinfx.trans_amino(WTseq[1:])
WTMat, MutCor, AATable = proteinfx.WTReader(WTseq, WTMat, MutCor, AATable)

#Test sequences in pairs for each set of files
for h in range(0,len(FileList), 2):
    i =  FileList[h]
    iR = FileList[h+1]
    FailC, FailSeq = 0, [] #List of failed sequences
    
    #Track the processing time for the files
    StartT = time.time()
    WTMat, MutCor, AATable = proteinfx.WTReader(WTseq, WTMat, MutCor, AATable)  #Resets the records
    SuccessNum, Failseq, AATable, MutCor = sequencefx.LineReader(FILETARGET, i, iR, WTMat, MutCor, FailSeq, FailC, AATable)
    
    #Generate the reports
    sequencefx.MutReport(WTMat, MutCor, i)
    if len(FailSeq) > 0:
        proteinfx.FastaSaver('%s.csv' % (i), FailSeq)
    sequencefx.AAHandler(AATable)
    
    EndT = time.time()
    print('Finished all %s sequences in %s seconds.' % (SuccessNum, (EndT - StartT)))
    print('Found %s failed alignments and %s indel suspects.\n' % (len(FailSeq), indelC))

#Matrix for the compiled error rates
Comp_Error = []
[Comp_Error.append([i]) for i in WTMat[0] ]
Comp_AA = []

#Compiles the Data after the analysis
for i in os.listdir():
    #Process the Error Reports
    if i.endswith('csv') and i.startswith('Report'):
        file = open(i)

        line = file.readline()
        parts = line.split(',')
        Comp_Error[0].append('%s-%s' % (i.strip(), parts[8].strip())) 
        Comp_Error[0].append('%s-Ratio' % (i.strip())) 

        for j in range(1,len(WTMat[0])):
            line = file.readline()
            parts = line.split(',')

            sum = 0
            for k in range(2,9):
                sum +=  int(parts[k].strip())

            try:
                Comp_Error[j].append(parts[k].strip())
                Comp_Error[j].append( (sum - int(parts[8].strip())) / int(parts[8].strip())) 
            except ZeroDivisionError:
                None 
        file.close()

    #Process the AA Mutation reports
    if i.endswith('csv') and i.startswith('AA'):
        file = open(i)

        for j in range(0,len(Atable)+1):
            line = file.readline().rstrip()
            parts = line.split(',')
            parts[0] = i
            
            if j > 22:
                Comp_AA.append(parts)
        file.close()

mypyfunctions.saver('Compiled Errors.csv', Comp_Error)
mypyfunctions.saver('Compiled AA.csv', Comp_AA)

print('Finished!')
