#Code by Zach Frye to align and compile residue mutation rates for HTS data.
#Portions of the output analysis are designed based on enrichment/depletion 
#calculations from (Whitehead et. al, PLOS one, 2015)

import itertools
import operator
import time
import os
import sys
import datetime

#Customized modules
import mypyfunctions
import proteinfx
import sequencefx

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
