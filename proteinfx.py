#Code by Zach Frye for basic protein functions

import itertools
import operator
import time
import os
import sys
import datetime

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
