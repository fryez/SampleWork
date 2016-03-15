# SampleWork
Sample programs for my graduate work.  See PDF for description of their use.

The files for the HTS experiment inlcude:
  HTS Align.py
  mypyfunctions.py
  proteinfx.py
  sequencefx.py
  
  Please use Alternative_AlignFile.py if the separated file and modules above don't work.

Reading the data and trimming the lower quality reads
This python script reads the fastq files one sequence at time and processes it before moving to the next. This allows the rapid scanning of incredibly large files without taxing the computer’s storage capacity with maintaining large amount of data.
 
Determination of alignment
The major working assumption of the program is the low error rate through the millions of sequences recovered and the coverage they provide. The workflow of the program is as follows:
1.      Read sequence from fastq file
2.      Remove nucleotides with quality scores with p < 0.05 (high probability of error)
3.      Compare the sequence to a set of keys made from the WT gene
4.      If the sequence aligns to 3 keys or more then the alignment is recorded
5.      The sequence is compared to the WT gene base by base saving any mutations in the nucleotide and residue sequence.
6.      If the number of mutations surpasses a threshold value (10 in this case), the program removes the sequence. Sequences with large errors have either InDels or something else wrong.
7.      If the mutations in a given sequence are below the threshold, the program records the location of the alignment and mutations found in the sequence.
8.      After processing all the sequences in a file, the program records the data in a file and moves onto the next fastq file.
9.      After processing all the fastq files, the program compiles the AA errors into a ‘Compiled AA.csv’ file where the user can analysis the most prevalent mutations for consensus. 
