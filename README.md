# SampleWork\n
Sample programs for my graduate work.  See PDF for description of their use.\n

The files for the HTS experiment inlcude:\n
  HTS Align.py\n
  mypyfunctions.py\n
  proteinfx.py\n
  sequencefx.py\n

Reading the data and trimming the lower quality reads\n
This python script reads the fastq files one sequence at time and processes it before moving to the next. This allows the rapid scanning of incredibly large files without taxing the computer’s storage capacity with maintaining large amount of data.\n
 
Determination of alignment\n
The major working assumption of the program is the low error rate through the millions of sequences recovered and the coverage they provide. The workflow of the program is as follows:\n
1.      Read sequence from fastq file\n
2.      Remove nucleotides with quality scores with p < 0.05 (high probability of error)\n
3.      Compare the sequence to a set of keys made from the WT gene\n
4.      If the sequence aligns to 3 keys or more then the alignment is recorded\n
5.      The sequence is compared to the WT gene base by base saving any mutations in the nucleotide and residue sequence.\n
6.      If the number of mutations surpasses a threshold value (10 in this case), the program removes the sequence. Sequences with large errors have either InDels or something else wrong.\n
7.      If the mutations in a given sequence are below the threshold, the program records the location of the alignment and mutations found in the sequence.\n
8.      After processing all the sequences in a file, the program records the data in a file and moves onto the next fastq file.\n
9.      After processing all the fastq files, the program compiles the AA errors into a ‘Compiled AA.csv’ file where the user can analysis the most prevalent mutations for consensus. \n
