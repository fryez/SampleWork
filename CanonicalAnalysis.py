#This program generates the RMSD tables based on the alignments of the TCR
#CDRs based on sequences in the mastertable. Need to add a combinatorial
#alignment scheme to align different residue selections of the loops

#From pymol import cmd, CmdException, stored;
from operator import itemgetter;
import os;
import string;
import datetime;
import re;
import types;
import random;
import csv;
 
FILETARGET = 'FILEPATH/';
PDBTARGET = 'FILEPATH/PDB/'; 

#Opens the table of crystal structures data to compare and the sequences to find 
#and stores it in seqlist. This list also includes other information.
#Could use this list to make a dictionary with the masterlist TCR information
def CDRLIST(): 
    file = open( 'TCRMasterlist.csv') #Spreadsheet with all the TCR information on it for pulling CDRs
    #Saves the TCRMasterlist header in case the spreadsheet gets manipulated later
    line = file.readline()
    parts = line.strip().split(',')
    for i in range(0, len(parts)):
        header[parts[i]] = i 
    #Compiles a list for the TCR properties
    line = file.readline()
    count = 0
    while (line != ''):
        seqlist.append([])
        parts = line.strip().split(',')
        for part in parts:
            seqlist[count].append( part )
        count += 1
        line = file.readline()
    file.close()
    return count;

#Used to interface with the 3-letter and 1-letter AA names
#Used to determine the various definitions of atoms read from the PDB files.
def OneLetter():    
    d, R = {}, {}
    file = open("oneletter.csv")
    for line in file:
        (key, val) = line.strip().split(',')
        d[key] = val
    R = {'A': 'ALA', 'C': 'CYS',  'D': 'ASP',  'E': 'GLU',  'F': 'PHE',  
         'G': 'GLY', 'H': 'HIS',  'I': 'ILE',  'K': 'LYS',  'L': 'LEU',
         'M': 'MET', 'N': 'ASN',  'P': 'PRO',  'Q': 'GLN',  'R': 'ARG',  
         'S': 'SER', 'T': 'THR',  'V': 'VAL',  'W': 'TRP',  'Y': 'TYR'}

    return d, R

#Returns all the items in a column of a matrix
def column(matrix, i):
    return [row[i] for row in matrix]

#Formats the first row and column of the output file with the structurenames
#First row and column are the pdb names according to how the comparisons are made
def recorder ():
    record.append([])
    record[0].append('CDRs')
    for i in range (0,count):
        record[0].append(seqlist[(count - i-1)][1])
    for i in range (1,(count+1)):
        record.append([])
        record[i].append( seqlist[(i - 1)][1] )

#Makes a list of all the PDB files in the target dir.
def checkParams(needle,haystack,selName,het,firstOnly):
	#This is just a helper function for checking the user input
	
	# check Needle
	if len(needle)==0 or type(needle)!=types.StringType:
		print ("Error: Please provide a string 'needle' to search for.")
		print ("Error: For help type 'help motifFinder'.")
		return False

	# check Haystack
	if len(haystack)==0 or type(haystack)!=types.StringType:
		print ("Error: Please provide valid PyMOL object or selection name")
		print ("Error: in which to search.")
		print ("Error: For help type 'help motifFinder'.")
		return False

	# check het
	try:
		het = bool(int(het))
	except ValueError:
		print ("Error: The 'het' parameter was not 0 or 1.")
		return False
	if type(selName)!=types.StringType:
		print ("Error: selName was not a string.")
		return False
	return True

#Generates separate objects for each peptide chain in a PDB file
def split_chains(selection='(all)', prefix=None):
    count = 0
    models = cmd.get_object_list('(' + selection + ')')

    for model in models:
        chains = cmd.get_chains('(%s) and model %s' % (selection, model))
        for chain in chains:
            if chain == '':
                chain = "''"
            count += 1
            if not prefix:
                name = '%s%s' % (model, count)

            else:
                name = '%s%04d' % (prefix, count)
            cmd.create(name, '(%s) and model %s and chain %s' % (selection, model, chain))
#        cmd.disable(model)
    return count;
cmd.extend('split_chains', split_chains)

#Renumbers the residues in a molecule object
def renumber(selection='all', start=1, startsele=None, quiet=1):
    start, quiet = int(start), int(quiet)
    model = cmd.get_model(selection)
    cmd.iterate(selection, 'atom_it.next().model = model',
            space={'atom_it': iter(model.atom)})
    if startsele is not None:
        startidx = cmd.index('first (' + startsele + ')')[0]
        for atom in model.atom:
            if (atom.model, atom.index) == startidx:
                startatom = atom
                break
        else:
            print(' Error: startsele not in selection')
            raise CmdException
    else:
        startatom = model.atom[0]
    for atom in model.atom:
        atom.adjacent = []
        atom.visited = False
    for bond in model.bond:
        atoms = [model.atom[i] for i in bond.index]
        atoms[0].adjacent.append(atoms[1])
        atoms[1].adjacent.append(atoms[0])
    minmax = [start, start]
    def traverse(atom, resi):
        atom.resi = resi
        atom.visited = True
        for other in atom.adjacent:
            if other.visited:
                continue
            if (atom.name, other.name) in [('C','N'), ("O3'", 'P')]:
                minmax[1] = resi+1
                traverse(other, resi+1)
            elif (atom.name, other.name) in [('N','C'), ('P', "O3'")]:
                minmax[0] = resi-1
                traverse(other, resi-1)
            elif (atom.name, other.name) not in [('SG', 'SG')]:
                traverse(other, resi)
    traverse(startatom, start)
    cmd.alter(selection, 'resi = atom_it.next().resi',
            space={'atom_it': iter(model.atom)})
    if not quiet:
        print(' Renumber: range (%d to %d)' % tuple(minmax))
cmd.extend('renumber', renumber)

#Secondary renumber function for protein chains
def renumber2(selection='all', start=1, startsele=None, quiet=1):
    start, quiet = int(start), int(quiet)
    atom = []
    resi = start - 1
    cmd.iterate('%s' % selection, 'atom.append([chain, resi, name])', space = locals())
    for i in atom:
        if i[2] == 'N':
            resi += 1
        cmd.alter('%s and (chain %s and resi %s)' % (selection, i[0], i[1]), 'resi = %s' % resi)
    return

#Sorts the table based on the TCR names before doing the comparisons
def sorter(rec, item):
    ind = 0
    for i in range (0,count):
        if (rec[i][0] == item):
            ind = i
            x = count + 1
    rec = sorted(rec, key = itemgetter(ind))

    i = 1
    while (i <= count):
        if (rec[i][0] == 'CDRs'):
            buffer = rec[i]
            for x in range ((i), 0, -1):
                rec[x] = rec[x-1]
            rec[0] = buffer
            i = count
        i +=1
    rec = zip(*rec)
    for x in range(0,len(rec)):
        rec[x] = list(rec[x])
    rec = sorted(rec, key = itemgetter(1))
    i = 1
    while (i <= count):
        if (rec[i][0] == 'CDRs'):
            buffer = rec[i]
            for x in range ((i), 0, -1):
                rec[x] = rec[x-1]
            rec[0] = buffer
            i = count
        i +=1
    printer(rec)

def printer(m):
    print('%s' % (m))
    for x in range(0,len(m)):
        print(m[x])
    print('')

#Fetches the structure file and isolates only the mainchain backbone atoms for pdbs
def opener(x):
    cmd.select('mainchain', None) #Make null selection for the mainchain
    
    os.chdir(PDBTARGET)        
    cmd.fetch ('%s' % (seqlist[x][1]))  #Open the crystal file
    os.chdir('%s%s' % (FILETARGET, outDir))

    cmd.remove('!(alt "+A")') #Removes the alternative atomic positions
    cmd.select('mainchain', 'n. N+CA+C+O in %s' % (seqlist[x][1]))
    cmd.select('side', '%s be. 0 of mainchain' % (seqlist[x][1]))
    cmd.remove('side')
    cmd.delete('side and mainchain')
    cmd.hide('nonbonded')
    cmd.hide('all')
    
#Finds specific AA sequence in a protein chain (from GitHub)
def findseq(needle, haystack, selName=None, het=0, firstOnly=1):
    # set the name of the selection to return.
    if selName == None:
        rSelName = "foundSeq" + str(random.randint(0,32000))
        selName = rSelName
    elif selName == "sele":
        rSelName = "sele"
    else:
        rSelName = selName

    # input checking
    if not checkParams(needle, haystack, selName, het, firstOnly=1):
        print ("There was an error with a parameter.  Please see")
        print ("the above error message for how to fix it.")
        return None

    # remove hetero atoms (waters/ligands/etc) from consideration
    cmd.select("__h", "br. " + haystack + " and not het")
        
    # get the AAs in the haystack
    aaDict = { 'aaList' : [] }
    cmd.iterate("(name ca) and __h","aaList.append((resi,resn))",space=aaDict)
    #modified the int(x[0]) for str(x[0])
    IDs = map( lambda x: str(x[0]), aaDict['aaList'] )
    AAs = ''.join(map( lambda x: one_letter[x[1]], aaDict['aaList'] ))

    reNeedle = re.compile( needle.upper() )
    it = reNeedle.finditer( AAs )
    if it != None:
        cmd.select( rSelName, 'None') # make an empty selection to which we add residues
        for i in it:
            (start,stop) = i.span()
            cmd.select( rSelName, rSelName + " __h and i. " + str(IDs[start]) + "-" + str(IDs[stop-1]))
            break
        cmd.delete("__h")
    else:
        print ('Sequence not found!')
        rSelName = None
    return rSelName
cmd.extend("findseq", findseq )

#Fills out the rest of the alignment matrix making it diagonally symmetric
def process(trecord):
    x = 1
    while x <= (count - 1): #Fills in the rest of the matrix along the diagonal
        y = 1
        while y <= x:
            trecord[(count+1-y)].append(trecord[(count-x)][y])
            y += 1
        x += 1
    buffer = []
    trecord = list(zip(*trecord)[::-1])
    for i in range((len(trecord)-1),0,-1):
        buffer = trecord[i]
        trecord[i] = trecord[i-1]
        trecord[i - 1] = buffer
    for i in range(0,len(trecord)):
        trecord[i] = list(trecord[i])  
    

    return trecord

#Saves the alignment tables into individual loops in the new folder
def saver(rec,label):
    os.chdir('%s%s' % (FILETARGET, outDir))
    outFile = open('%s.csv' % (label), 'w');
    for i in rec:
        try:
            i.append('\n')
            outFile.write(','.join(i));
        except AttributeError:
            outFile.writelines('%s\n' % (i))
    outFile.close();
    os.chdir(FILETARGET);
    return
#Makes the output file for the PDB
def output(cname, ave, std, AAseq):
    outFile = open('%s.pdb' % (cname), 'w')
    recordCounter = 0;
    recordNumber = 1;
    counter = 0;
    residueSequenceNumber = 0;
    for i in range(0,len(ave)):
        atom = None
        symbol = None
        if recordCounter == 0:
            atom = 'N';
            symbol = 'N';
            resName = Rone_letter[AAseq[residueSequenceNumber]]
            residueSequenceNumber += 1;
        elif recordCounter == 1:
            atom = 'CA';
            symbol = 'C';
        elif recordCounter == 2:
            atom = 'C';
            symbol = 'C';
        else:
            atom = 'O';
            symbol = 'O';
        temperature = std[i]
        outFile.write(recordline(recordNumber, atom, ave[i], symbol, residueSequenceNumber, resName, temperature));        
        recordNumber += 1;
        recordCounter += 1;
        if recordCounter > 3:
            recordCounter = 0;
    outFile.close();
#Writes the PDB file for the ouput using bfactor as the stdeviation
def recordline(j, atom, coord, symbol, residueSequenceNumber, resName, temperature=0.0):
    residue= resName
    chain='A'
    occupancy=1
    segmentIdentifier='A1'
    charge=''

    buffer = [];
    # Record name
    buffer.append('ATOM'.ljust(6, ' '));
    # Atom serial number
    buffer.append(str(j).rjust(5, ' '));
    # Space
    buffer.append(' ');
    # Atom name
    buffer.append(atom.ljust(4, ' '));
    # Alternate location indicator
    buffer.append(' ');
    # Residue name
    buffer.append(residue.ljust(3, ' '));
    # Space
    buffer.append(' ');
    # Chain identifier
    buffer.append(chain);
    # Reside sequence number
    buffer.append(str(residueSequenceNumber).rjust(4, ' '));
    # Code for insertion of residues
    buffer.append(' ');
    # Space
    buffer.append(''.ljust(3, ' '));
    # Orthogonal coords for X
    buffer.append(('%.3f' % coord[0]).rjust(8, ' '));
    # Orthogonal coords for Y
    buffer.append(('%.3f' % coord[1]).rjust(8, ' '));
    # Orthogonal coords for Z
    buffer.append(('%.3f' % coord[2]).rjust(8, ' '));
    # Occupancy
    buffer.append(('%.2f' % occupancy).rjust(6, ' '));
    # Temperature factor
    buffer.append(('%.2f' % temperature).rjust(6, ' '));
    # Space
    buffer.append(''.ljust(6, ' '));
    # Segment identifier
    buffer.append(segmentIdentifier.ljust(4, ' '));
    # Element symbol
    buffer.append(symbol.rjust(2, ' '));
    # Charge on atom
    buffer.append(charge.ljust(2, ' '));
    buffer.append('\n');
    return ''.join(buffer);

#This loads the Chain ID table
def ch_id_loader():
    data = []

    os.chdir(FILETARGET)
    file = open( 'Chain_ID.csv') #Spreadsheet with all the TCR information on it for pulling CDRs
    #Saves the TCRMasterlist header in case the spreadsheet gets manipulated later
    line = file.readline()
    while line != '':
        parts = line.strip().split(',')
        buffer = []
        for part in parts:
            buffer.append(part)

        data.append(buffer)

        line = file.readline()
            
    return data

#This function finds the sequence of interest and isolates it as a selection
#Both the CDR loop and selection exist as separate entities. 
#This is for possible expansion for the selections to inculde multiple CDRs
def isolater(TCRnum, LoopNum, NonCon, ch_id):
    ch_holder, mhc_holder = [], []
    del_list = [] #List of chains to delete
    TCR = seqlist[TCRnum][1]
    CDR = seqlist[TCRnum][LoopNum]
    x = split_chains(TCR) #Sets the total number of chains to search through.
    TCRaList, TCRbList,  = [], []

    iC = ''
    i = 1
    chct = 1 #Counts the number of Identical CDRs with the sequence.
    j = 0
    while (i <= x):
        #Find the Sequence in the chain
        iC = findseq( '%s' % (CDR), '%s%s' % ( TCR, i ), '%sCDR%s-%s' % (TCR, LoopNum, chct), 0, 1)

        #If the CDR sequence is there with all the appropriate atoms in the backbone initiate 
        if (cmd.select('%s' % iC) != 0) and (cmd.select('%s and name C+CA+N+O' % iC) == len(CDR)*4):
            cmd.extract('%sCDR%s-%s' % (TCR, LoopNum, chct), '%sCDR%s-%s' % (TCR, LoopNum, chct))
            #cmd.select('rSI', '%s%s and !(br. %sCDR%s-%s)' % ( TCR, i, TCR, LoopNum, chct ))    
            cmd.remove('%s%s' % ( TCR, i ))
            cmd.delete('%s%s' % ( TCR, i ))      
            #renumber2('%sCDR%s-%s' % (TCR, LoopNum, chct))

            if LoopNum == 18:
                TCRaList.append(cmd.get_chains('%sCDR%s-%s' % (TCR, LoopNum, chct))[0])
            if LoopNum == 21:
                TCRbList.append(cmd.get_chains('%sCDR%s-%s' % (TCR, LoopNum, chct))[0])
            
                #Detect the pMHC if analyzing CDR3beta and the PDB is bound
            if LoopNum == 23 and 'Bound' == seqlist[TCRnum][header['B/UB']]:          
                #Find the center of the CDR3beta loop 
                M_iter = []
                cmd.iterate('%sCDR%s-%s' % ( seqlist[TCRnum][1], LoopNum, chct), 'M_iter.append([chain ,resi, name, b])', space = locals())
                center = M_iter[int(len(M_iter)/2)]
                cmd.select('bC104_selector','resi 1 and %sCDR%s-%s' % ( seqlist[TCRnum][1], LoopNum, chct))
                cmd.select('alphachain', '(bc. (all w. 12 of bC104_selector)) and (!%sCDR%s-%s and (!%s and !(chain M_iter[0][0])))' % (seqlist[TCRnum][1], LoopNum, chct, TCR))

                #Detect proximity to each chain since they are separated.
                try:
                    a_chain = cmd.get_chains('br. ((%s and !(%s and chain %s)) w. 5 of ((chain %s and resi %s) and (name %s and %sCDR%s-%s)))'
                                             % (seqlist[TCRnum][1], seqlist[TCRnum][1], M_iter[0][0], M_iter[-1][0], M_iter[-1][1], 
                                                M_iter[-1][2], seqlist[TCRnum][1], LoopNum, chct))[0]
                except IndexError:
                    a_chain = chr(ord(center[0]) - 1).upper()
                cmd.select('alpha', 'br. ((%s and !(%s and chain %s)) w. 8 of ((chain %s and resi %s) and (name %s and %sCDR%s-%s)))' 
                                            % (seqlist[TCRnum][1], seqlist[TCRnum][1], M_iter[0][0], M_iter[-1][0], M_iter[-1][1], 
                                            M_iter[-1][2], seqlist[TCRnum][1], LoopNum, chct))
                #Determine the molecules nearest to the CDR3beta loop
                cmd.select('%s-MHC' % seqlist[TCRnum][1], '(bc. (%s w. 14 of (resi %s and %sCDR%s-%s))) and !(alphachain)' 
                            % (seqlist[TCRnum][1], center[1], seqlist[TCRnum][1], LoopNum, chct))
                cmd.select('%s-MHC' % seqlist[TCRnum][1], '%s-MHC and (%s and !chain %s+%s)' % (seqlist[TCRnum][1], seqlist[TCRnum][1], center[0], a_chain))
                        
                for gc in cmd.get_chains('%s-MHC' % seqlist[TCRnum][1]):
                    if gc != '':
                        mhc_holder.append(gc)
                
            chct += 1
            j = x
        elif LoopNum != 23:
            cmd.delete('%sCDR%s-%s' % (TCR, LoopNum, chct))
            cmd.remove('%s%s' % (TCR, i))
            cmd.delete('%s%s' % (TCR, i))
        else:
            del_list.append([i, chct])
        i += 1

    chct -= 1

    if LoopNum == 23 and 'Bound' == seqlist[TCRnum][header['B/UB']]:
        for i in del_list:
            #cmd.delete('%sCDR%s-%s' % (TCR, LoopNum, i[1]))
            cmd.remove('%s%s' % (TCR, i[0]))
            cmd.delete('%s%s' % (TCR, i[0]))
            cmd.delete('%s-MHC' % (TCR))
       
    if (j == 0):
        print('Wrong CDR sequence in %s-%s' % (TCR, loop[CDRLoop]))
        NonCon.append('%s-%s' % (loop[CDRLoop], TCR))

    cmd.delete('sele')
    
    #if 'Bound' == seqlist[TCRnum][header['B/UB']]:
    if len(TCRaList) != 0 and LoopNum == 18:
        ch_id.append(['%s' % TCR, str(len(TCRaList)), '%s' % (','.join(TCRaList))])
    elif LoopNum == 18:
        ch_id.append(['%s' % TCR, '0'])

    if len(TCRbList) != 0 and LoopNum == 21:
        try:
            TCRindex = column(ch_id, 0).index(TCR)
            ch_id[TCRindex].append(str(len(TCRbList)))
            ch_id[TCRindex].append('%s' % (','.join(TCRbList)))
        except ValueError:
            print('Alpha missing for %s' % TCR)
            ch_id.append(['%s' % TCR, '0', str(len(TCRbList)), '%s' % (','.join(TCRbList))])
    elif LoopNum == 21:
        try:
            TCRindex = column(ch_id, 0).index(TCR)
            ch_id[TCRindex].append('0')
        except ValueError:
            print('error finding the beta chain for %s' % TCR)

    if len(mhc_holder) != 0 and LoopNum == 23:
        try:
            TCRindex = column(ch_id, 0).index(TCR)
            ch_id[TCRindex].append(str(len(mhc_holder)))
            ch_id[TCRindex].append('%s' % (','.join(mhc_holder)))
        except ValueError:
            print(TCR, 'Both chains missing?')
            ch_id.append(['%s (MHC only)' % TCR, str(len(mhc_holder)), '%s' % (','.join(mhc_holder))])
    elif LoopNum == 23:
        try:
            TCRindex = column(ch_id, 0).index(TCR)
            ch_id[TCRindex].append('0')
        except ValueError:
            print('error finding the pMHC chain for %s' % TCR)
            ch_id.append(['%s' % TCR, '0','0', '0'])
    
    cmd.remove('%s' % (TCR))
    cmd.delete('%s' % (TCR))
    
    return chct, NonCon, ch_id

#Used to analyse the identical loops in the PDB file
def IDLooper(chct, i):
    maximum = 0.0
    div = 1
    A_iter, B_iter, dist_iter, Ave_iter, Aveb_iter, xyz = [],[],[],[],[],[]
        
    for k in range (1,chct): #Compares all the identical sequence CDRs for variations
        A_iter = []
        cmd.iterate('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, k), 'A_iter.append([chain ,resi, name, b])', space = locals())
        
        if k == 1: #Only does this to initialize the Average holder on the first pass
            CDRL = len(A_iter) #The total number of atoms in the CDR
            xyz = cmd.get_model('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, k), 1).get_coord_list()
            for m in range(0,CDRL):
                dist_iter.append(0)
                Ave_iter.append(xyz[m])
                Aveb_iter.append(A_iter[m][3])
        
        for l in range ((k + 1), (chct + 1)): #selected the mobile state for alignment and comparison
            fitter = cmd.align('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, l), '%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, k ),100.0, 5, -100, quiet = 1)
            if fitter[0] > maximum:
                maximum = fitter[0]
                
            B_iter = []
            cmd.iterate('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, l), 'B_iter.append([chain ,resi, name, b])', space = locals())
            xyz = cmd.get_model('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, k), 1).get_coord_list()

            for m in range(0,CDRL):
                dist_iter[m] += cmd.get_distance('(%sCDR%s-%s and chain %s) and (resi %s and name %s)' % ( seqlist[i][1], CDRLoop, k, A_iter[m][0],A_iter[m][1],A_iter[m][2]),
                                                '(%sCDR%s-%s and chain %s) and (resi %s and name %s)' % ( seqlist[i][1], CDRLoop, l, B_iter[m][0],B_iter[m][1],B_iter[m][2]))
                for n in range(0,3):
                    Ave_iter[m][n] += xyz[m][n]
                Aveb_iter[m] += B_iter[m][3]
            div += 1
            
    for m in range(0,CDRL):
        dist_iter[m] = dist_iter[m]/div
        
    if maximum > IDcutoff:
        for k in range (1,(chct + 1)):
            cmd.show('sticks', '%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, k))

            A_iter = []
            cmd.iterate('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, k), 'A_iter.append([chain ,resi, name])', space = locals())

            for l in range(0,len(A_iter)):
                cmd.alter('(%sCDR%s-%s and chain %s) and (resi %s and name %s)' % ( seqlist[i][1], CDRLoop, k, A_iter[l][0],A_iter[l][1],A_iter[l][2]), 'b = %s' % dist_iter[l])
            cmd.spectrum ('b', 'green_yellow_red', '%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, k))

        #print('printing variable CDRs for %s with a %s max' % (seqlist[i][1], maximum))
        cmd.orient('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, 1))
        cmd.zoom('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, 1), complete=1)
        cmd.ray(1000,1000)
        os.chdir('%s%s' % (FILETARGET, outDir))
        cmd.set('ray_shadow', 'off')
        cmd.bg_color('white')
        cmd.png('%s %s' % ( loop[CDRLoop], seqlist[i][1]))
        os.chdir(PDBTARGET)
    elif maximum <= IDcutoff:
        for m in range(0, CDRL):
            for n in range(0,3):
                Ave_iter[m][n] = Ave_iter[m][n] / div
        os.chdir('%s%s' % (FILETARGET, outDir))
        output('%s-ave-%s' % ( loop[CDRLoop], seqlist[i][1]), Ave_iter, Aveb_iter, seqlist[i][CDRLoop])

    return maximum
    
#Opens and isolates CDR loops
def CompileLoops(IDloop, Nonconverging, ch_id):
    os.chdir(PDBTARGET);
    for i in range (0,count):
        if (seqlist[i][CDRLoop] != ''):
            opener(i) #Opens the individual files
            chct = isolater(i, CDRLoop, Nonconverging, ch_id)[0] #Isolates the CDR loops into their own objects
            
            #cmd.remove('all') #Just for Chain ID
            #cmd.delete('all') #Just for Chain ID

            #''' #Just for Chain ID
            #If there are more than one identical loop then process it further
            if chct == 1:
                cmd.set_name('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, 1), '%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, 1))
            elif chct > 1:
                Max_loop = IDLooper(chct, i)
            
                if Max_loop <= IDcutoff: #If the identical loops are pretty much the same (within the 1A cutoff) then find the average structure
                    #print('TCR %s converges and the average is used' % seqlist[i][1])
                    for j in range(1,(chct+1)):
                        cmd.remove('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, j))
                        cmd.delete('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, j))
                        cmd.delete('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, j))
                        #print('removed %sCDR%s-%s' % ( seqlist[i][1], CDRLoop, j))

                    os.chdir('%s%s' % (FILETARGET, outDir))
                    cmd.load('%s-ave-%s.pdb' % ( loop[CDRLoop], seqlist[i][1]))
                    cmd.set_name('%s-ave-%s' % (loop[CDRLoop], seqlist[i][1]), '%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, 1))
                    #print('set name to %sCDR%s-%s' % ( seqlist[i][1], CDRLoop, 1))
                else:
                    print('TCR %s does not have a converging loop' % seqlist[i][1])
                    Nonconverging.append('%s-%s' % (loop[CDRLoop], seqlist[i][1]))    
                    for j in range(1,(chct+1)):
                        cmd.remove('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, j))
                        cmd.delete('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, j))
                        cmd.delete('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, j))
                IDloop.append([loop[CDRLoop], seqlist[i][1], str(Max_loop)])
            #''' #Just for Chain ID
    
    os.chdir(FILETARGET);
    return Nonconverging, ch_id

#To clear out all the information from the previous program Run    
cmd.delete('all')
cmd.reinitialize()
#Make the target filename fo the record output
os.chdir(FILETARGET);
outDir = datetime.datetime.now().strftime('%Y-%m-%dT%H-%M-%S');
os.mkdir(outDir);  #Only to be used if there will be a folder from the processes
one_letter, Rone_letter = OneLetter() #sets up the dictionaries for res/atoms

count = 0
CDRLoop = 18 #1 to 6 is a1,a2,a3,b1,b2,b3 respectively
cutoff = 1.0
IDcutoff = 1.0
seqlist, header, record, newCDR = [], {}, [], []
loop = {18: 'A1', 19: 'A2', 20: 'A3', 21: 'B1', 22: 'B2', 23: 'B3'}
LoopResi = {1: '23-41', 2: '53-68', 3: '104-118'}

alpha_chain = ''

#Establishes the list of properties to use in the program from the masterlist
count = CDRLIST()
count = 20  #delete this when ready 249 is the total

print(len(seqlist), count)

IDloop, NonConverge, ch_id = [], [], []
ch_idIMGT = ch_id_loader()

#Iterate through all the TCRs and all 6 CDRs
for x in range (0,6):
    poorfit = 0
    CDRLoop = 18 + x
    record = []
    recorder() #sets up the title row and column of the list
    CompileLoops(IDloop, NonConverge, ch_id)
    
    #''' #Just for Chain ID
    print(x, NonConverge)
    i = 0
    while (i <= count-1):
        j = count - 1


        try:    #If the i-loop is on the nonconverging list then skip it
            NonConverge.index('%s-%s' % (loop[CDRLoop], seqlist[i][1]))
            skip = 0
        except ValueError:
            skip = 1

        if ((seqlist[i][CDRLoop] != '') and (skip == 1)):
            while (j > i):
                try:
                    NonConverge.index('%s-%s' % (loop[CDRLoop], seqlist[j][1]))
                    skip = 0
                except ValueError:
                    skip = 1

                if ((seqlist[j][CDRLoop] != '') and (skip == 1)):
                    tL = cmd.select('%sCDR%s-%s' % (seqlist[i][1], CDRLoop, 1)) #Finds the CDR length of the target selection
                    mL = cmd.select('%sCDR%s-%s' % (seqlist[j][1], CDRLoop, 1)) #Finds the CDR length of the mobile selection
                    #if ((len(seqlist[j][CDRLoop])*1.0 != mL/4.0) or (len(seqlist[i][CDRLoop]) * 1.0 != tL/4.0)) or tL != mL:
                        #print(len(seqlist[j][CDRLoop]), mL/4.0, len(seqlist[i][CDRLoop]), tL/4.0)
                    if len(seqlist[j][CDRLoop])*1.0 != mL/4.0:
                        print('Mobile selection off')
                    if len(seqlist[i][CDRLoop])*1.0 != tL/4.0:
                        print('Target selection off')    
                
                if ((seqlist[j][CDRLoop] != '') and (len(seqlist[j][CDRLoop]) == len(seqlist[i][CDRLoop]))) and (skip == 1):
                    fitter = cmd.align('%sCDR%s-%s' % (seqlist[i][1], CDRLoop, 1), '%sCDR%s-%s' % (seqlist[j][1], CDRLoop, 1),100.0, 5, -100, quiet = 1)
                    minimum = min(cmd.select('%sCDR%s-%s' % (seqlist[i][1], CDRLoop, 1)), cmd.select('%sCDR%s-%s' % (seqlist[j][1], CDRLoop, 1)))/4
                    #For perfect fits of the entire loops
                    if(fitter[6] == minimum):
                        record[(i+1)].append(str(fitter[0])[:5])
                    elif(abs(fitter[6] - minimum) <= 4):
                        record[(i+1)].append('LLL%s' %(str(fitter[0])[:5]))                                             
                        poorfit += 1
                        
                    else: #Triggered if the loops are different lengths by over 4 residues
                        record[(i+1)].append('poor')

                        cmd.hide('everything', 'all')
                        cmd.show('sticks', '%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, 1))
                        cmd.color('red', '%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, 1))
                        cmd.show('sticks', '%sCDR%s-%s' % ( seqlist[j][1], CDRLoop, 1))
                        cmd.color('blue', '%sCDR%s-%s' % ( seqlist[j][1], CDRLoop, 1))
                        cmd.orient('%sCDR%s-%s' % ( seqlist[i][1], CDRLoop, 1))

                elif ((seqlist[j][CDRLoop] != '') and (skip == 1)):
                    fitter = cmd.align('%sCDR%s-%s' % (seqlist[i][1], CDRLoop, 1), '%sCDR%s-%s' % (seqlist[j][1], CDRLoop, 1),100.0, 5, -100, max_skip = 0, quiet = 1)
                    record[(i+1)].append('LLL%s' % str(fitter[0])[:5])
                    
                else: 
                    record[(i+1)].append('Null')                  
                j -= 1

            i += 1
            record[(i)].append(str(0))            
        else:
            while (j >= i):
                record[(i+1)].append('Null')
                j -= 1
            i += 1
    record = process(record)
    
    print("There were %s poor fits" % (poorfit))

    cmd.remove('all')
    cmd.delete('all')
    saver(record, '%s' % loop[CDRLoop])

saver(IDloop, '%s ID loops' % loop[CDRLoop])
saver(NonConverge, 'Nonconverging loops')
saver(ch_id, 'Chain_ID')

cmd.orient()
print('Finished!')
