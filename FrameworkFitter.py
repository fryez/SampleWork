@@ -0,0 +1,283 @@
#This program is to find an appropriate framework for peptide binder
#validation work.  PDB files are compared to the model based on both
#framework and later antigen-binding similarity
#Will need to modify the alignment parameters to find a suitable match

#from pymol import cmd, CmdException, stored;
from operator import itemgetter;
import os;
import string;
import datetime;
import re;
import types;
import random;
import csv;

FILETARGET = 'FILEPATH';

#Open the list of PDB files to fetch from the server to analyze
def opener(filename):
    flist = []
    file = open( '%s.csv' % (filename))
    line = file.readline()
    count = 0
    while (line != ''):
        flist.append([])
        parts = line.strip().split(',')
        for part in parts:
            flist[count].append(part)
        count += 1
        line = file.readline()
    file.close()

    return flist

#Fetches the PDB files and prepares them for alignment
def fetch(j):
    cmd.fetch('%s' % (j))
    cmd.hide('everything', '%s' % (j))
    cmd.select('nonbond', '%s and !(n. N xt. 20)' % (j))
    cmd.remove('nonbond')
    cmd.delete('nonbond')
    cmd.show('ribbon', '%s' % (j))
    cmd.delete('side')
    cmd.util.cbc('%s' % (j))

def printer(m):
    for x in m:
        print(x)
    print('')

#Generates separate objects for each peptide chain in a PDB file
def split_chains(selection='(all)', prefix=None):
    count = 0
    models = cmd.get_object_list('(' + selection + ')')
    for model in models:
        for chain in cmd.get_chains('(%s) and model %s' % (selection, model)):
            if chain == '':
                chain = "''"
            count += 1
            if not prefix:
                name = '%s%s' % (model, count)
            else:
                name = '%s%04d' % (prefix, count)
            cmd.create(name, '(%s) and model %s and chain %s' % (selection, model, chain))
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
            print (' Error: startsele not in selection')
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
        print (' Renumber: range (%d to %d)' % tuple(minmax))

cmd.extend('renumber', renumber)

#Isolated the sequences between the C and the W designating
#CDR1 in both alpha and beta chains
def sequencer(i, CDRLoop):
    cmd.select("__h", "br. " + '%sCDR%s' % (seqlist[i][1], CDRLoop) + " and not het")
    AAs = ''
    aaDict = { 'aaList' : [] }
    cmd.iterate("(name ca) and __h","aaList.append((resi,resn))",space=aaDict)
    IDs = map( lambda x: str(x[0]), aaDict['aaList'] )
    AAs = ''.join(map( lambda x: one_letter[x[1]], aaDict['aaList']))
    while (AAs[:1] != 'C'):
        AAs = AAs[1:]
    AAs = AAs[::-1]
    while (AAs[:1] != 'W'):
        AAs = AAs[1:]
    AAs = AAs[::-1]

    newCDR.append([])
    newCDR[i].append(seqlist[i][0])
    newCDR[i].append(AAs)
    return

#Returns all the items in a column of a matrix
def column(matrix, i):
    return [row[i] for row in matrix]

#Saves the alignment tables into individual loops in the new folder
def saver(rec, label):
    outFile = open('%s.csv' % (label), 'w');
    for i in rec:
        temp = ''
        for j in i:
            temp = temp + ',' + str(j)
        temp += ',\n'    
        outFile.write('%s, ' % (temp))
    outFile.close();
    os.chdir(FILETARGET);
    return

#Used to save an image of the candidate for a quick visual check later
def image(name, chain):
    cmd.select('X', 'bc. (%s w. 0.5 of (target and chain %s))' % (name, chain))
    cmd.save('%s-%s.pdb' % (i[0], chain), 'X')
    cmd.select('X', '!X and (%s and chain %s)' % (name, chain))
    cmd.remove('X')
    cmd.delete('X')
    #cmd.color('green', '%s' % (name))
    cmd.color('red', 'target and chain H')
    cmd.color('blue', 'target and chain L')
    cmd.show('ribbon', '%s' % (name))
    #cmd.show('sticks', 'target')
    cmd.orient('I')
    cmd.zoom('target and chain %s' % (chain))
    cmd.ray(1000,1000)
    cmd.png('%s-%s' % (name, chain), dpi = 1000)
    return

#Used to compare features of antibodies from the PDB to the model loop framework and CDR loops.
def align(i, chain, ResC):
    #The target in this function is the model structure of complex 14 and the selected chains.
    RMS = cmd.align('%s' % (i[0]), 'target and chain %s' % (chain),3.0,10,-10.0, -0.5, 50, 'Alignment', quiet = 1)
    RMSD = [RMS[0], RMS[1]]

    #If positive hit based on overall framework alignment and position
    if RMSD[0] < 1.0 and RMSD[1] >= (ResC / 3.0):       
        if chain == 'H+L':
            Inum = 36
        else:
            Inum = 18
            
        #Determine if the interface is positionally similar to the target interface
        d = 2.0
        X_I = cmd.select('target_I', 'n. CA and (%s w. %s of (I and chain %s))' % (i[0], d, chain))
        while X_I > Inum and d > 0:
            X_I = cmd.select('target_I', 'n. CA and (%s w. %s of (I and chain %s))' % (i[0], d, chain))
            d -= 0.05
        iRMSD = [10,10]

        #Only does the comparision if the interface is fairly similar
        if X_I == Inum:
            iRMS = cmd.align('target_I', 'I and chain %s' % (chain),4.0,10,-10.0, -0.5, 50, 'I_alignment', quiet = 1)
            iRMSD = [iRMS[0], iRMS[1]]
        if iRMSD[0] < 1.0:
            image(i[0], chain)           
            cmd.save('%s-%s.pse' % (i[0], chain))
            if (chain == 'H+L'):
                cmd.remove('%s' % (i[0]))
                cmd.delete('%s' % (i[0]))

            i.append('Full')
            for j in RMSD:
                i.append(j)
            for j in iRMSD:
                i.append(j)            
        else:
            if (chain == 'H+L'):
                cmd.remove('%s' % (i[0]))
                cmd.delete('%s' % (i[0]))
            i.append('Partial')
            for j in RMSD:
                i.append(j)
            for j in iRMSD:
                i.append(j)       
    else:
        i.append('No')
        for j in RMSD:
            i.append(j)
        i.append(10)
        i.append(10)
        if (chain == 'H+L'):
            cmd.delete('%s' % (i[0]))
            os.remove('%s.pdb' % (i[0]))
    cmd.delete('Alignment')
    cmd.delete('target_I')
    return i
        
#To clear out all the information from the previous program Run
cmd.delete('all')
cmd.reinitialize()
#Make the target filename fo the record output

os.chdir(FILETARGET);
outDir = datetime.datetime.now().strftime('FrameHits-%Y-%m-%dT%H-%M-%S');
os.mkdir(outDir);  #Only to be used if there will be a folder from the processes
cutoff = 1.0

#Load the candidate loop and interface to compare the PDB structures against
cmd.load('Loop.pdb')
cmd.hide('everything', 'Loop')
cmd.show('ribbon', 'Loop')
cmd.util.cbc('Loop')
cmd.load('Interface.pdb')
cmd.hide('everything', 'Interface')
cmd.show('ribbon', 'Interface')
cmd.select('I', 'n. CA in Interface')

#cmd.select('complex', 'n. CA in (chain A and Chain B)')
ResCountH = cmd.select('target', 'chain H and n. CA in Loop')
ResCountL = cmd.select('target', 'chain L and n. CA in Loop')
ResCount = cmd.select('target', 'n. CA in Loop')

flist = opener('pdblist')
os.chdir('%s/%s' %(FILETARGET, outDir))

positive = []
NotFound = []
for i in flist:
    Fit = 0
    try:
        try:
            fetch(i[0])
            align(i, 'H', ResCountH)
            if i[2] < 1.0 and i[3] >= (ResCountH / 3.0):
                Fit = 1
            align(i, 'L', ResCountL)
            if i[4] < 1.0 and i[5] >= (ResCountL / 3.0):
                Fit = 1
            align(i, 'H+L', ResCount)
            if i[6] < 1.0 and i[7] >= (ResCount / 3.0):
                Fit = 1
            if Fit == 1:
                positive.append(i[0])
        except Selector-Error:
            None
    except NameError:
        print('Selector-Error on %s' % (i[0]))

saver(flist, 'FrameHits')

print('Finished!')
