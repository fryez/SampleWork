#Code by Zach Frye with general python functions

import itertools
import operator
import time
import os
import sys
import datetime


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
