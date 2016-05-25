#! /usr/bin/env python2

#
# Copyright (C) Telecom ParisTech
# 
# This file must be used under the terms of the CeCILL. This source
# file is licensed as described in the file COPYING, which you should
# have received as part of this distribution. The terms are also
# available at:
# http://www.cecill.info/licences/Licence_CeCILL_V1.1-US.txt
#

import sys
import argparse

import des
import km
import operator
import numpy as np


import math




def equal_length(slow, fast):
    if len(slow) < len(fast):
        pass     
def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

def pearson_def(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    return diffprod / math.sqrt(xdiff2 * ydiff2)
    
def test():


    start = 0x041041041041
    step  = 0x041
    end   = 0x100000000000
    for i in range(0, end , start ):
        print bin(i)[2:].zfill(48)        
#print bin(i)[2:8].zfill(6), bin(i)[9:15].zfill(6), bin(i)[16:22].zfill(6), bin(i)[23:23+6].zfill(6), bin(i)[30:36].zfill(6), bin(i)[37:43].zfill(6), bin(i)[44:50].zfill(6), bin(i)[51:57].zfill(6)

def avg2(n1,n2):
    l = len(n1) + len(n2)
    res = float(sum(n1) + sum(n2))

    return res/l
def avg (n):
    l = float(len(n))
    return sum(n)/l
def index_max (n):
    return max(xrange(len(n)), key=values.__getitem__)

def new_attem():

    if not des.check ():
        sys.exit ("DES functional test failed")

    # *************************************
    # * Check arguments and read datafile *
    # *************************************
    argparser = argparse.ArgumentParser(description="Apply P. Kocher's TA algorithm")
    argparser.add_argument("datafile", metavar='file', 
                        help='name of the data file (generated with ta_acquisition)')
    argparser.add_argument("n", metavar='n', type=int,
                        help='number of experiments to use')
    args = argparser.parse_args()

    if args.n < 1:                                      # If invalid number of experiments.
        sys.exit ("Invalid number of experiments: %d (shall be greater than 1)" % args.n)

#    print args.n
    # Read encryption times and ciphertexts. n is the number of experiments to use.
    read_datafile (args.datafile, args.n)
    """ # *****************************************************************************
    # * Compute the Hamming weight of output of first (leftmost) SBox during last *
    # * round, under the assumption that the last round key is all zeros.         *
    # *****************************************************************************
    rk = 0x000000000000
    # Undoes the final permutation on cipher text of n-th experiment.
    r16l16 = des.ip (ct[args.n - 1])
    # Extract right half (strange naming as in the DES standard).
    l16 = des.right_half (r16l16)
    # Compute output of SBoxes during last round of first experiment, assuming
    # the last round key is all zeros.
    sbo = des.sboxes (des.e (l16) ^ rk)  # R15 = L16, K16 = rk
    # Compute and print Hamming weight of output of first SBox (mask the others).
    print >> sys.stderr, "Hamming weight: %d" % hamming_weight (sbo & 0xf0000000)

    # ************************************
    # * Compute and print average timing *
    # ************************************
    print >> sys.stderr, "Average timing: %f" % (sum (t) / args.n)

    # ************************
    # * Print last round key *
    # ************************
    print >> sys.stderr, "Last round key (hex):"
    print ("0x%012X" % rk)
    print sbo

    """ 
    slow = []
    fast = []
    sbo  = []
    for i in range(0,8):
        slow.append([]);
        fast.append([])
        for j in range(0,64):
            fast[i].append([])
            slow[i].append([])
            sbo.append([]) 
    mask = 0xf0000000
    subkey = 0x0
    for key in range(64):
        for k in range(0,args.n):
            sbo[key].append( des.sboxes(des.e(des.right_half(des.ip(ct[k])))^subkey))
        subkey  += 0x041041041041
    print args.n 
    for key in range(64):
        for k in range(args.n):
            mask = 0xf0000000
            for sbox in range(8):
                H = hamming_weight(sbo[key][k]&mask)
                if H==0 or H==1:
                    fast[sbox][key].append(t[k])
                    
                if H==4 or H==3:
                    slow[sbox][key].append(t[k])

                mask = mask >> 4            
    print '\n'
    losning = []
    differance = np.zeros(shape=(8,64))
    

    #for elem in slow:
    #    print slow


    print slow[0][0]
    for i in range(8):
        for j in range(64):
            differance[i][j] = avg(slow[i][j]) - avg(fast[i][j])
            print differance[i][j] 
        md = max(differance[i])
        print md
    for elem in differance:
        print elem




    """
    for j in range(64): 
        pass


    losning =[]
    for sbox in range(8):
        print 'testing for box number: '  + str(sbox)
        for key in range(64):
            pass        

    """
    #x = slow.sum(axis=1)
    #y = fast.sum(axis=1)
    #print min(x)
    #print 'x= ', x
    

    #for j in range(x.size):
    #    diff[0][j] = x[j]-y[j]

    #print diff
#    for sbox in range(8):
#        for key in range(64):
#            print slow[sbox,key]
#            diff[sbox][key] = slow.sum(axis=sbox) 
#        print '\n'    
        #diff[sbox][key] = sum(sum(slow[sbox,key] ))/ float(len(slow[sbox,key])) - sum(sum(fast[sbox,key]))/len(fast[sbox,key]) 
        
    #print diff

def main ():
    # ************************************************************************
    # * Before doing anything else, check the correctness of the DES library *
    # ************************************************************************
    if not des.check ():
        sys.exit ("DES functional test failed")

    # *************************************
    # * Check arguments and read datafile *
    # *************************************
    argparser = argparse.ArgumentParser(description="Apply P. Kocher's TA algorithm")
    argparser.add_argument("datafile", metavar='file', 
                        help='name of the data file (generated with ta_acquisition)')
    argparser.add_argument("n", metavar='n', type=int,
                        help='number of experiments to use')
    args = argparser.parse_args()

    if args.n < 1:                                      # If invalid number of experiments.
        sys.exit ("Invalid number of experiments: %d (shall be greater than 1)" % args.n)

#    print args.n
    # Read encryption times and ciphertexts. n is the number of experiments to use.
    read_datafile (args.datafile, args.n)
    """ # *****************************************************************************
    # * Compute the Hamming weight of output of first (leftmost) SBox during last *
    # * round, under the assumption that the last round key is all zeros.         *
    # *****************************************************************************
    rk = 0x000000000000
    # Undoes the final permutation on cipher text of n-th experiment.
    r16l16 = des.ip (ct[args.n - 1])
    # Extract right half (strange naming as in the DES standard).
    l16 = des.right_half (r16l16)
    # Compute output of SBoxes during last round of first experiment, assuming
    # the last round key is all zeros.
    sbo = des.sboxes (des.e (l16) ^ rk)  # R15 = L16, K16 = rk
    # Compute and print Hamming weight of output of first SBox (mask the others).
    print >> sys.stderr, "Hamming weight: %d" % hamming_weight (sbo & 0xf0000000)

    # ************************************
    # * Compute and print average timing *
    # ************************************
    print >> sys.stderr, "Average timing: %f" % (sum (t) / args.n)

    # ************************
    # * Print last round key *
    # ************************
    print >> sys.stderr, "Last round key (hex):"
    print ("0x%012X" % rk)
    print sbo

    """ 
    
    diff = {0:{},1:{},2:{},3:{},4:{},5:{},6:{},7:{}}

    slow = []
    fast = []

    for i in range(0,8):
        slow.append([]); 
        fast.append([])
        for j in range(0,64):
            fast[i].append([])
            slow[i].append([])
    for elem in slow:
        print elem

    print len(slow[0])
    start = 0x000000000000
    end   = 0xfc0000000000
    step = 0x040000000000
    counter=[]
    la=0
    mask = 0xf0000000
    cnt = 0
    subkey = 0x0 


        
    for key in range(64):
    #for sk in range(start, end+step,step):
        for j in range(args.n): 
            l16 = des.right_half(des.ip(ct[j]))
            R48  = des.e(l16)
            sbo = des.sboxes(R48^subkey)
            mask = 0xf0000000
            for sbox in range(0,8):
                H   = hamming_weight(sbo&mask)
                if H ==0 or H==1:
                    fast[sbox][key].append(t[j])
                    #print 'added fast'
                elif H ==4 or H == 3:
                    slow[sbox][key].append(t[j])
                mask = mask >> 4
        subkey  += 0x041041041041
    nokkel = []


    for elem in slow: 
        print elem
    print len(slow)
    print len(slow[0])
    print len(slow[0])
    for sbox in range(8):
        print 'testing for bxnr ' + str(sbox)

        for key in range(64):
            diff[sbox][key] = (sum (slow[sbox][key]) / len(slow[sbox][key]) - sum (fast[sbox][key]) / len(fast[sbox][key]))
        maxdiff = max(diff[sbox].values())
        for subk in diff[sbox].keys():
            if diff[sbox][subk] == maxdiff:
                nokkel.append(subk)
    for n in range(8):
        print hex(nokkel[n])





        
        subkey += 0x041041041041
        #f = open('workfile', 'w')
        #for i in range(len(score)):
        #    print score[i], counter[i]
        #    f.write(str(score[i]))
        #    f.write('\t')
        #    f.write(str(counter[i]))
        #    f.write('\n')
            
        #max_idx ,max_val = max(enumerate(score), key=operator.itemgetter(1)) #find max index
        #print max_idx, max_val, hex(max_idx)
        #print '\n\n'
        #slow = []
        #fast = []
        #score=[]
        #counter = []
        #la   = 0
        #end  = 0xfc0000000000 >> sbx
        #step = 0x040000000000 >> sbx
        #mask = mask >> 4
 
# Open datafile <name> and store its content in global variables
# <ct> and <t>.
def read_datafile (name, n):
    global ct, t

    if not isinstance (n, int) or n < 0:
        raise ValueError('Invalid maximum number of traces: ' + str(n))

    try:
        f = open (str(name), 'rb')
    except IOError:
        raise ValueError("cannot open file " + name)
    else:
        try:
            ct = []
            t = []
            for _ in xrange (n):
                a, b = f.readline ().split ()
                ct.append (int(a, 16))
                t.append (float(b))
				
        except (EnvironmentError, ValueError):
            raise ValueError("cannot read cipher text and/or timing measurement")
        finally:
            f.close ()

# ** Returns the Hamming weight of a 64 bits word.
# * Note: the input's width can be anything between 0 and 64, as long as the
# * unused bits are all zeroes.
# See: http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
def hamming_weight (v):
    v = v - ((v>>1) & 0x5555555555555555)
    v = (v & 0x3333333333333333) + ((v>>2) & 0x3333333333333333)
    return (((v + (v>>4) & 0xF0F0F0F0F0F0F0F) * 0x101010101010101) >> 56) & 0xFF

if __name__ == "__main__":
    #main ()
    new_attem()
    #test()
    print avg([1,2,3,4,5,6])
