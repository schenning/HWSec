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


def avg2(n1,n2):
    l = len(n1) + len(n2)
    res = float(sum(n1) + sum(n2))

    return res/l
def avg (n):
    l = float(len(n))
    return sum(n)/l
def index_max (n):
    return max(xrange(len(n)), key=values.__getitem__)


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
    

    b='110101'
    slow = []
    fast = []
    score = []
    start = 0x000000000000
    end   = 0xfc0000000000
    step = 0x040000000000
    counter=[]
    la=0
    mask = 0xf0000000
    key = []
    dic = {}
    ki = 0
    #       0x0000000f;
    for sbx in range(6,43,6):
        ki=0
        for sk in range(start, end,step):
            for j in range(args.n):
       
                l16 = des.right_half(des.ip(ct[j]))
                R48  = des.e(l16)
                sbo = des.sboxes(R48^sk)
                #print bin(sbo&mask)
                H   = hamming_weight(sbo & mask)

                if H ==0:
                    fast.append(t[j])
                elif H ==4:
                    slow.append(t[j])	
            score.append(avg(slow) - avg(fast))
            dic[ki] = sum(slow)/float(len(slow)) - sum(fast)/float(len(fast))
            ki +=1
            counter.append(la)
            la+=1
        maxdiff = max(dic.values())
        for sla in dic.keys():
            if dic[sla] == maxdiff:
                print hex(sla)



      #  for i in range(len(score)):
      #      print score[i], counter[i]
      #  max_idx ,max_val = max(enumerate(score), key=operator.itemgetter(1)) #find max index
      #  print max_idx, max_val, hex(max_idx)

        print '\n\n'
        slow = []
        dic={}
        fast = []
        score=[]
        counter = []
        la   = 0
        end  = 0xfc0000000000 >> sbx
        step = 0x040000000000 >> sbx
        mask = mask >> 4 
#        mask = mask>>4  
#        end = end >> 6
#        step = step >>5
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
    main ()
    print avg([1,2,3,4,5,6])
