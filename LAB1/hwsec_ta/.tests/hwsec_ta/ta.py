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
import time
import threading




def avg (n):
    l = float(len(n))
    return sum(n)/l




def m_th(sbox, slow, fast,n,sbo, mask): 
    

    for key in range(64):
        for k in range(n):
            H = hamming_weight(sbo[key][k]&mask)
            if H==0 or H==1:
                fast[sbox][key].append(t[k])                    
            if H==4 or H==3:
                slow[sbox][key].append(t[k])


class myThread (threading.Thread):
    def __init__(self, sbox, slow,fast,n,sbo, mask):
        threading.Thread.__init__(self)
        self.sbox = sbox
        self.slow = slow
        self.fast = fast
        self.n    = n
        self.sbo  = sbo
        self.mask = mask



  
    def run(self):
        m_th(self.sbox, self.slow, self.fast, self.n, self.sbo, self.mask)
def main():
    
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
    differance = {0:{},1:{},2:{},3:{},4:{},5:{},6:{},7:{}}
    global slow
    global fast 
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
    mask = 0xf0000000
    t1 = myThread(0, slow,fast,args.n,sbo, mask)
    t2 = myThread(1, slow,fast,args.n,sbo, mask >>4)
    t3 = myThread(2, slow,fast,args.n,sbo, mask >>8)
    t4 = myThread(3, slow,fast,args.n,sbo, mask >>12)
    t5 = myThread(4, slow,fast,args.n,sbo, mask >>16)
    t6 = myThread(5, slow,fast,args.n,sbo, mask >>20)
    t7 = myThread(6, slow,fast,args.n,sbo, mask >>24)
    t8 = myThread(7, slow,fast,args.n,sbo, mask >>28)



    t1.start()
    t2.start()
    t3.start()
    t4.start()
    t5.start()
    t6.start()
    t7.start()
    t8.start()

    threads = []

    threads.append(t1)
    threads.append(t2)
    threads.append(t3)
    threads.append(t4)
    threads.append(t5)
    threads.append(t6)
    threads.append(t7)
    threads.append(t8)
    for t in threads: 
        t.join()
    res = ''
    losning = []
    for sbox in range(8):
        for key in range(64):
            differance[sbox][key] = (avg(slow[sbox][key]) - avg(fast[sbox][key]))
        md = max(differance[sbox].values())
        for subk in differance[sbox].keys():
            if differance[sbox][subk] == md:
                losning.append(subk)
    for i in range(len(losning)):
        res = res + bin(losning[i])[2:].zfill(6)
    print hex(int(res,2))









 
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
    main()
