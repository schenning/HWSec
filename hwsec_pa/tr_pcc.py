#!/usr/bin/env python2

"""The tr_pcc library, a software library dedicated to the computation of
Pearson Correlation Coefficients (PCC).

:Date: 2009-07-29
:Authors:
    - Renaud Pacalet, renaud.pacalet@telecom-paristech.fr

Defines a data structure and a set of functions used to compute and manage
Pearson Correlation Coefficients (PCC) between a floating point (float),
one-dimension, vector random variable and a set of floating point scalar
random variables.  These coefficients are a statistical tool to evaluate the
correlation between random variables. In this library, the vectors are
considered as a collection of scalar floating point random variables and each
component of the vectors is considered as independant. In the following,
vectors are denoted Z[.] and the i-th component of vector Z[.] is denoted
Z[i]. The formula of a PCC between the vector random variable X[.] and the
scalar random variable Y is::

PCC(X[.], Y)[.] = {E(X[.]*Y) - E(X[.])*E(Y)[.]} / {s(X[.]) * s(Y)[.]}

where E(Z) is the expectation (average value) of random variable Z and s(Z)
is its unbiased standard deviation. E(Z[.])[.] is the expectation vector of
the vector random variable Z[.]: each component of E(Z[.])[.] is the
expectation of the corresponding component of Z[.]: E(Z[.])[i] = E(Z[i]).
Similarly, s(Z[.])[.] is the unbiased standard deviation vector of the vector
random variable Z[.]: s(Z[.])[i] = s(Z[i]). The product between a vector and
a scalar is defined as: (X[.]*Y)[i] = X[i]*Y. The PCC between a vector random
variable X and a scalar random variable Y is a vector defined by::

PCC(X[.], Y)[i] = PCC(X[i], Y)

The tr_pcc library can be used to compute a set of PCCs between one vector
random variable (denoted X[.]), common to all PCCs, and a set of scalar
random variables (denoted Y0, Y1, ..., Yn-1). To compute such a set of PCCs
one must first instantiate a PCC context (`pccContext`), indicating the
number ny of Y random variables. Then, realizations of the random variables
must be accumulated into the context: first a realization of the X[.]
variable (`pccContext.insert_x()`), followed by realizations of each of the ny Y
variables (`pccContext.insert_y()`). Once a sufficient number of realizations are
accumulated in the PCC context, a call to `pccContext.consolidate()` computes the
ny PCCs. Calls to `pccContext.get_pcc()` return the different vector PCCs. Note:
more realizations can be accumulated after a call to `pccContext.consolidate()`,
another call to `pccContext.consolidate()` will take them into account altogether
with the previous ones and compute the new PCC values. Example of use with
10-components vectors and ny=4 Y scalar variables, in which `get_next_x` and
`get_next_y` are two functions returning realizations of the random variables::

import tr_pcc
import sys
...
ctx = tr_pcc.pccContext (10, 4)     # 10-samples traces, 4 PCCs
for i in range (nexp):              # For experiments
    x = get_next_x ()
    ctx.insert_x (x)
    for j in range (4):             # For PCCs
        y = get_next_y (j)
        ctx.insert_y (j, y)

ctx.consolidate ()

for j in range (4):                 # For PCCs
    pcc = ctx.get_pcc (j)
    sys.stdout.write ("PCC(X[.], Y%d)[.] =" % j)
    for i in range(10):             # For samples
        sys.stdout.write (" %lf" % pcc[i])
    sys.stdout.write ("\n")

ctx = tr_pcc.pccContext (100, 12);  # 100-samples traces, 12 PCCs
...

Attention
=========

It is an error to break the realization insertion scheme: if you initialized
your PCC context for ny Y variables, first insert a realization of X[.],
followed by one and only one realization of each of the ny Y variables. Then,
insert a new realization of X[.] and ny new realizations of the ny Y
variables, and so on.  Consolidate only after inserting the realization of
the last Y variable. Do not consolidate when in an intermediate state.

"""

import math

class pccContext:
    """The data structure used to compute and manage a set of Pearson correlation
    coefficients.

    Attributes:
        ny (int): Number of Y random variables
        nr (int): Current number of realizations of the random variables
        l (int): Size of the vector random variable (number of components)
        pcc (List[List[float]]): List of the PCCs

    """

    def __init__ (self, l, ny):
        """Initializes a PCC context.

        Args:
            l (int): Length of the X[.] vector random variable (number of components)
            ny (int): Number of Y random variables to manage

        """
        
        if not isinstance (l, int) or not isinstance (ny, int):
            raise TypeError ('parameters l and ny should be numbers')
            
        if l < 1 or ny < 1:
            raise ValueError('Invalid parameters value: l=%d, ny=%d' % (l, ny))

        self.ny = ny
        self.nr = 0
        self.l = l
        self.__x = [0.0] * l
        self.__x2 = [0.0] * l
        self.__y = [0.0] * ny
        self.__y2 = [0.0] * ny
        self.__xy = [[0.0] * l] * ny
        self.pcc = [[0.0] * l] * ny
        self.__state = 0
        self.__flags = [0] * ny
        self.__rx = None

    def insert_x (self, x):
        """Insert a realization of X[.] in a PCC context

        Args:
            x (List[float]): The realization of X[.]

        """

        if len(x) < self.l:
            raise ValueError ('size of realization of X[.] should be %d but is only %d' % (self.l, len(x)))

        if any ([d != self.__state for d in self.__flags]):
            raise ValueError ('missing realizations %s' % ', '.join(['Y%d' % d for d in range(ny) if self.__flags[d] != self.__state]))

        self.__rx = x[:self.l]
        self.__x = [a+b for a,b in zip (self.__x, self.__rx)]
        self.__x2 = [a+b*b for a,b in zip (self.__x2, self.__rx)]
        #print "x : " + str(self.__x)
        #print "x2 : " + str(self.__x2)
        #print "rx : " + str(self.__rx)
        self.__state = 1 - self.__state
        self.nr += 1

    def insert_y (self, ny, y):
        """Insert a new Y realization in a PCC context

        Args:
            ny (int): Index of the Y random variable (0 to self.ny - 1)
            y (float): Realization of the Y random variable

        """

        
        if not isinstance (ny, int) or not (isinstance (y, float) or isinstance (y, int)):
            raise TypeError('parameters ny and y should be numbers')

        if ny < 0 or ny >= self.ny:
            raise ValueError ('Invalid Y index: %d' % ny)
        if self.__rx is None:
            raise ValueError ('a realization of X[.] should be first inserted')
        if self.__flags[ny] == self.__state:
            raise ValueError ('Y realization #%d inserted twice' % ny)

        self.__y[ny] += y
        self.__y2[ny] += y*y
        self.__xy[ny] = [a+b*y for a,b in zip (self.__xy[ny], self.__rx)]
        self.__flags[ny] = self.__state
        #print " ===== insert_y (%d, %f) ===== " % (ny, y)
        #print "y : " + str(self.__y)
        #print "y2 : " + str(self.__y2)
        #print "xy : " + str(self.__xy)

    def consolidate (self):
        """Consolidates a set of PCCs (computes all the PCCs from the already
        inserted realizations)

        """

        if any ([d != self.__state for d in self.__flags]):
            raise ValueError ('missing realizations %s' % ', '.join(['Y%d' % d for d in range(ny) if self.__flags[d] != self.__state]))

        if self.nr < 2:
            raise ValueError ('not enough realizations (%d, min 2)' % self.nr)

        tmp_v = [math.sqrt (self.nr * a - b*b) for a,b in zip (self.__x2, self.__x)]
        self.pcc = [[(self.nr*xy - x*y)/tmp/math.sqrt (self.nr*y2 - y*y) for xy,x,tmp in zip (xy_v, self.__x, tmp_v)] for y,y2,xy_v in zip (self.__y, self.__y2, self.__xy)]

    def get_pcc (self, ny):
        """Returns a pointer to the PCC vector number ny

        Args:
            ny (int): Index of the PCC to get (0 to self.ny - 1)

        Returns:
            a pointer to the PCC vector number ny

        """

        if not isinstance (ny, int):
            raise TypeError ('parameter ny should be a number')
        if ny < 0 or ny >= self.ny:
            raise ValueError ('Invalid Y index: %d' % ny)

        return self.pcc[ny]
