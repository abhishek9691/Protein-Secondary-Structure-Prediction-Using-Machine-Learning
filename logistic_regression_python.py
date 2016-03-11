#!/usr/bin/env python


import sys
import math
import operator


ITERATIONS = 100
LEARNING_RATE = 0.05


#
# Load the data, which must be lines of the form x1 x2 ... xn y
#
D = [[1.] + map(float, line.split()) for line in open(sys.argv[1])]


def p1(d, w):
    return 1.0 / (1.0 + math.exp(-sum(map(operator.mul, w, d[:-1]))))

    
n = len(D[0]) - 1
w = [0.0] * n

for iter in range(ITERATIONS):
    g = [0.0] * n
    for d in D:
        error = d[-1] - p1(d, w)
        for j in range(n):
            g[j] += error * d[j]
    for j in range(n):
        w[j] += LEARNING_RATE * g[j]


for d in D:
    print str(d) + ', p1 = ' + str(p1(d, w))

print '\nWeight vector: ' + str(w)
