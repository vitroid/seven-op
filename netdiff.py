#!/usr/bin/env python

#read two @NGPH and one coordinate
#usage : cat a.ngph a.ar3a | this b.ngph

import sys

b = sys.argv[1]

def loadNGPH(fh):
    hb = set()

    line = fh.readline()
    nmol = int(line)
    while True:
        line = fh.readline()
        x,y = map(int,line.split())
        if x < 0:
            break
        bond = (x,y)
        hb.add(bond)
    return hb


def loadAR3A(fh):
    coord = []
    line = fh.readline()
    nmol = int(line)
    for i in range(nmol):
        line = fh.readline()
        xyz = [float(i) for i in line.split()]
        coord.append(xyz)
    return coord



#given arg
fh = open(b,"r")
while True:
    line = fh.readline()
    if len(line) == 0:
        break
    columns = line.split()
    if columns[0] == "@NGPH":
        hbb = loadNGPH(fh)


#stdin

while True:
    line = sys.stdin.readline()
    if len(line) == 0:
        break
    columns = line.split()
    if columns[0] == "@NGPH":
        hba = loadNGPH(sys.stdin)
    elif columns[0] == "@AR3A":
        coord = loadAR3A(sys.stdin)

diff = hba^hbb
#print len(diff)
print diff
