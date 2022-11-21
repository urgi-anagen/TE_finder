#!/usr/bin/env python

import os
import sys
import numpy
import shutil

input_file = sys.argv[1]
#output_file = sys.argv[2]


align_file = open(input_file)
frag_len=[]
for l in align_file:
    tok=l.split('\t')
    frag_len.append(abs(int(tok[1])-int(tok[2]))+1)

amin=numpy.amin(frag_len)
amax=numpy.amax(frag_len)
q25=numpy.quantile(frag_len,0.25)
q50=numpy.quantile(frag_len,0.5)
q75=numpy.quantile(frag_len,0.75)
print(input_file,'-> read ',len(frag_len),'fragments : min={} Q25={} Q50={} Q75={} max={}'.format(amin,q25,q50,q75,amax))