#Defne Surujon  20 March 2017
#Use this to spike protein fasta sequences with mutations

import os
import random
from optparse import OptionParser

options = OptionParser(usage='%prog -i [inputfile] -o [outputfile] -m [mutation rate]',
                       description="Specify input fasta, output fasta and mutation rate (between 0 and 1)")

options.add_option("-i","--input",dest="inputfile",
                   help="Input .fasta file")
options.add_option("-o","--output",dest="outputfile",
                   help="Output .fasta file")
options.add_option("-m","--mutate",dest="mutationrate",
                   help="Mutation rate (between 0 and 1)")

aas=["A","V","M","S","T","Y","F","K","D","E",
     "Q","W","R","I","L","P","G","H","C","N"]

#change input amino acid into something else
#all amino acids have equal probability
def swap_single_aa(aa):
    newaa=aa
    while newaa==aa:
        newaa=aas[int(random.random()*20)]
    return newaa

def mutate_protein(seq,mutrate):
    for i in range(0,len(seq)):
        if random.random()<=mutrate:
            seq=seq[:i]+swap_single_aa(seq[i])+seq[i+1:]
        else: pass
    return seq

def process_fasta(filename,mutrate,outfilename):
    f=open(filename)
    g=open(outfilename,'w')
    for line in f:
        if line[0]==">":g.write(line)
        else:
            g.write(mutate_protein(line.strip().upper(),mutrate)+"\n")
    f.close()
    g.close()

#main
def main():
	opts, args = options.parse_args()
	inputfilename=opts.inputfile
	outputfilename=opts.outputfile
	mutationrate=float(opts.mutationrate)
	process_fasta(inputfilename,mutationrate,outputfilename)
##/main

if __name__ == '__main__':
    main()
