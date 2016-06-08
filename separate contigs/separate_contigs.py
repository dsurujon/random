##Defne Surujon
##8 June 2016

##Read a fasta file, and make new fasta files with each individual sequence

import os
from optparse import OptionParser

options = OptionParser(usage='%prog input output ',
                       description="Specify input fasta file and output directory")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.fasta)")
options.add_option("-o","--outdir",dest="outputdir",
                   help="output directory")

#read fasta file, return dictionary
def read_sequences(filename):
    f=open(filename)
    lines=f.readlines()
    f.close()
    titles={}
    for i in lines:
        if i[0]==">":
            title=i[:-1]
            titles[title]=""
        else:
            titles[title]+=(i.upper()).replace("\n","")
    return titles

#write a new fasta file for each contig
def write_contigs(c,filename):
    os.makedirs(filename)
    for i in c:
        thisfile=filename+"/"+i[1:]+".fasta"
        f=open(thisfile,'w')
        f.write(i+"\n")
        f.write(c[i])
        f.close()


def main():
    opts, args = options.parse_args()
    allcontigs=read_sequences(opts.inputfile)
    write_contigs(allcontigs,opts.outputdir)

if __name__ == '__main__':
    main()
