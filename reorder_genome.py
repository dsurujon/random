## Defne Surujon
## January 15 2019

## Find dnaA in a genbank, make a fasta file
## where the first nucleotide is the start of gene dnaA
## in the forward orientation.



from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser
import re
import os


options = OptionParser(usage='%prog -i [inputdir] ',
                       description="Specify input (containing .gbk files)")

options.add_option("-i","--indir",dest="inputdir",
                   help="input directory of genbank files")

def reorder_gbk(infile, outfile):
    mystrain=SeqIO.read(infile,"genbank")
    strainseq = mystrain.seq
    strainID = mystrain.id
    for feature in mystrain.features:
        try:
            if feature.type == "CDS" and feature.qualifiers['gene'][0]=='dnaA':
                startloc = feature.location.nofuzzy_start
                strand = feature.location.strand
        except KeyError: pass
    if strand == 1:
        newseq = strainseq[startloc:] + strainseq[:startloc]
    else:
        newseq = strainseq[:startloc].reverse_complement() + strainseq[startloc:].reverse_complement()
    newrecord = SeqRecord(newseq, id = strainID)
    with open(outfile, 'w') as f:
        SeqIO.write(newrecord, f, 'fasta')
    
    

def main():
    opts, args = options.parse_args()
    inputdir = opts.inputdir

    gbkfiles = os.listdir(inputdir)

    for infile in gbkfiles:
        print("processing ", infile)
        
        fullinputfile = os.path.join(inputdir, infile)
        outfile = '.'.join(infile.split('.')[:-1])+".fasta"
        fulloutputfile = os.path.join(inputdir, outfile)

        reorder_gbk(fullinputfile, fulloutputfile)
        



if __name__ == '__main__':
    main()
