# DS 6/4/18
# This script will take in a fasta file and add unique ID's to
# the fasta header.
# Going from
# ">Myheader"
# To
# ">[NN]Myheader"
# where the number NN is inserted right before the existing header
# for each sequence in the series.

from optparse import OptionParser


options = OptionParser(usage='%prog -i [input] -o [output]',
                       description="Specify input fasta and output file")

options.add_option("-i","--infile",dest="inputfile",
                   help="input file (.fasta)")

options.add_option("-o","--outfile",dest="outputfile",
                   help="output file (.fasta)")


def add_IDs(infile, outfile):
    with open(infile) as f:
        lines = f.readlines()
        c=1
        with open(outfile,'w') as g:
            for line in lines:
                if line[0] == ">":
                    modline = ">"+"["+str(c)+"]"+line[1:]
                    c+=1
                else:
                    modline = line
                g.write(modline)
                
                
            


def main():
    opts, args = options.parse_args()
    infile = opts.inputfile
    outfile = opts.outputfile

    add_IDs(infile, outfile)


if __name__ == '__main__':
    main()
