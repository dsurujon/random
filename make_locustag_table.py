## Defne Surujon
## August 27, 2019

## Make a correspondence table of old vs new locus tags for each gene
## in a gbk file


import pandas as pd
from Bio import SeqIO
from Bio import Seq
from optparse import OptionParser



options = OptionParser(usage='%prog -i [inputfile] -o [outputfile]',
                       description="Specify input file(.gbk) and output file (for .csv) ")

options.add_option("-i","--input",dest="inputfile",
                   help="input genbank file")
options.add_option("-o","--output",dest="outputfile",
                   help="output csv file")

def make_LT_table(inputfile, outputfile):
    mystrain = SeqIO.read(inputfile, 'genbank')
    oldLT_list = []
    newLT_list = []
    for feature in mystrain.features:
        if feature.type=="gene":
            newlocus = feature.qualifiers['locus_tag']
            try:
                oldlocus = feature.qualifiers['old_locus_tag']
            except KeyError:
                oldlocus = [""]
            oldLT_list.append(oldlocus[0])
            newLT_list.append(newlocus[0])
    LT_table = pd.DataFrame({'Old': oldLT_list, 'New':newLT_list})
    LT_table.to_csv(outputfile)

def main():
    opts, args = options.parse_args()
    inputfile = opts.inputfile
    outputfile = opts.outputfile

    make_LT_table(inputfile, outputfile)



if __name__ == '__main__':
    main()
