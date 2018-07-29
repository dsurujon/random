## DS  28 Aug 2018
## Add KEGG pathway information and annotation to each gene
## on an RNAseq dataset

import pandas as pd
import os
from optparse import OptionParser

options = OptionParser(usage='%prog -i [input directory] -a [annot]',
                       description="Specify input directory that has all RNAseq data "
                       "files to be merged with and annotation file. By default, this "
                       "script uses the Gene column in the data file and the TIGR4.old "
                       "column in the annotation file")

options.add_option("-i","--input",dest="inputdir",
                   help="input directory")
options.add_option("-a","--annot",dest="annotfile",
                   help="annotation file")

def merge_single(RNAseqfilename, annot_DF):
    rnadata = pd.read_csv(RNAseqfilename)
    mergeddata = rnadata.merge(annot_DF, how = "outer", left_on = 'Gene', right_on = 'TIGR4.old')
    outfilename = os.path.splitext(RNAseqfilename)[0]+"_Annotated.csv"
    mergeddata.to_csv(outfilename)



def main():
    opts, args = options.parse_args()
    inputdir = opts.inputdir
    annotfile = opts.annotfile
    annot_DF = pd.read_csv(annotfile)
    datafiles = os.listdir(inputdir)

    for file in datafiles:
        fullfilename = os.path.join(inputdir,file)
        print(file)
        try:
            merge_single(fullfilename,annot_DF)
        except (RuntimeError, TypeError, NameError,ValueError, KeyError):
            print("File "+file+" could not be processed\n")

if __name__ == '__main__':
    main()
