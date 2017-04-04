##Defne Surujon
##1 April 2017
##Split a fasta file conaining protein sequences into
##groups based on the protein name in the fasta header
##annotation

import os
from optparse import OptionParser
import re

options = OptionParser(usage='%prog -i input -o output ',
                       description="Specify input file and output directory")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.fasta)")
options.add_option("-o","--outdir",dest="outputdir",
                   help="output directory")

#find string between two specified strings
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""
    
#read fasta file, return dictionary
def read_sequences(filename):
    f=open(filename)
    titles={}
    clusters={}
    first=True
    for i in f.readlines():
        if i[0]==">":
            if first==False and protein_name in clusters:
                clusters[protein_name]+=title+"\n"+titles[title]+"\n"
            elif first==False:
                clusters[protein_name]=title+"\n"+titles[title]+"\n"
            first=False
            title=i.strip()
            protein_name=find_between(title, "[protein=","]")
            titles[title]=""
        else:
            titles[title]+=(i.upper()).strip()
    f.close()
    return titles,clusters

def make_write_friendly(fname):
    return re.sub('[/\:|*?<>]','',fname)

def make_clusters(inputfasta, outputdirectory):
    myseqs,myclusts=read_sequences(inputfasta)
    for cluster in myclusts:
        clusterfile=make_write_friendly(cluster)
        outfilename=os.path.join(outputdirectory,clusterfile+".fasta")
        cfile=open(outfilename,"w")
        cfile.write(myclusts[cluster])
        cfile.close()

def main():
    opts, args = options.parse_args()
    allproteins=opts.inputfile
    clustersdir=opts.outputdir
    if not os.path.exists(clustersdir):
        os.makedirs(clustersdir)
    make_clusters(allproteins,clustersdir)

if __name__ == '__main__':
    main()
