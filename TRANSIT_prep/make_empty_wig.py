## Defne Surujon
## April 4, 2018

## Given a genbank file for a bacterial genome, make an empty
## .wig file, which contains the coordinates of TA sites that
## are unique (ie. cannot be double mapped during TnSeq data
## processing). The output .wig file is a template with no
## insertion data. 

from Bio import SeqIO
from Bio import Seq
from optparse import OptionParser
import regex
import re


options = OptionParser(usage='%prog -g [genome] -o [output]',
                       description="Specify input genome and output file")

options.add_option("-g","--genome",dest="genomefile",
                   help="genome file (.gbk)")

options.add_option("-o","--outfile",dest="outputfile",
                   help="output file (.wig)")



def get_ta(strain):
    seq = strain.seq
    seqstr = str(seq).upper()
    l = len(seqstr)
    ta_ALL = [i.start()+1 for i in re.finditer("TA",seqstr)]
    ta_ALL.sort()
    downstream = [seqstr[i+1:i+17] for i in ta_ALL if i<l-17]
    downstream = downstream + [seqstr[i+1:]+seqstr[:i+17-l] for i in ta_ALL if i>=l-17]
    upstream = [seqstr[i-17:i-1] for i in ta_ALL if i>17]
    upstream = [seqstr[l-(17-i):]+seqstr[:i-1] for i in ta_ALL if i<=17] + upstream
    pos = {i:0 for i in ta_ALL}
    return (pos, downstream,upstream, ta_ALL)


def get_repeats_RE(upstream, downstream, ta_sites):
    repeat_sites= []
    for i in range(0,len(ta_sites)):
        queryF = regex.compile("(" +downstream[i]+ "){e<=1}")
        queryR = regex.compile("(" +str(Seq.Seq(upstream[i]).reverse_complement())+ "){e<=1}")
        
        if len(list(filter(queryF.match,downstream)))>1:
            #print(queryF)
            repeat_sites += [ta_sites[i]]
        elif len(list(filter(queryR.match,downstream)))>1:
            #print(queryR)
            repeat_sites += [ta_sites[i]]

    return(repeat_sites)

def main():
    opts, args = options.parse_args()
    genomefile = opts.genomefile
    outfile = opts.outputfile
    # get all TA sites
    mystrain=SeqIO.read(genomefile,"genbank")
    ta_table, downstream, upstream, ta_sites = get_ta(mystrain)
    
    print("Total number of TA sites: ", len(ta_table))
    repsRE = get_repeats_RE(upstream, downstream, ta_sites)
    print("Number of duplicate TA sites: ", len(set(repsRE)))
    
    talist = list(set(ta_sites) - set(repsRE))
    talist.sort()
    g = open(outfile,'w')
    g.write('variableStep chrom='+mystrain.name+'\n')
    for tasite in talist:
        linestr = [str(tasite), "0"]
        g.write("\t".join(linestr)+"\n")
    g.close()



if __name__ == '__main__':
    main()
