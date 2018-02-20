## Defne Surujon
## Feb 20, 2018

# make a wig file from a genbank and map file
# map file is in the output of TnSeq

from Bio import SeqIO
from optparse import OptionParser
import re


options = OptionParser(usage='%prog input output',
                       description="Specify input directory and output file")

options.add_option("-g","--genome",dest="genomefile",
                   help="genome file (.gbk)")
options.add_option("-i","--infile",dest="inputfile",
                   help="input map file (.map)")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output file (.prot_table)")


def make_map_dict(mapfilename):
    mapdict = {}
    f = open(mapfilename)
    flines = f.readlines()
    f.close()
    for i in flines:
        linelist = i.split()
        thiscount = int(linelist[0])
        thisstrand = linelist[1]
        thispos = int(linelist[2])
        thislen = int(linelist[3])
        if thisstrand == "+":
            thispos += thislen - 2
        
        if thispos in mapdict:
            mapdict[thispos] += thiscount
        else:
            mapdict[thispos] = thiscount
    return(mapdict)

def get_ta(strain):
    seq = strain.seq
    seqstr = str(seq).upper()
    taF = [i.start()+1 for i in re.finditer("TA",seqstr)]
    taR = [i.start()+1 for i in re.finditer("AT",seqstr)]
    ta_ALL = taF+taR
    ta_ALL.sort()
    return {i:0 for i in ta_ALL}

def main():
    opts, args = options.parse_args()
    genomefile = opts.genomefile
    infile = opts.inputfile
    outfile = opts.outputfile
    # get all TA sites
    mystrain=SeqIO.read(genomefile,"genbank")
    ta_table = get_ta(mystrain)

    # get occupied TA sited
    mdict = make_map_dict(infile)

    # add occupied ta sites to the master list
    for tasite in mdict:
        try:
            ta_table[tasite] = mdict[tasite]
        except KeyError:
            print("WARNING: Site "+str(tasite)+"does not appear in the genome!!")
    

    g = open(outfile,'w')
    g.write('variableStep chrom='+mystrain.name+'\n')
    for tasite in ta_table:
        linestr = [str(tasite), str(ta_table[tasite])]
        g.write("\t".join(linestr)+"\n")
    g.close()



if __name__ == '__main__':
    main()
