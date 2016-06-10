##Defne Surujon
##10 June 2016

##Read a fasta file, and reorder the sequences based on a preferences file

from optparse import OptionParser

options = OptionParser(usage='%prog input preferences',
                       description="Specify input fasta file and preferences")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.fasta)")
options.add_option("-p","--pref",dest="preferences",
                   help="Ordering and revcomp information (.txt)")
#reverse complement of a string
def revcomp(s):
    sc=""
    for i in range(1,len(s)+1):
        if s[-i]=="T":sc+="A"
        if s[-i]=="A":sc+="T"
        if s[-i]=="G":sc+="C"
        if s[-i]=="C":sc+="G"
    return sc

#read fasta file, return list of contigs
def read_sequences(filename):
    f=open(filename)
    lines=f.readlines()
    f.close()
    titles=[]
    myseq=""
    for i in lines[1:]:
        if i[0]==">":
            titles.append(myseq)
            myseq=""
        else:
            myseq+=(i.upper()).replace("\n","")
    #don't forget to grab the last seq
    titles.append(myseq)
    return titles

#write a new fasta file with the appropriate combination of contigs
def write_contigs(c,order,rev,filename):
    thisfile=filename[:-6]+"_COMPLETE.fasta"
    f=open(thisfile,'w')
    
    f.write(">"+filename[:-6]+".complete\n")
    for i in range(0,len(c)-1):
        
        if rev[i]==1:
            f.write(c[order[i]]+"NNNNNNNNNNNNNNNNNNNNNNNNN")
        else:
            f.write(revcomp(c[order[i]])+"NNNNNNNNNNNNNNNNNNNNNNNNN")
    #we don't want trailing N's after the last contig
    if rev[-1]==1:
        f.write(c[order[-1]])
    else:
        f.write(revcomp(c[order[-1]]))
                
    f.close()

#read the preferences file, return two lists
def read_pref(preffile):
    f2=open(preffile)
    f2lines=f2.readlines()
    f2.close()
    o=[int(x)-1 for x in f2lines[0].split()]
    r=[int(x) for x in f2lines[1].split()]
    return o,r

def main():
    opts, args = options.parse_args()
    allcontigs=read_sequences(opts.inputfile)
    ordering, rc=read_pref(opts.preferences)
    write_contigs(allcontigs,ordering,rc,opts.inputfile)

if __name__ == '__main__':
    main()
