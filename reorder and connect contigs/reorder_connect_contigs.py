##Defne Surujon
##10 June 2016

##Read a fasta file, and reorder the sequences based on a preferences file

from optparse import OptionParser

options = OptionParser(usage='%prog input',
                       description="Specify input text file")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.txt) that contains strain names, reference genomes and their sizes")

#calculate the number of Ns coming after contig i
def numberN(order,c,cs,ce,gs,ge,f,i):
    n=0
    if f[i]==1 and f[i+1]==1:
        n=(gs[i+1]-ge[i])-cs[i+1]-(len(c[order[i]])-ce[i])
    elif f[i]==0 and f[i+1]==1:
        n=(gs[i+1]-gs[i])-cs[i]-cs[i+1]
    elif f[i]==1 and f[i+1]==0:
        n=(ge[i+1]-ge[i])-(len(c[order[i]])-ce[i])-(len(c[order[i+1]])-ce[i+1])
    elif f[i]==0 and f[i+1]==0:
        n=(ge[i+1]-gs[i])-cs[i]-(len(c[order[i+1]])-ce[i+1])

    return n

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
def write_contigs(c,order,cstart,cend,gstart,gend,rev,filename,gen):
    thisfile=filename+"_COMPLETE.fasta"
    f=open(thisfile,'w')

    f.write(">"+filename+".complete\n")
    #leading N's:
    if rev[0]==1:
        f.write("N"*(gstart[0]-cstart[0]))
    else:
        f.write("N"*(gstart[0]+cend[0]-len(c[order[0]])))
    nums=0
    for i in range(0,len(c)):
        
        if rev[i]==1:
            numslast=(gen-gend[i])-(len(c[order[i]])-cend[i])
            f.write(c[order[i]])
        elif rev[i]==0:
            numslast=(gen-gend[i])-cstart[i]
            f.write(revcomp(c[order[i]]))
        
        if i==len(c)-1:nums=numslast
        else:nums=numberN(order,c,cstart,cend,gstart,gend,rev,i)
        f.write("N"*max(nums,25))

    f.close()

#read the preferences file, return two lists
def read_pref(preffile):
    f2=open(preffile)
    f2lines=f2.readlines()
    f2.close()
    o=[int(x)-1 for x in f2lines[0].split()]
    cstart=[int(x) for x in f2lines[1].split()]
    cend=[int(x) for x in f2lines[2].split()]
    gstart=[int(x) for x in f2lines[3].split()]
    gend=[int(x) for x in f2lines[4].split()]
    r=[int(x) for x in f2lines[5].split()]

    return o,cstart,cend,gstart,gend,r

def process_fasta(inputfile, preferences,genome):
    allcontigs=read_sequences(inputfile)
    ordering, cstart,cend,gstart,gend,rc=read_pref(preferences)
    write_contigs(allcontigs,ordering,cstart,cend,gstart,gend,rc,inputfile,genome)

def main():
    opts, args = options.parse_args()

    f=open(opts.inputfile)
    lines=f.readlines()
    f.close()
    prefs=[x for x in lines[0].split()]
    genomelens=[int(x) for x in lines[2].split()]

    for i in range(0,len(prefs)):
        fastaname=prefs[i]+".gff.fasta"
        prefname=prefs[i]+"_results_positions.txt"
        process_fasta(fastaname,prefname,genomelens[i])

if __name__ == '__main__':
    main()
