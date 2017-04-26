
#find string between two specified strings
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return (s[start:end])
    except ValueError:
        return ("")

#read fasta file, return dictionary of sequences
def read_sequences(filename):
    f=open(filename)
    lines=f.readlines()
    f.close()
    titles={}
    mytitle=""
    for i in lines:
        if i[0]==">":
            genename=find_between(i,"[gene=","]")
            mytitle=genename
            titles[mytitle]=""
        else:
            titles[mytitle]+=(i.upper()).strip()

    return titles

def get_wanted_proteins(outfile,maxprots,gen):
    f=open(outfile,'w')

    for i in range(1,gen):
        counter=0
        for j in range(1,20000):
            strnum='%02d' %i
            genenum='%05d' %(j*5)
            genename="SPS"+strnum+"_"+genenum
            try:
                f.write(">"+genename+"\n"+allseqs[genename]+"\n")
                counter+=1
                if counter>maxprots:
                    break
            except KeyError:pass
                #print (genename+" not found")
    f.close()


xStr="C:/Users/defne/Documents/BC/TVO/pangenome/50Str.Combo.fasta"
allseqs=read_sequences(xStr)

firstx=100
numgenomes=30
filename="C:/Users/defne/Documents/BC/TVO/pangenome/first"+str(firstx)+"_"+str(numgenomes)+"strains"+".fasta"
get_wanted_proteins(filename,firstx,numgenomes)


