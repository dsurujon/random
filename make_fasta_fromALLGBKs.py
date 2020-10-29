
import os
import numpy
from Bio import SeqIO
from optparse import OptionParser

options = OptionParser(usage='%prog -i [inputdir] -o [outputfile]',
                       description="Specify input directory (for .gbk) and output file (for .fasta)")

options.add_option("-i","--indir",dest="inputdir",
                   help="input directory of gff files")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output fasta file")


def make_fasta_from_GBKs(thisdir,outfilename):
    n=1
    gbks = os.listdir(thisdir)
    with open(outfilename,'w') as f:
        for gbkfile in gbks:
            fullfilename = os.path.join(thisdir,gbkfile)
            mystrain = SeqIO.read(fullfilename,"genbank")
            strainname = mystrain.name
            for thisfeature in mystrain.features:
                if thisfeature.type == "CDS":
                    try:
                        aaseq = thisfeature.qualifiers['translation'][0]
                        try:
                            m=aaseq.index('/translation=')
                            aaseq = aaseq[m+13:]
                        except ValueError:
                            pass
                        locustag = thisfeature.qualifiers['locus_tag'][0]
                        annot = thisfeature.qualifiers['product'][0]
                        headerlist = ['>[',str(n),']tvo|',strainname,'.[gene=',locustag,'].[protein=',annot,']\n']
                        f.write(''.join(headerlist))
                        f.write(aaseq+'\n')
                        n += 1

                    except KeyError: pass

def make_fasta_from_contigGBKs(thisdir,outfilename):
    n=1
    gbks = os.listdir(thisdir)
    with open(outfilename,'w') as f:
        for gbkfile in gbks:
            fullfilename = os.path.join(thisdir,gbkfile)
            # the assembled chromosome is not available, so we need to
            # parse into a list of contigs
            mystrain = list(SeqIO.parse(fullfilename,"genbank"))

            strainname = gbkfile[:15]

            allfeatures = []
            for contig in mystrain:
                allfeatures = allfeatures+contig.features
                
            for thisfeature in allfeatures:
                if thisfeature.type == "CDS":
                    try:
                        aaseq = thisfeature.qualifiers['translation'][0]
                        try:
                            m=aaseq.index('/translation=')
                            aaseq = aaseq[m+13:]
                        except ValueError:
                            pass
                        locustag = thisfeature.qualifiers['locus_tag'][0]
                        annot = thisfeature.qualifiers['product'][0]
                        headerlist = ['>[',str(n),']tvo|',strainname,'.[gene=',locustag,'].[protein=',annot,']\n']
                        f.write(''.join(headerlist))
                        f.write(aaseq+'\n')
                        n += 1

                    except KeyError: pass

def main():
    opts, args = options.parse_args()
    inputdir = opts.inputdir
    outputfile = opts.outputfile
    make_fasta_from_contigGBKs(inputdir, outputfile)

if __name__ == '__main__':
    main()
