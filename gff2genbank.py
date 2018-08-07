## DS 8/7/18
## Convert a GFF3 file into genbank

from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from optparse import OptionParser
import re
import os


options = OptionParser(usage='%prog -i [inputdir] -o [outputdir]',
                       description="Specify input (for .gff) and output (for .gbk) directories")

options.add_option("-i","--indir",dest="inputdir",
                   help="input directory of gff3 files")
options.add_option("-o","--outdir",dest="outputdir",
                   help="output directory of genbank files")

myfile = "D:\\defne\\Documents\\BC\\TVO\\clustering\\GFF20_TESTS\\DS-00-000-0.gff"

def get_genome(inputfile):
    with open(inputfile) as f:
        lines = f.readlines()
    sepline = lines.index("##FASTA\n")
    genome_id = lines[sepline+1][1:].strip()

    genome_features = []
    for i in range(2,sepline):
        splitline = lines[i].split()
        startloc = int(splitline[3])
        endloc = int(splitline[4])
        strand = splitline[6]
        locustag = splitline[8].split('.')[1]
        thisfeature = SeqFeature(FeatureLocation(start=startloc, end=endloc),
                                 type=splitline[2],
                                 qualifiers = {'locus_tag':[locustag]})
        genome_features.append(thisfeature)
    
    genome_seq = ""
    
    for i in range(sepline+2, len(lines)):
        genome_seq = genome_seq + lines[i].strip()
    return(genome_id, genome_seq, genome_features)

def make_gbk_from_gff(infile, outfile):
    g_id, g_seq, g_features = get_genome(infile)
    myseq = Seq.Seq(g_seq, IUPAC.unambiguous_dna)
    myrecord = SeqRecord(myseq, id = g_id, name = g_id)
    for feature in g_features:
        myrecord.features.append(feature)
    SeqIO.write(myrecord, outfile, 'genbank')

def main():
    opts, args = options.parse_args()
    inputdir = opts.inputdir
    outputdir = opts.outputdir

    gfffiles = os.listdir(inputdir)
    gfffiles = [g for g in gfffiles if g[-4:]==".gff"]

    for infile in gfffiles:
        print("processing ", infile)
        
        try:
        	fullinputfile = os.path.join(inputdir, infile)
        	outfile = '.'.join(infile.split('.')[:-1])+".gbk"
        	fulloutputfile = os.path.join(outputdir, outfile)

        	make_gbk_from_gff(fullinputfile, fulloutputfile)
        except OverflowError: pass
        



if __name__ == '__main__':
    main()
