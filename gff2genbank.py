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
import copy


options = OptionParser(usage='%prog -i [inputdir] -o [outputdir]',
                       description="Specify input (for .gff) and output (for .gbk) directories")

options.add_option("-i","--indir",dest="inputdir",
                   help="input directory of gff3 files")
options.add_option("-o","--outdir",dest="outputdir",
                   help="output directory of genbank files")

stranddict = {'+':1, '-':-1}

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
        featstrand = stranddict[splitline[6]]
        locustag = splitline[8].split('=')[1]
        thisfeature = SeqFeature(FeatureLocation(start=startloc, end=endloc, strand = featstrand),
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
    f1 = SeqFeature(FeatureLocation(0, len(myseq)), type="source")
    myrecord.features.append(f1)
    for feature in g_features:
        feature.qualifiers['transl_table']=[11]
        if feature.location.strand == 1:
            aaseq = myseq[feature.location.start-1:feature.location.end-1].translate()
        else:
            aaseq = myseq[feature.location.start-1:feature.location.end-1].reverse_complement().translate()
        feature.qualifiers['translation']=[str(aaseq)]
        feature.type="CDS"
        feature.qualifiers['product'] = feature.qualifiers['locus_tag']
        
        genefeature = copy.deepcopy(feature)
        genefeature.type = 'gene'
        genefeature.qualifiers = {'locus_tag':feature.qualifiers['locus_tag']}
        
        myrecord.features.append(genefeature)
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
        	outfile = outfile.replace('-', '_')
        	fulloutputfile = os.path.join(outputdir, outfile)

        	make_gbk_from_gff(fullinputfile, fulloutputfile)
        except OverflowError: pass
        



if __name__ == '__main__':
    main()
