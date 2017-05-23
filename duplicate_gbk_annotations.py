## Defne Surujon
## May 2017

#input a genbank file with only CDS/tRNA/rRNA/ncRNA
#annotatoins, write a new file with all such features
#duplicated as "gene" features

import os
from optparse import OptionParser

options = OptionParser(usage='%prog input output',
                       description="Specify input gbk file and output gbk file")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk)")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output file (.csv)")
				   

def make_new_gbk(infilename,outfilename):
	f=open(infilename)
	lines=f.readlines()
	f.close()
	
	tags=["CDS","tRNA","rRNA","ncRNA"]
	geneflags=["/locus_tag","/old_locus_tag","/gene"]\
	
	g=open(outfilename,'w')
	
	curr_line_ix=0
	curr_feature_start_ix=0
	curr_feature_end_ix=0
	inFeatures=False
	
	for i in range(0,len(lines)):
		curr_line_ix=i
		curr_line=lines[i]
		if "ORIGIN" in curr_line:
			inFeatures==False
		if inFeatures==False:
			g.write(curr_line)
		if "FEATURES" in curr_line:
			inFeatures=True
		#when we are in FEATURES
		elif inFeatures:
			#find the start of the feature
			if any([t in curr_line for t in tags]):
				curr_feature_start_ix=i
				curr_feature_end_ix=i+1
				#find the end of the feature
				while any([t in lines[curr_feature_end_ix] for t in tags])==False and "ORIGIN" not in lines[curr_feature_end_ix]:
					curr_feature_end_ix+=1
				genestartline=lines[curr_feature_start_ix].replace("     CDS             ","     gene            ")
				g.write(genestartline)
				for jj in range(curr_feature_start_ix+1,curr_feature_end_ix):
					if any([gf in lines[jj] for gf in geneflags]):
						g.write(lines[jj])
			g.write(curr_line)

	g.close()



def main():
	opts, args = options.parse_args()
	
	infile = opts.inputfile
	outfile = opts.outputfile
	
	make_new_gbk(infile,outfile)

if __name__ == '__main__':
    main()
		