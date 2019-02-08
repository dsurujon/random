## Defne Surujon 06/08/18

# Remove insertions with <k reads from a wig file

from optparse import OptionParser
import re
import os


options = OptionParser(usage='%prog -i [input dir] -k [number of reads]',
                       description="Specify input directory and minimum number of reads")


options.add_option("-i","--infile",dest="inputdir",
                   help="input directory containing .wig files")
options.add_option("-k","--kmin",dest="kmin",
                   help="minimum number of reads")

				   
				   
def filter_wig(infile,k):
	with open(infile) as f:
		lines = [line.split() for line in f.readlines()]
		outfile = infile[:-4]+"_Filter.wig"
		with open(outfile,'w') as g:
			g.write("\t".join(lines[0])+"\n")
			#print(lines[1])
			for line in lines[1:]:
				if int(line[1])<k:
					line[1] = "0"
				g.write("\t".join(line)+"\n")

def main():
	opts, args = options.parse_args()
	inputdir = opts.inputdir
	k = int(opts.kmin)
	os.chdir(inputdir)
	wigfiles = os.listdir('./')
	for wigfile in wigfiles:
		print(wigfile)
		if wigfile[-4:]==".wig":
			filter_wig(wigfile,k)
	
if __name__ == '__main__':
    main()
