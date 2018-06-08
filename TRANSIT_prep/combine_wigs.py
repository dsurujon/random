## DS 6/6/18
## Combine wig files into one 
## 

import glob
from optparse import OptionParser


options = OptionParser(usage='%prog -i [input prefix]',
                       description="Specify input file prefix")

options.add_option("-i","--infile",dest="inputfile",
                   help="input file prefix")
				   
def main():
	opts, args = options.parse_args()
	infile = opts.inputfile
	outfile = infile+"_ALL.wig"
	allinputfiles = glob.glob(infile+"*")
	with open(allinputfiles[0]) as f:
		wiglines = [line.split() for line in f.readlines()]
	for i in range(1,len(allinputfiles)):
		with open(allinputfiles[i]) as f2:
			wiglines2 = [line.split() for line in f2.readlines()]
			for j in range(1,len(wiglines2)):
				wiglines[j][1] = str(int(wiglines[j][1])+int(wiglines2[j][1]))
	with open(outfile,'w') as g:
		for l in wiglines:
			g.write("\t".join(l)+"\n")
			
if __name__ == '__main__':
    main()

