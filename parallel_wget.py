# DS 01/31/19
# Download files from ftp

import multiprocessing
from subprocess import call
import sys
from optparse import OptionParser


options = OptionParser(usage='%prog -i [inputfile] -o [outputdir]',
                       description="Specify input file(.txt) and output directory (for .bam) ")

options.add_option("-i","--input",dest="inputfile",
                   help="input file with SRA accession numbers")

def process_sample(sample_input):
	samplename = sample_input[0]
	dlcommand = "wget "+ samplename	
	call(dlcommand, shell=True)

if __name__ == '__main__':
	opts, args = options.parse_args()
	inputfile = opts.inputfile
	
	with open(inputfile) as f:
		samplenames = [i.strip() for i in f.readlines()]

	nthreads = 50
	# use a list of lists to iterate through
	sample_inputs = [[i] for i in samplenames]
	
	p = multiprocessing.Pool(nthreads)
	try: 
		p.map(process_sample, sample_inputs)
	except AttributeError:
		print(sys.exc_info()[0])
		p.terminate()
