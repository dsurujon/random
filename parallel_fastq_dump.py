# DS 10/04/20
# Download fastq files from SRA

import multiprocessing
from subprocess import call
import sys
from optparse import OptionParser


options = OptionParser(usage='%prog -i [inputfile]',
                       description="Specify input file(.txt) with accession list ")

options.add_option("-i","--input",dest="inputfile",
                   help="input file with SRA accession numbers")


def process_sample(sample_input):
	samplename = sample_input[0]
	dlcommand = "/store/home/surujon/SRA/sratoolkit.2.9.2-centos_linux64/bin/fastq-dump "+ samplename
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
