# DS 12/19/19
# Download SAM files from SRA, convert to BAM, remove SAM

import multiprocessing
from subprocess import call
import sys
from optparse import OptionParser


options = OptionParser(usage='%prog -i [inputfile] -o [outputdir]',
                       description="Specify input file(.txt) and output directory (for .bam) ")

options.add_option("-i","--input",dest="inputfile",
                   help="input file with SRA accession numbers")
options.add_option("-o","--output",dest="outputdir",
                   help="output directory for downloaded files")

def process_sample(sample_input):
	samplename = sample_input[0]
	outdir = sample_input[1]
	outsam = outdir+"/"+samplename+".sam"
	outbam = outdir+'/'+samplename+".bam"
	dlcommand = "sam-dump "+ samplename+" --output-file "+outsam
	convertcommand = "samtools view -S -b "+outsam+" > "+outbam
	removecommand = "rm "+outsam
	
	call(dlcommand, shell=True)
	call(convertcommand, shell=True)
	call(removecommand, shell=True)

if __name__ == '__main__':
	opts, args = options.parse_args()
	inputfile = opts.inputfile
	outputdir = opts.outputdir
	
	with open(inputfile) as f:
		samplenames = [i.strip() for i in f.readlines()]

	nthreads = 50
	# use a list of lists to iterate through
	sample_inputs = [[i, outputdir] for i in samplenames]
	
	p = multiprocessing.Pool(nthreads)
	try: 
		p.map(process_sample, sample_inputs)
	except AttributeError:
		print(sys.exc_info()[0])
		p.terminate()