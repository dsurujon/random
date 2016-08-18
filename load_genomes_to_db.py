from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
import os
from _mysql_exceptions import *

server = BioSeqDatabase.open_database(driver="MySQLdb", user="surujon",passwd = "rna314rulz", host = "prince.bc.edu", db="pangenome")


refpath="/cluster/sequencer/TVOLab/20Strains_withproteins/"
newpath="/cluster/sequencer/TVOLab/30_Selected_Strains/FinalGBK/"
referencestrains=[f for f in os.listdir(refpath) if f.endswith(".gbk")]
newstrains=[f for f in os.listdir(newpath) if f.endswith(".gbk")]

server.new_database("tvogenomes")
refdb=server["tvogenomes"]

for strain in referencestrains:
	myfile=refpath+strain
	handle=open(myfile,"r")
	print("parsing "+strain)
	refdb.load(SeqIO.parse(handle,"genbank"))
	handle.close()

server.new_database("newgenomes")
newdb=server["newgenomes"]

for strain in newstrains:
	try:
		myfile=newpath+strain
		handle=open(myfile,"r")
		print("parsing "+strain)
		newdb.load(SeqIO.parse(handle,"genbank"))
		handle.close()
	except IntegrityError as e:
		print(e)


server.commit()
