## Add annotation to RNAseq data

Use script ```add_annotation_RNAseq.py``` as follows:

```
python add_annotation_RNAseq.py -i [inputdirectory] -a [annotationfilename]
```

Where ```inputdirectory``` is the directory that contains all RNAseq data files (usually fold change expression) that will be annotated. ```annotationfilename``` is the file that contains the annotations to be added. For TIGR4 we use the file ```KEGG_and_Annot.csv```.
Make sure the data files all have a column called ```Gene``` with the old locus tags, and the annotation file has the column ```TIGR4.old```. 

Once you run this, each file in ```inputdirectory``` should have a corresponding ```[filename]_Merge.csv``` with the RNAseq data and added annotations. 

Required packages: 
```pandas```

----------------
## Get Normalization Genes for Tn-Seq Analysis

This program takes an input GBK file and writes a list of noralization genes.
Normalization genes for TnSeq are those that have the product annotation "Transposase" or "Mobile element"

Use script ```get_norm_genes.py``` as follows:
```
python get_norm_genes.py -i [strain].gbk -o [outputfile] (optional: --old)
```

Use the flag ```--old``` to retrieve old locus tags instead of new ones. 

----------------
## Filter background mutations

Filter population mutations by comparing against adaptation experiments done in CDM (control for background mutations).
Retain only mutations that are present in the experimental condition.

The experiment sheet provides a list of gd files to be analyzed. It has two columns labeled "File" and "Group". The entries under file should correspond to the .gd file names (.gd extension included), and Group could be either "E" (experimental) or "C" (control).

The genbank reference genome is used for adding locus tags to mutations if they appear in a gene. 

The filtering is done by keeping all mutations that satisfy the following two conditions:
1. The control populations do not have >(lower cutoff) frequency
2. There is at least one experimental population that has >(higher cutoff) frequency

Example usage: 
```
python filter_gd.py -i testdir -s Filter_Kan.csv -g NC_012469.gbk -o testout_10_50.csv -l 10 -u 50
```
Required packages: 
```pandas```
```Biopython```
