## Add annotation to RNAseq data

Use script ```add_annotation_RNAseq.py``` as follows:

```
python add_annotation_RNAseq.py -i [inputdirectory] -a [annotationfilename]
```

Where ```inputdirectory``` is the directory that contains all RNAseq data files (usually fold change expression) that will be annotated. ```annotationfilename``` is the file that contains the annotations to be added. For TIGR4 we use the file ```KEGG_and_Annot.csv```.
Make sure the data files all have a column called ```Gene``` with the old locus tags, and the annotation file has the column ```TIGR4.old```. 

Once you run this, each file in ```inputdirectory``` should have a corresponding ```[filename]_Merge.csv``` with the RNAseq data and added annotations. 

Required packages: ```pandas```

----------------
