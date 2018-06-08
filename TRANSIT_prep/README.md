These scripts help prepare files for [TRANSIT](http://saclab.tamu.edu/essentiality/transit/), a software package for gene essentiality analysis. 
-------------------
Requirements:
1. Annotated genome of the organism (as a .gbk file)
2. Map files from Tn-Seq experiments

-------------------

1. Make prot_table from a reference genome (.gbk file format) using the make_prot_table script. This step prepares the coding sequence annotations as the appropriately formatted protein table for TRANSIT. You only need to do this once per strain. 

```python make_prot_table.py -i [referencegenome].gbk -o [strainname].prot_table```

2. Make an empty .wig file from a reference genome (.gbk). This step prepares a table of available and unique insertion sites from a reference genome. You only need to do this once per strain. (NB: this step is quite time consuming, since it checks whether each insertion site is within a unique sequence region. The insertion sites that are in repeated regions are discarded)

```python make_empty_wig.py -g [referencegenome].gbk -o [strainname].wig```

3. Add the observed insertions to the empty wig file for each experiment (each individual map file). 

```python populate_empty_wig.py -i [emptywigfile].wig -m [mapfile].map -o [finalwigfile].wig```

4. [Optional] Combine multiple .wig files. If there are multiple libraries prepared for the same experiment, it is possible to aggregate these together. The files to be aggregated should have the same prefix and be in the same directory. 

```python combine_wigs.py -i [file prefix]```

5. [Optional] Filter out insertion sites with few (<k) insertions from all .wig files in a directory. 

```python filter_insertions.py -i [input directory] -k [minimum insertion count]```
