# taxamatcher
taxamatcher is an easy to use tool for assigning taxon IDs of reference taxonomy to the desired list of taxa or consensus lineages. This tool is useful when there is a need to assign taxon identifiers of reference taxonomy map file to a list of available taxa (i.e. list of consensus lineages from desired OTU table). Reference taxonomy can be 16S ribosomal RNA database like Greengenes and SILVA with corresponding map files like "97_otu_taxonomy.txt" and "tax_slv_ssu_132.txt". Latest versions of these files are readily available and can be downloaded from respective online websites. Moreover, after assignment of reference taxon identifiers, it is possible to easily obtain truncated trees from newick tree files of processed reference database taxonomy  

# Currently available taxonomic classifications
- GreenGenes
- SILVA

# Installation
    $ pip install taxamatcher
    
# Usage
    $ taxamatcher -r [gg or silva] -m [ref-map-file] -t [my-lineages].csv -o [output].csv
- [ref-map-file] - path to map file. Example is "97_otu_taxonomy.txt" for greengenes taxonomy
- [my-lineages] - path to target taxonomy. Target file must be CSV file with taxon IDs in first column and consensus lineages in second column
- [output].csv - path to output file. Output will be generated as CSV file

# Author
- Farid MUSA (farid.musa.h@gmail.com)
