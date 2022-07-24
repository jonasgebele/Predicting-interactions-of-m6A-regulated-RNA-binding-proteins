# Predicting-interactions-of-m6A-regulated-RNA-binding-proteins

## Data Processing
For the protein-RNA interactions we used the pre-processed PAR-CLIP dataset from Mukherjee et al. 2019, which analyzed 66 proteins in HEK293 cells. The methylation data came from a miCLIP dataset that has been also done on HEK293 cells
### Miclip Data
Available in data folder.

### Generate Positive Samples
As a fist step we need to generate positive and negative samples, that we can use as classified inputs to train our model on.
The positive samples are available in the parclip dataset but we have to extend them to generate sequences with a length of 200 amino acids (AA) as the length of the original bed files would vary.
```
cat proteinFile.bed | awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' \
> proteinFile_extendedSLength.bed
```
Command to replace the spaces with tabs (Bedtools standard).
```
sed 's/ /\t/g' proteinFile_extendedSLength.bed > proteinFile_extendedSLength.bed
```
Dynamically determine the center of the sequence and extend the length from that position in both direction with a length of 100.
This ensures that the original sequence is centered around the chromStart-chromEnd.

### Negative Samples 1
```
cat <proteinFile> | bedtools slop -i stdin -g XPO5 -l -250 -r 250 | \
awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' \
> bedfile-extended-ranges
```
```
sed 's/ /\t/g' bedfile-extended-ranges > bedfile-extended-ranges
```
200nt from the right
```
cat <proteinID> | bedtools slop -i stdin -g XPO5 -l 250 -r -250 | \
awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' 
>  bedfile-extended-ranges
sed 's/ /\t/g' bedfile-extended-ranges > bedfile-extended-ranges
```
200nt from the left

### Negative Samples 2
```
cat <proteinID> | awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' \
>  bedfile-extended-ranges
```
```
sed 's/ /\t/g' bedfile-extended-ranges >  bedfile-extended-ranges
```

### Genome Data
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz
gzip -d GRCh38.primary_assembly.genome.fa.gz
bedtools getfasta -fi GRCh38.primary_assembly.genome.fa -bed bedfile-extended-ranges > sequence-data
```
```
cat -n sequence-data | sort -uk2 | sort -nk1 | cut -f2- > sequence-data-unique
```
Remove duplicate lines from a file assuming you don't mind that lines are sorted.
### Methylation sites
```
bedtools intersect -a parclip_data/sequence-data \
-b miclip_data/GSE63753_hek293.abcam.CIMS.m6A.9536.bed \
> sequence-data-intersections
```

