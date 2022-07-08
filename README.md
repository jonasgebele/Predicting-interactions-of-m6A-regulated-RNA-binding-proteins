# Predicting-interactions-of-m6A-regulated-RNA-binding-proteins
MOTIVATION: m6A can affect protein-RNA interactions
METHODS: integrate m6A and RNA sequence data into a model to predict binding sites
INTERPRETATION: determine what drives binding

## Pre-Processing
### Miclip Data
Available in data folder.

### Parclip Data
Get via gdown.
```console
pip install gdown
gdown --folder https://drive.google.com/drive/folders/1RRQYkRJWxy6C4dldDzd2WjEjUkeM6Qsy
tar -xf parclip_data.tar.gz
```

### Generate Positive Samples
```console
cat <proteinID> | awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' > bedfile-extended-ranges
```
Handle tab and whitespace sperator.
```console
sed 's/ /\t/g' bedfile-extended-ranges > bedfile-extended-ranges
```
Calculates difference, roundes off, and shifts cromEnd right with [100-difference] and shifts cromStart left starting from cromEnd of -200 to the left.
This ensures that the sample is centered around the chromStart-chromEnd.

### Negative Samples 1
```console
cat <proteinID> | bedtools slop -i stdin -g XPO5 -l -250 -r 250 | awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' > bedfile-extended-ranges
```
```console
sed 's/ /\t/g' bedfile-extended-ranges > bedfile-extended-ranges
```
200nt from the right
```console
cat <proteinID> | bedtools slop -i stdin -g XPO5 -l 250 -r -250 | awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' >  bedfile-extended-ranges
sed 's/ /\t/g' bedfile-extended-ranges > bedfile-extended-ranges
```
200nt from the left

### Negative Samples 2
```console
cat <proteinID> | awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' >  bedfile-extended-ranges
```
```console
sed 's/ /\t/g' bedfile-extended-ranges >  bedfile-extended-ranges
```
