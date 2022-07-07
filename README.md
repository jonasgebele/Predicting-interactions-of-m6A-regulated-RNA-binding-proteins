# Predicting-interactions-of-m6A-regulated-RNA-binding-proteins
MOTIVATION: m6A can affect protein-RNA interactions
METHODS: integrate m6A and RNA sequence data into a model to predict binding sites
INTERPRETATION: determine what drives binding

## Pre-Processing

### Parclip Data
```console
pip install gdown
gdown --folder https://drive.google.com/drive/folders/1RRQYkRJWxy6C4dldDzd2WjEjUkeM6Qsy
tar -xf parclip_data.tar.gz
```

### Positive Samples
```console
cat SRR488727 | awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' > SRR488727_positive_extended
sed 's/ /\t/g'  SRR488727_positive_extended > SRR488727_positive_extended_tab
```
constant 100 & -200 is correct since the chromStart one is included and the chromEnd isn’t<br>
calculates difference, roundes off, and shifts cromEnd right with [100-difference] and shifts cromStart left starting from cromEnd of -200 to the left. This ensures that the sample is centered around the cromStart-cromEnd

### Negative Samples 1
```console
cat SRR488727 | bedtools slop -i stdin -g XPO5 -l -250 -r 250 | awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' > SRR488727_negative1_right_extended
sed 's/ /\t/g'  SRR488727_negative1_left_extended > SRR488727_negative1_left_extended_tab
```
200nt [or 201-depending on rounding] from the right
```console
cat SRR488727 | bedtools slop -i stdin -g XPO5 -l 250 -r -250 | awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' > SRR488727_negative1_left_extended
sed 's/ /\t/g'  SRR488727_negative1_right_extended > SRR488727_negative1_right_extended_tab
```
200nt [or 201-depending on rounding] from the left

### Negative Samples 2
```console
cat SRR488741 | awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' > SRR488741_positive_extended
sed 's/ /\t/g'  SRR488741_positive_extended > SRR488741_positive_extended_tab
```
