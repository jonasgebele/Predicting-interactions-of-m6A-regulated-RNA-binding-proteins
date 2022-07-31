# Predicting Interactions of m6A regulated RNA-Binding Proteins

The rising success of mRNA vaccines has had a major role in showing the importance of chemical
modifications and how they can influence interactions drastically. We put our focus on a specific
modification and how it can affect interactions between RNA and proteins. It has been shown for the
often-occurring m6A modification that it seems to play a part in regulating the interaction between
RNA and proteins, specifically RNA-binding proteins (RBPs). As this is a well-researched topic, it
provides us with many crosslinking immunoprecipitation (CLIP) datasets for RBPs binding sites and the
methylation of those, which we used to train a convolutional neural network (CNN). We then assessed
different encodings of methylation, models and dataset distribution regarding their predictions on m6A
regulated RBPs.

![RBP_methylation](https://user-images.githubusercontent.com/41921050/182041627-81cfbd84-915a-4922-bf8d-f140cdb0ef19.png)

## Data Processing
### PAR-CLIP Data
For the protein-RNA interactions we used the pre-processed PAR-CLIP dataset from Mukherjee et al. 2019, which analyzed 66 proteins in HEK293 cells. The methylation data came from a miCLIP dataset that has been also done on HEK293 cells.

### Generate Positive Samples
As a fist step we need to generate positive and negative samples, that we can use as classified inputs to train our model on.
The positive samples are available in the parclip dataset but we have to extend them to generate sequences with a length of 200 amino acids (AA) as the length of the original bed files would vary.
```
$ cat proteinFile.bed | awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' \
> positives.bed
```
Command to replace the spaces with tabs (Bedtools standard).
```
$ sed 's/ /\t/g' positives.bed > positives.bed
```
Dynamically determine the center of the sequence and extend the length from that position in both direction with a length of 100.
This ensures that the original sequence is centered around the chromStart-chromEnd.

### Generate Negative Samples
Negative samples of type 1 for the model are generated by sampling randomly in a window around the positives e.g. at distance 200nt from both left and from right.</br>
From the right:
```
$ cat proteinFile.bed | bedtools slop -i stdin -g XPO5 -l -250 -r 250 | \
awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' \
> negativesType1Right.bed
```
```
$sed 's/ /\t/g' negativesType1Right.bed > negativesType1Right.bed
```
From the left:
```
$ cat proteinFile.bed | bedtools slop -i stdin -g XPO5 -l 250 -r -250 | \
awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' \
>  negativesType1Left.bed
```
```
$ sed 's/ /\t/g' negativesType1Left.bed > negativesType1Left.bed
```
Negative samples of type 2 are positives of other RBPs.
```
$ cat proteinFile.bed | awk '$3 = $3 + 100 - sprintf("%.0f", (($3 - $2)/2)), $2 = $3 - 200' \
>  negativesType2.bed
```
```
$ sed 's/ /\t/g' negativesType2.bed > negativesType2.bed
```
Remove duplicate lines from a file assuming you don't mind that lines are sorted.
```
$ cat -n positives.bed| sort -uk2 | sort -nk1 | cut -f2- > positives.fasta
```
### Gencode Data
```
$ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.primary_assembly.genome.fa.gz
$ gzip -d GRCh38.primary_assembly.genome.fa.gz
```
Intersections done via bedtools.
```
$ bedtools getfasta -fi gencode_data/GRCh38.primary_assembly.genome.fa -bed positives.bed > positives.fasta
```
### Generate Data-Set out of MI-CLIP Data and processed PAR-CLIP Data
Generate labeled and encoded sequence-data out of modified positives, negatives and different encodings.
```
$ python EncodingPreprocessing.py miclip.bed positives.fasta "positives" "methylationRate" 0
```

## Network Architecture
![image](https://user-images.githubusercontent.com/41921050/182041673-99f8467b-1ff8-4a57-8113-4ba017c1cfdf.png)


