import sys

def getUniqueMethylationSitesFromIntersectionsFileOfProtein(file_path) -> set:
    file = open(file_path, 'r')
    lines = file.readlines()
    uniqueMethylationSites = set()
    for line in lines:
        line = line.split("\t")
        uniqueMethylationSites.add((str(line[0]), str(line[1]),str(line[8]).rstrip()))
    file.close()
    return uniqueMethylationSites

def getStatusCode(status) -> int:
    if status == "positive":
        return 1
    elif status == "negative":
        return 0
    else:
        raise ValueError('Status-Code for protein-label not *positive* or *negative*. Wrong format.')
        
def getEncodingCode(encoding) -> str:
    if encoding == "oneHot":
        return 0
    elif encoding == "methylationRate":
        return 1
    else:
        raise ValueError('Encoding-Code not *oneHot* or *methylationRate*. Wrong format.')
        
def createProteinSequenceEncodingAnnotation(file_path, status, encodingCode, cutoff, methylationSites) -> None:
    chromosome = ""
    seqBeg = 0
    seqEnd = 0
    label = getStatusCode(status)
    encodingCode = getEncodingCode(encodingCode)
    methylationRate = 0
    
    encodedSequences = open("p_mr.fasta", "w")
    inputSequences = open(file_path, 'r')
    lines = inputSequences.readlines()
    
    for line in lines:
        if line[0] == ">":
            chromosome = line[1:].partition(":")[0].strip()
            seqList = (line.partition(":")[2]).partition("-")
            seqBeg = int(seqList[0])
            seqEnd = int(seqList[2])
            encodedSequences.write(">" + str(label) + "\n")
            # print(">" + str(label), end="\n") - debugging
        else:
            encoding = "0" * 200
            for methylationSite in methylationSites:
                position = int(methylationSite[1])
                if seqBeg <= position < seqEnd and chromosome == methylationSite[0]:
                    if encodingCode == 0:
                        encoding = encoding[:position-seqBeg] + "1" + encoding[position-seqBeg + 1:]
                    elif encodingCode == 1:
                        rate = str(int(float("{:.1f}".format(float(methylationSite[2])))*10))
                        # rate = int(rate) if int(rate) >= int(cutoff) else 0
                        # uncomment for cut-off methylation-rate encoding
                        encoding = encoding[:position-seqBeg] + str(rate) + encoding[position-seqBeg + 1:]
            encodedSequences.write(line + encoding + "\n")    
            # print(line + encoding, end="\n") - debugging
    encodedSequences.close()
    inputSequences.close()
    return

def createProteinSequnceEncodingHighMethylationRateOnly(file_path, status, encodingCode, cutoff, methylationSites) -> None:
    chromosome = ""
    seqBeg = 0
    seqEnd = 0
    label = getStatusCode(status)
    encodingCode = getEncodingCode(encodingCode)
    methylationRate = 0
    
    encodedSequences = open("fileName.fasta", "w")
    inputSequences = open(file_path, 'r')
    lines = inputSequences.readlines()
    output = ""
    
    for line in lines:
        if line[0] == ">":
            chromosome = line[1:].partition(":")[0].strip()
            seqList = (line.partition(":")[2]).partition("-")
            seqBeg = int(seqList[0])
            seqEnd = int(seqList[2])
            output = output + (">" + str(label) + "\n")
        else:
            encoding = "0" * 200
            for methylationSite in methylationSites:
                position = int(methylationSite[1])
                if seqBeg <= position < seqEnd and chromosome == methylationSite[0]:
                    if encodingCode == 0:
                        encoding = encoding[:position-seqBeg] + "1" + encoding[position-seqBeg + 1:]
                    elif encodingCode == 1:
                        rate = str(int(float("{:.1f}".format(float(methylationSite[2])))*10))
                        rate = int(rate) if int(rate) >= int(cutoff) else 0
                        encoding = encoding[:position-seqBeg] + str(rate) + encoding[position-seqBeg + 1:]
            
            output = output + (line + encoding + "\n")
            if encoding != ("0" * 200):
                encodedSequences.write(output)
            output = ""
    encodedSequences.close()
    inputSequences.close()
    return

def createProteinSequenceBaseline(file_path, status) -> None:
    label = getStatusCode(status)
    encodedSequences = open("pos_baseline.fasta", "w")
    inputSequences = open(file_path, 'r')
    lines = inputSequences.readlines()
    for line in lines:
        if line[0] == ">":
            encodedSequences.write(">" + str(label) + "\n")
        else:
            encodedSequences.write(line)    
    encodedSequences.close()
    inputSequences.close()
    return

if __name__ == "__main__":
    bedFile = "bedfile.bed"
    intersectionSet = getUniqueMethylationSitesFromIntersectionsFileOfProtein(bedfile)

    fastaFile = "fastafile.fasta"
    status = "negative"
    encoding = "methylationRate"
    cutOff = 0
    createProteinSequnceEncodingHighMethylationRateOnly(fastaFile, status, encoding, cutOff, intersectionSet)
