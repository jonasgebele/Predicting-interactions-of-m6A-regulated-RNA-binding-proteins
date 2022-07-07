def getUniqueMethylationSitesFromIntersectionsFileOfProtein(file_path) -> set:
    # file with intersection with the protein needs to be in the same directory
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
        
def createProteinSequenceEncodingAnnotation(file_path, status, encodingCode, methylationSites) -> None:
    chromosome = ""
    seqBeg = 0
    seqEnd = 0
    label = getStatusCode(status)
    encodingCode = getEncodingCode(encodingCode)
    methylationRate = 0
    
    encodedSequences = open("positives.fasta", "w")
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
                        encoding = encoding[:position-seqBeg] + rate + encoding[position-seqBeg + 1:]
                        
            encodedSequences.write(line + encoding + "\n")    
            # print(line + encoding, end="\n") - debugging
    
    encodedSequences.close()
    inputSequences.close()
    return
    
if __name__ == '__main__':
    path = 'GSE63753_hek293.abcam.CIMS.m6A.9536.bed'
    s = getUniqueMethylationSitesFromIntersectionsFileOfProtein(path)

    path2 = "SRR488727_positive_extended_fasta_unique - Kopie"
    createProteinSequenceEncodingAnnotation(path2, "positive", "methylationRate", s)
