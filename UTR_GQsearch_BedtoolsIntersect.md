

# 1. Creating a list of all UTRs

I started by using Rockhopper predictions of transcription/translation start/stop sites to define the genomic locations. I selected only the UTRs that were at least 10 bp long, becasue the is the minimum size of a 2G-GQ sequence. I did all of these data manipulations in R, then saved the dataframe as a csv file.

```
library(IRanges)

## Importing data from Rockhopper RNA-seq analysis
RNA_seq <- read.table("~/R/R_projects/G-Quadruplexes/data/Transcript_boundaries.txt", header=TRUE, sep = ",")
## Removes rows with no entries
RNA_seq2 <- na.omit(RNA_seq)
## Subsets the data to include only the columns we're interested in
UTRs <- RNA_seq2[c("Gene", "Transcription.Start", "Translation.Start", "Strand")]
## Re-orders the data columns in a way that makes sense
UTRs <- UTRs[, c("Transcription.Start", "Translation.Start", "Strand", "Gene")]
## New variable created for dataframe containing UTRs of length >= 10
newUTRs <- data.frame(Transcription.Start = numeric(0), Translation.Start = numeric(0), Strand = character(0), Gene = str(0), y = numeric(0))

## Creates a new vector, "UTR_length", that contains the differences between transcription start and translation start
for(x in UTRs){
  UTR_length = c(abs(UTRs[x, "Translation.Start"] - UTRs[x, "Transcription.Start"]))
}

## Adds the column "UTR_length" to the UTRs dataframe
UTRs <- cbind(UTRs, UTR_length)

## Goes through the UTRs dataframe row by row and adds rows that have length >= 10 to the dataframe newUTRs
for(x in 1:nrow(UTRs)){
  if(UTRs[x, "UTR_length"] >= 10){
    newUTRs <- rbind(newUTRs, UTRs[x,])
  }
}

## Creates new empty dataframes plus_UTRs and neg_UTRs to store UTRs found on the positive and negative strands separately
plus_UTRs <- data.frame(Transcription.Start = numeric(0), Translation.Start = numeric(0), Strand = character(0), Gene = str(0), y = numeric(0))
neg_UTRs <- data.frame(Transcription.Start = numeric(0), Translation.Start = numeric(0), Strand = character(0), Gene = str(0), y = numeric(0))

## Goes through the dataframe newUTRs and adds all rows with genes on the positive strand to plus_UTRs and all the ones on the negative strand to neg_UTRs
for(x in 1:nrow(newUTRs)){
  if(newUTRs[x, "Strand"] == "+"){
    plus_UTRs <- rbind(plus_UTRs, newUTRs[x,])
  }else{
    neg_UTRs <- rbind(neg_UTRs, newUTRs[x,])
  }
}

## Writing both dataframes to their own csv files
write.csv(plus_UTRs, file = "plus_UTRs.csv",row.names=FALSE)
write.csv(neg_UTRs, file = "neg_UTRs.csv",row.names=FALSE)
```
# 2. Creating a list of all GQs, their genomic postions, and their strand

Used the re.finditer() function to find all non-overlapping instances of the search pattern, and to extract information about them. This information was stored in a csv file that could later be manipulated. 
```
#!/usr/bin/python
import sys, fileinput, re, csv, ushuffle
sequence = ""
file = fileinput.input()

for line in file:
    if line[0] == ">":
        title = line[1:]
    else:
        sequence = sequence + line
sequence = sequence.upper().replace("\n", "")

p1 = re.compile("GGG[ATGCN]{1,7}GGG[ATGCN]{1,7}GGG[ATGCN]{1,7}GGG")
p2 = re.compile("CCC[ATGCN]{1,7}CCC[ATGCN]{1,7}CCC[ATGCN]{1,7}CCC")

with open ('/home/gradstd4/GQ_project/Outputs/sven_shuff_GQ.csv', 'wb') as file:
    writer = csv.writer(file)
    writer.writerow(['Chromosome', 'Start', 'End', 'Strand'])

for m in p1.finditer(shuff):
    with open('/home/savannah/GQ_project/Outputs/sven_shuff_GQ.csv', 'a') as file:
        writer = csv.writer(file)
        writer.writerow(['sven', m.start(), m.end(), '+'])
for m in p2.finditer(shuff):
    with open('/home/savannah/GQ_project/Outputs/sven_shuff_GQ.csv', 'a') as file:
        writer = csv.writer(file)
        writer.writerow(['sven', m.start(), m.end(), '-'])
```

# 3. Converting these to bed files

In bed files, the format for the column order is:
1. Chromosome
2. Start
3. End
4. Name
5. Score
6. Strand
...

For our purposes, we need to ahve data in rows 1-4 and 6. The reason we need to have separate files for UTRs on the positive and on the negative strand is that with UTRs on the negative strand, the transcription start site will actually be less than the translation start site, and we can't have a negative range. There a couple simple formatting steps that can be done in Excel. First, add a column to all data files with the chrosome name for each entry. Then, for the negative UTRs, place the translation start column before the transcription start column. Now, we will use a few simple commands to transform this into a bed file.

Remove the header:
```
tail -n +2 data.csv > noheader_data.csv
```
Translate ',' to tabs:
```
cat noheader_data.csv | tr ',' '\t' > noheader_data.bed
```

# 4. Running bedtools intersect

This handy program finds regions of overlap between two files and returns the data in several different useful formats depending on your needs. Read more about bedtools interesect [here](http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html).
```
bedtools interesect -a <GQ file> -b <UTR file> -wa -wb > output.txt
```
The output file was then manipulated in Excel. The output is a little messy, so I organized it by sorting it based on the second column. This put all of the entries for UTRs below the entries for GQs in the appropriate order. Then I found the middle and moved the UTR entries so that they were next to the GQ entries. I added header names. I used VLOOKUP to find expression level of the genes and the product of the gene so that I could sort them based on expression level and start my analysis with the most highly expressed genes. 