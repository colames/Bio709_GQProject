# Finding GQs located in UTRs

## Defining transcript boundaries

RNA-seq data generated from our lab was run through Rockhopper. This program takes fastq formated files and aligns them to the reference genome, then determines transcription start and stop sites as well as translation start ands stop sites. It also gives information about the expression level of the gene, but this is not important for this analysis.

I wrote a script that finds all UTRs using transcriptome data, then pulls out the sequences of these UTRs from the fasta genome file, and searches these sequences for GQ sequences. It then matches up these GQ sequences to a list of known GQ sequences to extract information on their genomic location and associated gene. 

The script requires as input three files: 
1. A fasta file of the reference genome.
2. Rockhopper output csv file.
3. Results from GQ search

For the fasta genome file, I removed the title line and stored the sequence in a string variable called `sequence`. I also removed carriage returns with the `Str.replace('\n', '')` function and converted all characters to uppercase with the `Str.upper()` function.

The Rockhopper file that was used as input contains information on transcript boundaries (transcription/translation start/stop sites). My script converts this csv file to a pandas dataframe. The column names can be renamed if needed using the `DataFrame.rename(columns = {'OldName': 'NewName'})` function. You can see the type of variable in each column of the dataframe by using the `print(DataFrame.dtypes)` command. In doing so, I saw that some of the Transcription/Translation Start/Stop positions were being taken as float varialbes and I wanted them to be seen as integers, so I changed them with the `DataFrame[columns].applymap(int64)` command.

The results from the GQ search file contains information on the genomic location, the associated gene, the sequence, and the length of each GQ sequence. It was also imported as pandas dataframe.

The first block of script deals with importing all of this data and formatting it properly.

```
#!/usr/local/bin/python3
import sys, fileinput, re
import numpy as np
import pandas as pd

## Saving sequence from the input fasta file as a string object
sequence = ""
file = fileinput.input()
for line in file:
    if line[0] == ">":
        title = line[1:]
    else:
        sequence = sequence + line
sequence = sequence.replace('\n', '').upper()

## Importing data from RNA seq and GQ searching as pandas dataframes, and changing variables to the correct type
RNA_seq = pd.DataFrame.from_csv("/home/gradstd4/RNA_seq_Rockhopper_Results.csv", header = 0, sep = ",", index_col=0)
RNA_seq = pd.DataFrame.dropna(RNA_seq)
#RNA_seq = RNA_seq.rename(columns={"Transcription Start": "Transcription.Start", "Translation Start": "Translation.Start", "Translation Stop": "Translation.Stop"})
#print(RNA_seq.dtypes)
cols = ['Transcription Start', 'Translation Start', 'Translation Stop', 'Transcription Stop']
RNA_seq[cols] = RNA_seq[cols].applymap(np.int64)
#print(RNA_seq)
GQ_seq = pd.DataFrame.from_csv("/home/gradstd4/Streptomyces_venezuelae.csv", header = 0, sep = ",", index_col=0)
cols2 = ['Start']
GQ_seq[cols2] = GQ_seq[cols2].applymap(np.int64)
#print(GQ_seq)
#print(GQ_seq.dtypes)
```

After this, I defined all UTRs from the `RNA_seq` data as anything where the transcription start site and translation start site are separated by 1 or more nucleotides and the translation stop and transcription stop were separated by more than 1 nucleotide. I then saved the start and stop of these UTRs as a list in the variable `UTRs`.

```
## Defining all UTRs in RNA seq data as a list of genomic positions
UTRs = []
for i, row in RNA_seq.iterrows():
    if RNA_seq.loc[i, "Translation Start"] > RNA_seq.loc[i, "Transcription Start"]:
        UTR = [RNA_seq.loc[i, "Transcription Start"], RNA_seq.loc[i, "Translation Start"]]
        UTRs.append(UTR)
    elif RNA_seq.loc[i, "Transcription Start"] > RNA_seq.loc[i, "Translation Start"]:
        UTR = [RNA_seq.loc[i, "Translation Start"], RNA_seq.loc[i, "Transcription Start"]]
        UTRs.append(UTR)
    elif RNA_seq.loc[i, "Transcription Stop"] > RNA_seq.loc[i, "Translation Stop"]:
        UTR = [RNA_seq.loc[i, "Translation Stop"], RNA_seq.loc[i, "Transcription Stop"]]
        UTRs.append(UTR)
    elif RNA_seq.loc[i, "Translation Stop"] > RNA_seq.loc[i, "Transcription Stop"]:
        UTR = [RNA_seq.loc[i, "Transcription Stop"], RNA_seq.loc[i, "Translation Stop"]]
        UTRs.append(UTR)
```

Once I had this list of positions of UTRs, I extracted the sequences from the `sequence` variable and saved them as list in `UTR_seq`. There were many empty strings, so I removed these using the `List.remove('')` function.

```
## Getting the UTR sequences by searhcing the object sequence and extracting the ones between positions listed in the UTR list
UTR_seq = []
for x in UTRs:
    for y in x:
        for z in x:
            UTR_seq.append(sequence[y:z])
while "" in UTR_seq:
    UTR_seq.remove("")
```
Once I had this list of UTR sequences, I searched for GQs in these sequences with the regular expression `GGG[ATGC]{1,7}GGG[ATGC]{1,7}GGG[ATGC]{1,7}GGG` and its complement. The results from the search were saved as list items in the variable `UTR_GQ`. The total number of GQs found in UTRs was defined as the length of this list.

```
## Searching for GQs in the list of UTR sequences
UTR_GQ = []
for j in UTR_seq:
    GQ_for = re.findall("GGG[ATGC]{1,7}GGG[ATGC]{1,7}GGG[ATGC]{1,7}GGG", j)
    if len(GQ_for) >= 1:
        UTR_GQ.append(GQ_for)
    GQ_rev = re.findall("CCC[ATGC]{1,7}CCC[ATGC]{1,7}CCC[ATGC]{1,7}CCC", j)
    if len(GQ_rev) >= 1:
        UTR_GQ.append(GQ_rev)
#print(UTR_GQ)
print('There are', len(UTR_GQ), 'GQs detected  in UTRs')
```

Finally, I wanted to match these hits up to the list of GQs we already have in the pandas dataframe `GQ_seq`. To do this, I first had to convert the `UTR_GQ` object from a list of lists of strings to a list of strings. I did this using the `', '.join()` function in a for loop. I checked the length of the object to make sure it hadn't changed, then used the `DataFrame.isin()` function to find all of the UTR GQ sequences in the list of all GQ sequences. The output was saved to a csv file.

```
## Converts this object from a list of lists to a list of strings
UTR_GQ = [', '.join(x) for x in UTR_GQ]
#print(len(UTR_GQ))

## Extracting information about the GQs found in UTRs by matching them to the GQ search data
final = GQ_seq[GQ_seq['Sequence'].isin(UTR_GQ)]
final.to_csv('/home/gradstd4/UTR_GQ_search_output.csv')
```
An alternate way of doing this analysis is instead of searching for GQs in UTR sequences, UTRs can be defined as a range from the start position of the UTR to the end position of the UTR, and I can use the list of GQ sequences and locations that we already have to find the ones with start sites located within UTRs.

```
## Makes a list of where all GQs start
GQ_starts = []
for j, row in GQ_seq.iterrows():
    GQ_starts.append(GQ_seq.loc[j, 'Start'])

## Defines the function intersect as a function that returns everything found in two input lists
def intersect(a, b):
    return list(set(a) & set(b))

## Finds all GQs that are found in UTRs and creates a list of these and their genomic locations
UTR_GQs = []
for x in UTRs:
    UTR_GQs.append(intersect(x, GQ_starts))

## Removes empty list items
while [] in UTR_GQs:
    UTR_GQs.remove([])

## Changes UTR_GQs from a list of lists of integers to a list of integers
new_UTR_GQs = []
def my_fun(temp_list):
    for ele in temp_list:
        if type(ele) == list:
            my_fun(ele)
        else:
            new_UTR_GQs.append(ele)
my_fun(UTR_GQs)

## Number of GQs located in UTRs
print(len(new_UTR_GQs))

## Matches locations of UTR GQs to list of GQs and pulls out information of these and exports it all to a csv file
final = GQ_seq[GQ_seq['Start'].isin(new_UTR_GQs)]
final.to_csv('/home/gradstd4/UTR_GQ_search_output.csv')
```

