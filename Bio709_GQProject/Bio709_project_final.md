# Biology 709 final project

## Introduction

G-quadruplexes (GQs) are nucleic acid secondary structures that form in G-rich sequences. To form, they require 4 G-tracts of 3-6 G's separated by no more than 7 nucleotides each. The goal of my research project is to uncover the regulatory roles of these structures in *Streptomyces*. We are interested in studying GQs in *Streptomyces* becasue their extremely GC-rich genomes (>70%) means that they likely have many sequences with the potential to form GQs, yet nothing is known about how these structures might affect gene expression in these organisms. We have already searched the genomes of several *Streptomcyes* species and found ~2,500 of them in each species. For this project, I have three main research questions:
1. Is the actual number of GQ sequences found greater than what we would expect to find by chance given the G-C content of *Streptomyces* genomes?
2. How many of the GQ sequences that we've identified are found in UTRs?
3. Is the number of GQ sequences in UTRs greater than the number we would expect to find by chance in these regions?

## Methods/results

### Question 1: Are GQs enriched in *Streptomyces* genomes?

To address this question, I wrote a Python script (icnluded below) that shuffles an input fasta sequence using the program [uShuffle](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-192), and searches for GQ sequences in the shuffled genome. uShuffle is a program for shuffling biological data that allows you to conserve doublet/triplet frequencies by only shuffling the genome in blocks of length *k*. A wrapper module for using uShuffle in Python can be written by following the instructions found [here](http://digital.cs.usu.edu/~mjiang/ushuffle/python.html).

For building the module, instead of their commands enter:

```
gcc -c -I/usr/include/python2.6 -fPIC ushuffle.c ushufflemodule.c
gcc -shared ushuffle.o ushufflemodule.o -o ushuffle.so
```

#### My script

My script uses the `re.findall()` regular expression search function to find all non-overlapping GQ sequences in the shuffled genome. This fuction comes up with a list of all occurrences of the pattern, so the number of GQ sequences found in the shuffled genome is simply `len(re.findall(pattern))`. The number of GQ sequences found in the shuffled genome is recorded in the list `num_GQ`, and the genome is re-shuffled. This can be repeated as many times as you want by simply changing n in the line: `while count <= n` to the number of times you want to resample. In the end, the program can plot a histogram, or export the data to a csv file for further statiscal analysis. The result is a distribution of the number of GQ sequences found in the shuffled genomes. 

Note, for this to work the file must be in unicode. To change the file from dos to unicode, run the following command:

```
dos2unix <filename>
```

This will save over the existing file in the right format.

```
#!/usr/bin/python
import sys, fileinput, ushuffle, re, csv
sequence = ""
num_GQ= []

#separating sequence from title in input file
for line in fileinput.input():
    if line[0] == ">":
        title = line[1:]
    else:
        sequence = sequence + line

#chnage all characters to uppercase and remove carriage returns
sequence = sequence.upper().replace("\n", "")

#create an object "shuff" that contains the shuffled genome. I used k-let size of 6 to try to conserve codon usage to some extent.
shuff = ushuffle.shuffle(sequence, len(sequence), 6)

#re-shuffle the shuffled sequence and search for GQs each time
count = 1
while count <= n:
    shuff = ushuffle.shuffle(shuff, len(shuff), 6)
    GQ_for = re.findall("GGG[ATCGN]{1,7}GGG[ATGCN]{1,7}GGG[ATGCN]{1,7}GGG", shuff)
    GQ_rev = re.findall("CCC[ATCGN]{1,7}CCC[ATGCN]{1,7}CCC[ATGCN]{1,7}CCC", shuff)
    GQ = GQ_for + GQ_rev
    print(len(GQ))
#counts the number of GQs found and adds this number to the list "num_GQ"
    num_GQ.append(len(GQ))
    count = count + 1

#writes data to csv file
with open('shuffled_GQ_search.csv', 'w') as output:
    writer = csv.writer(output, lineterminator = '\n')
    for val in num_GQ:
        writer.writerow([val])

#plotting histogram
import matplotlib

#required to force matplotlib to not use any Xwindows backend
matplotlib.use("Agg")

import matplotlib.pyplot as plt

plt.hist(num_GQ)

plt.xlabel('number of GQs')
plt.ylabel('count')
plt.title('Histogram of GQs')
plt.grid(True)
plt.savefig(<filename>)

print("done")
```

When I ran this program on the *S. venezuelae* genome with n = 1,000, I got an average of approximately 1,250 GQ sequences (Figure 1) while the actual number in this genome is 2,999. Since the actual number is much higher than the average, we can say that the is some enrichment of GQs in this genome.

![Figure 1](https://cloud.githubusercontent.com/assets/26418440/25544579/117d36de-2c29-11e7-88eb-966155aa0c98.png)

**Figure 1:** Distribution of nubmer of GQs found in shuffled *S. venezuelae* genome (n = 1,000).

### Question 2: How many GQ sequences are in UTRs?

We know from the GQ search data that we already have that there are ~600 GQ sequences in intergenic regions. However, we don't konw how many of these are in the untrasnlated regions (UTRs) of transcribed genes. We have *in vivo* reporter assays that demonstrate that GQ sequences in these regions can influence gene expression, and we would like to find genes that have GQ sequences in the UTRs so that we can follow up with them experimentally to understand how these sequences might be affecting the expression of their associated genes. 

To do this, I first ran RNA-seq data from our lab through Rockhopper, which gave me information on predicted transcription/translation start/stop sites were for all genes that were expressed under our experimental conditions. I then wrote a Python script (included below) that took the Rockhopper output (a csv file) as well as the output from our previous GQ search (also a csv file) and converts them to Pandas dataframes. Once the data was in this format, I was able to define UTRs as ranges from transcription starts to translation starts and translation stops to transcription stops. Once I had defined UTRs in this way, I could find all GQs that start in UTRs by using the intersect function to return values that were found in both the list of UTRs and the list of GQ locations. I could then search for these values in the full dataset to return information on these GQ sequences. 

```
#!/usr/local/bin/python3
import sys, fileinput, re
import numpy as np
import pandas as pd

## Importing data from csv files as Pandas dataframes and changing required columns to integer values
#print(RNA_seq)
GQ_seq = pd.DataFrame.from_csv("/home/gradstd4/VenezuelaeSVEN.csv", header = 0, sep = ",", index_col=0)
cols2 = ['Start']
GQ_seq[cols2] = GQ_seq[cols2].applymap(np.int64)
#print(GQ_seq)
#print(GQ_seq.dtypes)

## Defining all UTRs in RNA seq data as ranges from genomic start to stop positions
UTRs = []
for i, row in RNA_seq.iterrows():
    if RNA_seq.loc[i, "Translation Start"] > RNA_seq.loc[i, "Transcription Start"]:
        UTR = range(RNA_seq.loc[i, "Transcription Start"], RNA_seq.loc[i, "Translation Start"])
        UTRs.append(UTR)
    elif RNA_seq.loc[i, "Transcription Start"] > RNA_seq.loc[i, "Translation Start"]:
        UTR = range(RNA_seq.loc[i, "Translation Start"], RNA_seq.loc[i, "Transcription Start"])
        UTRs.append(UTR)
    elif RNA_seq.loc[i, "Transcription Stop"] > RNA_seq.loc[i, "Translation Stop"]:
        UTR = range(RNA_seq.loc[i, "Translation Stop"], RNA_seq.loc[i, "Transcription Stop"])
        UTRs.append(UTR)
    elif RNA_seq.loc[i, "Translation Stop"] > RNA_seq.loc[i, "Transcription Stop"]:
        UTR = range(RNA_seq.loc[i, "Transcription Stop"], RNA_seq.loc[i, "Translation Stop"])
        UTRs.append(UTR)

## Creates a list of all GQ start positions form the Pandas dataframe
GQ_starts = []
for j, row in GQ_seq.iterrows():
    GQ_starts.append(GQ_seq.loc[j, 'Start'])

## Defines intersect as a function that takes two lists and returns values that are found in both lists
def intersect(a, b):
    return list(set(a) & set(b))

## Finds all GQs that start in UTRs
for x in UTRs:
    UTR_GQs.append(intersect(x, GQ_starts))

## Removes empty values
while [] in UTR_GQs:
    UTR_GQs.remove([])

## Flattens UTR_GQs and saves it as new_UTR_GQs so that it is now a list of integers instead of a list of lists of integers
new_UTR_GQs = []
def my_fun(temp_list):
    for ele in temp_list:
        if type(ele) == list:
            my_fun(ele)
        else:
            new_UTR_GQs.append(ele)
my_fun(UTR_GQs)

## Prints number of GQs that are found in UTRs
print(len(new_UTR_GQs))

## Takes values from the list of GQs in UTRs and finds them in the Pandas dataframe, then takes that row from the dataframe and saves it to a new file
final = GQ_seq[GQ_seq['Start'].isin(new_UTR_GQs)]
final.to_csv('/home/gradstd4/UTR_GQ_search_output.csv')
```

When I did this, I found 98 GQs in UTRs, which is a manageable list for me to look at manually using RNA-seq data to determine which ones to follow up on experimentally. 

### Question 3: Are GQs enriched in UTRs?

Now that I know that there are 98 GQs in UTRs, I wondered whether this represented an enrichment in these regions. To address this question, I combined my scripts from the first two questions (above) to shuffle the genome and search for GQs in the shuffled genome, but only in the regions that are defined as UTRs in the actual genome.

I started by using the same method as before to define all UTRs, only this time I saved this list as a Python object so that I could import it into another script. 

```
#!/usr/local/bin/python3
import sys, fileinput, re
import numpy as np
import pandas as pd
import csv

## Importing data from RNA seq and GQ searching as pandas dataframes, and changing variables to the correct type
RNA_seq = pd.DataFrame.from_csv("/home/gradstd4/RNA_seq_Rockhopper_Results.csv", header = 0, sep = ",", index_col=0)
RNA_seq = pd.DataFrame.dropna(RNA_seq)
cols = ['Transcription Start', 'Translation Start', 'Translation Stop', 'Transcription Stop']
RNA_seq[cols] = RNA_seq[cols].applymap(np.int64)

## Defining all UTRs in RNA seq data as a list of genomic positions
UTRs = []
for i, row in RNA_seq.iterrows():
    if RNA_seq.loc[i, "Translation Start"] > RNA_seq.loc[i, "Transcription Start"]:
        UTR = list(range(RNA_seq.loc[i, "Transcription Start"], RNA_seq.loc[i, "Translation Start"]))
        UTRs.append(UTR)
    elif RNA_seq.loc[i, "Transcription Start"] > RNA_seq.loc[i, "Translation Start"]:
        UTR = list(range(RNA_seq.loc[i, "Translation Start"], RNA_seq.loc[i, "Transcription Start"]))
        UTRs.append(UTR)
    elif RNA_seq.loc[i, "Transcription Stop"] > RNA_seq.loc[i, "Translation Stop"]:
        UTR = list(range(RNA_seq.loc[i, "Translation Stop"], RNA_seq.loc[i, "Transcription Stop"]))
        UTRs.append(UTR)
    elif RNA_seq.loc[i, "Translation Stop"] > RNA_seq.loc[i, "Transcription Stop"]:
        UTR = list(range(RNA_seq.loc[i, "Transcription Stop"], RNA_seq.loc[i, "Translation Stop"]))
        UTRs.append(UTR)
flat_UTRs = []
for x in UTRs:
    for y in x:
        flat_UTRs.append(y)
#print(flat_UTRs)

with open('UTRs.py', 'w') as f:
    f.write('UTRs = %s' % UTRs)
```

I then wrote a separate script that shuffles the genome and looks for GQs using the `re.finditer()` function, which searches for all non-overlapping matches and returns their positions within the string. Using these positions and the UTR positions (imported from Python object created above), I was able to determine the number of GQs in these same regions in the shuffled genome. The number of re-samplings can be easily varried in the same way as in my code from question 1.

```
#!/usr/bin/python

import sys, fileinput, re, ushuffle, csv
import numpy as np

sequence = ""
file = fileinput.input()

UTR_list = []
from UTRs import UTRs as UTR_list
#print(UTR_list)

## Removing title from fasta genome file, transforming to uppercase characters, and removing carriage returns
for line in file:
    if line[0] == ">":
        title = line[1:]
    else:
        sequence = sequence + line
sequence = sequence.upper().replace("\n", "")

shuff = ushuffle.shuffle(sequence, len(sequence), 6)

## Defines intersect as a function that takes two lists and returns items that are found in both lists
def intersect(a, b):
    return list(set(a) & set(b))

p1 = re.compile("GGG[ATGCN]{1,7}GGG[ATGCN]{1,7}GGG[ATGCN]{1,7}GGG")
p2 = re.compile("CCC[ATGCN]{1,7}CCC[ATGCN]{1,7}CCC[ATGCN]{1,7}CCC")

## Flattens the list from a list of lists of integers to a list of integers
def my_fun(temp_list):
    for ele in temp_list:
        if type(ele) == list:
            my_fun(ele)
        else:
            new_UTR_GQs.append(ele)

## Shuffling the genome, then searching for GQs and determining how many are in UTRs by comparing their start positions to the list of UTR positions
num_GQ = []
count = 1
while count <= n:
    shuff = ushuffle.shuffle(shuff, len(shuff), 6)
    GQs = []
    UTR_GQs = []
    for m in p1.finditer(shuff):
        GQs.append(int(m.start()))
    for m in p2.finditer(shuff):
        GQs.append(int(m.start()))
    for x in UTR_list:
        UTR_GQs.append(intersect(x, GQs))
    while [] in UTR_GQs:
        UTR_GQs.remove([])
    new_UTR_GQs = []
    my_fun(UTR_GQs)
    print(len(new_UTR_GQs))
    num_GQ.append(len(new_UTR_GQs))
    count = count + 1

#print(num_GQ)

## Saves number of GQs found in UTRs as a csv file
with open('shuffled_UTR_search.csv', 'w') as output:
    writer = csv.writer(output, lineterminator = '\n')
    for val in num_GQ:
        writer.writerow([val])

#plotting histogram of the nubmer of GQs found
import matplotlib

#required to force matplotlib to not use any Xwindows backend
matplotlib.use("Agg")

import matplotlib.pyplot as plt

plt.hist(num_GQ)

plt.xlabel('number of GQs')
plt.ylabel('count')
plt.title('Histogram of GQs')
plt.grid(True)
plt.savefig(<filename>)
```

When I chnaged the number of genome shuffles to 1,000, I got the following distribution with a mean around 20 GQs (Figure 2). This is much lower than the actual number of GQs found in UTRs (98), indicating that there is some enrichment of GQs in UTRs.

![Figure 2](https://cloud.githubusercontent.com/assets/26418440/25546064/6375095c-2c2f-11e7-8af6-10cdad923315.png)

**Figure 2:** Distribution of the number of GQs found in regions that were defined as UTRs in the actual *S. venezuelae* genome (n = 1,000).

## Conclusions and future directions

The findings from this project were important because they helped validate *Streptomcyes* as a model organism for studying GQs since they appear to be enriched in these species despite having such high G-C content. We also now have a list of GQ sequences in UTRs that we will be able to follow up on experimentally in oder to establish a role for these strucutres in gene regulation. We also showed that there are more of these sequences in UTRs than we did with our random shuffling approach, which means that many of the ones that are found there likely serve some regulatory funciton. Moving forward, I would like to broaden our study of GQ sequences to include ones that are found within coding regions. This is where most of them were found (~80%), so I am interested in how these sequences might affect gene expression at the translational level. 
