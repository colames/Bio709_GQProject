# Biology 709 final project

## Introduction

G-quadruplexes (GQs) are nucleic acid secondary structures that form in G-rich sequences. To form, they require 4 G-tracts of 3-6 G's separated by no more than 7 nucleotides each. The goal of my research project is to uncover the regulatory roles of these structures in *Streptomyces*. We are interested in studying GQs in *Streptomyces* because they are extremely GC-rich organisms (>70%), and they likely have many GQ-forming sequences in their genomes, but nothing is known about how these might affect the expression of the associated gene. We have already searched the genomes of several *Streptomcyes* species and found ~2,500 of them in each species. For this project, I have three main research questions:
1. Is the actual number of GQ sequences found greater than what we would expect to find by chance given the G-C content of *Streptomyces* genomes?
2. How many of the GQ sequences that we've identified are found in UTRs?
3. Is the number of GQ sequences in UTRs greater than the number we would expect to find by chance in these regions?

## Methods/results

### Question 1: Are GQs enriched in *Streptomyces* genomes?

To address this question, I wrote a Python script (icnluded below) that shuffles an input fasta sequence using the program [uShuffle](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-192), and searches for GQ sequences in the shuffled genome. uShuffle is a program for shuffling biological data that allows you to conserve nucleotide frequencies by only shuffling the genome in blocks of length *k*. A wrapper module for using uShuffle in Python can be written by following the instructions found [here](http://digital.cs.usu.edu/~mjiang/ushuffle/python.html).

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

#chnage all characters to uppercase, remove N's and carriage returns
sequence = sequence.upper().replace("\n", "")

#create an object "shuff" that contains the shuffled genome
shuff = ushuffle.shuffle(sequence, len(sequence), 6)

#re-shuffle the shuffled sequence and search for GQs each time
count = 1
while count <= n:
    shuff = ushuffle.shuffle(shuff, len(shuff), 6)
    GQ_for = re.findall("GGG[ATCGN]{1,7}GGG[ATGCN]{1,7}GGG[ATGCN]{1,7}GGG", shuff)
    GQ_rev = re.findall("CCC[ATCGN]{1,7}CCC[ATGCN]{1,7}CCC[ATGCN]{1,7}CCC", shuff)
    GQ = GQ_for + GQ_rev
    print(count, len(GQ))
#counts the number of GQs found and adds this number to the list "num_GQ"
    num_GQ.append(len(GQ))
    count = count + 1
print(num_GQ)

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
plt.savefig("/home/gradstd4/plot.png")

print("done")
```

When I ran this program on the S. venezuelae genome with n = 1,000, I got an average of approximately 1,250 GQ sequences (Figure 1) while the actual number in this genome is 2,999. Since the actual number is much higher than the average, we can say that the is some enrichment of GQs in this genome.

![Figure 1](C:/Users/Owner/Desktop/plot(1).png)