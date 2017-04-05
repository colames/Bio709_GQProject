# Are GQs enriched in the *S. venezuelae* genome?

## uShuffle

uShuffle is a program for shuffling biological data that allows you to conserve nucleotide frequencies by only shuffling the genome in blocks of length *k*. A wrapper module for using uShuffle in Python can be written by following the instructions at the website below.

http://digital.cs.usu.edu/~mjiang/ushuffle/python.html

For building the module, instead of their commands enter:

```
gcc -c -I/usr/include/python2.6 -fPIC ushuffle.c ushufflemodule.c
gcc -shared ushuffle.o ushufflemodule.o -o ushuffle.so
```
## A Python program that searches fasta files for GQs

Using fasta files as input, the following code will allow you to search for GQ sequences. 

```
#!/usr/bin/python
import sys, fileinput, re
sequence = ""
file = fileinput.input()

for line in file:
    if line[0] == ">":
        title = line[1:]
    else:
        sequence = sequence + line
sequence = sequence.upper().replace("\n", "")

GQ_for = re.findall("GGG[ATGCN]{1,7}GGG[ATGCN]{1,7}GGG[ATGCN]{1,7}GGG", sequence)
GQ_rev = re.findall("CCC[ATGCN]{1,7}CCC[ATGCN]{1,7}CCC[ATGCN]{1,7}CCC", sequence)
GQ = GQ_for + GQ_rev

print(len(GQ))
```
The first block of code separates the title from the DNA in the fasta file and stores the DNA sequence as a string object called "sequence". It then transforms this object into all uppercase characters and removes N's and carriage returns. 

Note, for this to work the file must be in unicode. To change the file from dos to unicode, run the following command:

```
dos2unix <filename>
```

This will save over the existing file in the right format. 

Next, the program uses the Python re.findall function to find all occurences of the regular expression of interest. The re.findall function searches a string from left to right and saves all non-overlapping matches to a list object in the order that they appear in the sequence. In this case, the pattern of interest is a GQ sequence. The Sequence file is searched twice: once for sequences containing G's and once for sequences containing C's (for ones that are found on the opposite strand of DNA).

The two lists are then combined into one list "GQ", which contains all GQ sequences found in the entire fasta file. The program then prints out the total number of GQ's found by finidng the length of this list object, which corresponds to the total number of GQ sequences. 

## Automating this to shuffle the genome repeatedly and search for GQs each time

```
#!/usr/bin/python
import sys, fileinput, ushuffle, re
sequence = ""
num_GQ= []

#separating sequence from title in input file
for line in fileinput.input():
    if line[0] == ">":
        title = line[1:]
    else:
        sequence = sequence + line

#chnage all characters to uppercase, remove N's and carriage returns
sequence = sequence.upper().replace("\n", "").replace("N", "")

#create an object "shuff" that contains the shuffled genome
shuff = ushuffle.shuffle(sequence, len(sequence), 6)

#re-shuffle the shuffled sequence and search for GQs each time
count = 1
while count <= 5:
    shuff = ushuffle.shuffle(shuff, len(shuff), 6)
    GQ_for = re.findall("GGG[ATCG]{1,7}GGG[ATGC]{1,7}GGG[ATGC]{1,7}GGG", shuff)
    GQ_rev = re.findall("CCC[ATCG]{1,7}CCC[ATGC]{1,7}CCC[ATGC]{1,7}CCC", shuff)
    GQ = GQ_for + GQ_rev
    print(count, len(GQ))
#counts the number of GQs found and saves this number of the list "num_GQ"
    num_GQ.append(len(GQ))
    count = count + 1
print(num_GQ)

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