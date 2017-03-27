# Bio709_GQProject

## Introduction

G-quadruplexes (GQs) are structures that form in G-rich nucleic acids. They require four G-tracts of three to six G's separated by short loops of one to seven nucleotides to form. Most of the work to date on the roles of these structures has been done in eukaryotes, but less is known about their functions in bacteria. Streptomyces present a great model system for studying this in bacteria because their extremely GC-rich genomes (>70%) mean they likely have many sequences capable of forming these structures. We have searched the genomes of three model streptomycetes for sequences containing these structures and found that there were approximately 2,500 in each species. In our model streptomycete, Streptomyces venezuelae, most of these were intragenic (~2,300) but there was still a large proportion located in intergenic regions (~600). Our work to date has focused on uncovering the regulatory roles of GQs on gene expression using a series of reporter constructs in which we cloned GQs upstream of a reporter gene to measure the impact of these sequences on reporter activity. We found that depending on the strand in which the GQs were found, they had different effects on reporter activity. This has led us to further investigate the roles of GQs in 3' and 5' UTRs through in silico analyses. 

## Research Questions

1. How frequently do GQs appear in 5' and 3' UTRs? We know that there are ~600 GQs in intergenic regions in S. venezuelae, but we don't know how many of them are in regions that are actually transcribed or where they are located relative to their associated gene (3' or 5').
2. Once we know how many GQs are in transcribed regions, our next question is: do they have any transcriptional effects? 
3. Given the high GC-content of Streptomyces genomes, is this more than the number of GQ sequences we would expect to find by chance in a genome of this size and this GC-content? Additionally, are they enriched in 5' and 3' UTRs compared to the rest of the genome? 

## Methods

### Mapping GQ sequences to 3' and 5' UTRs

We can address this question in three steps:

1. Identfy GQ sequences in the genome and pull out genomic locations of these sequences. 
2. Identify transcript boundaries using Rockhopper.
3. Write a Python script to identify GQ sequences that are located either between transcription and translation starts or between translation and transcription stops.

### Determine transcriptional effects of GQs

This is potentially outside of the scope for this course, but is something I will be interested in looking at in the future. 

### Determine whether GQs are enriched in *Streptomyces*

To answer this question, I will use a permutation approach to determine whether there are more GQ sequences than we would find by chance as follows:

1. Randomly shuffle the *S. venezuelae* genome using uShuffle to preserve doublet/triplet frequencies.
2. Search for number of GQs in randomly shuffled genome, and write this number to an output file.
3. Repeat steps 1-2 1000+ times to get a distribution of number of GQ sequences found. 
4. Determine where the number of actual GQ sequences found falls within this distribution.


