## Motif-Mark

### What is this?
This program takes DNA sequences and motif patterns and makes a picture showing where the motifs appear on each gene.It was made for a class assignment to practice Python, OOP, and visualization.

### What does it show?
For each gene, the output image shows:

- A black line for the gene

- Black boxes for exons

- Colored boxes for motifs

- Thin colored lines below for overlapping motifs

- A legend at the bottom

### What files do I need?

**1. FASTA file**

Contains gene sequences.

**Rules:**

- > starts a new gene

- Uppercase = exon

- Lowercase = intron

**Example:**
```
>Gene1
ATGCttaaATGC
```

**2. Motif file**

Text file with one motif per line.

**Example:**
```
YGCY
ATG
CGT
```


### How to run it

In the terminal:

```
python motif-mark.py -f genes.fasta -m motifs.txt
```

This creates:

```
genes.png
```
### How motif matching works

Motifs use IUPAC codes.

**Example:**
```
Y = C or T
```
The program turns motifs into regex patterns and searches the DNA.

It also finds overlapping matches.

### How overlaps are handled

If motifs overlap, they would cover each other.

**To fix this:**

- The program puts overlapping motifs into different “lanes”

- Each lane is drawn on a separate line

- This makes overlaps easy to see

### Tools used

- Python 3

- PyCairo

- Regular expressions