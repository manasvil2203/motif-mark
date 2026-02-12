#!/usr/bin/env python

import argparse
import os
import cairo


#Creating my varuables using Argparse
def get_args():
    parser = argparse.ArgumentParser(description= "Visualize sequence motifs on genomic regions using pycairo")
    parser.add_argument("-f", "--fasta", help="Input fasta containing sequences", required=True)
    parser.add_argument("-m", "--motifs", help="Input text file containing motifs (one per line)", required=True)    

    return parser.parse_args()  

args = get_args()

#Making my variable global
f: str = args.fasta
m: str = args.motifs

# Making a class called Segment
class Segment:
    """ Represents one exon or intron region"""

    # defining my constructor
    def __init__(self, start: int, end: int, kind: str): 
        # start positions
        self.start = start 
        #end position
        self.end = end 
        # is it a exon or intron
        self.kind = kind 
    # returns how long the segment is
    def length(self) -> int:
        return self.end - self.start 


class SequenceRecord:
    """ Represents one FATSA record with exon/intron segments"""
    
     #contructor
    def __init__(self, header: str, sequence: str):
        # stores the ehader of the entry
        self.header = header
        # stores the sequence of the entry
        self.seq = sequence
        # stores length once
        self.length = len(sequence) 

        # List of Segment objects
        self.segments = []

        # Build exon/intron blocks immediately
        self._build_segments()
    
    #function to split exons and introna
    def _build_segments(self):
        """Split sequence into exon/intron segments based on case. Uppercase = exon, lowercase = intron."""
        # if it is uppercase
        if self.seq[0].isupper():
            current_type = "exon"
        # if it is lowercase    
        else:
            current_type = "intron"

        # initialize to zero
        start = 0

        # Scan sequence
        for i in range(1, len(self.seq)):
            #if it is uppercase
            if self.seq[i].isupper():
                base_type = "exon"
            # if it is lowercase    
            else:
                base_type = "intron"

            # If type changes, close old segment
            if base_type != current_type:
                # then assign record to class segment, store exon, intron info
                seg = Segment(start, i, current_type)
                # Append continously to the list
                self.segments.append(seg)
                # Assign next record 
                start = i
                # change current type
                current_type = base_type

        # Add final segment
        final_seg = Segment(start, len(self.seq), current_type)
        self.segments.append(final_seg)


def read_fasta(path: str) -> list[SequenceRecord]:
    """Read a FASTA file and return a list of SequenceRecord objects."""
    # make empty list to store each record
    records = []
    # Header has no value, so assign it None for now
    header = None
    # Empty list to stpre sequence lines
    seq_lines = []
    
    # open file
    with open(path, "r") as fh:

        for line in fh:
            # basically gets rid of new line character  
            line = line.strip()
            # If the line is a header
            if line.startswith(">"):
                # If it is not the first one
                if header is not None:
                    #
                    sequence = "".join(seq_lines)
                    # Also initialize class SequenceRecord within the record list and add each record in that class format
                    records.append(SequenceRecord(header, sequence))
                #assign header the next line
                header = line
                # Empty the list
                seq_lines = []
            # if the line is a sequence
            else:
                #add the sequnce to the list of seq lines
                seq_lines.append(line)

        # Add last record
        if header is not None:
            sequence = "".join(seq_lines)
            records.append(SequenceRecord(header, sequence))

    return records

# function to read the list of motifs file
def read_motifs(path: str) -> list[str]:
    """Read a motifs text file (one motif per line) and return a list of motifs."""
    # empty list to store the motifs
    motifs = []
    
    # Open the file
    with open(path, "r") as fh:
        # go through every line
        for line in fh:
            #get rid of new line characters
            line = line.strip()
            # append each motif to our list
            motifs.append(line)

    return motifs


# dictionary of IUPAC ambiguous letters
IUPAC = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},  # treat U as T
    "Y": {"C", "T"},
}

# function to match motif charcters in FASTA files to the ambigous IUPAC charcters
def chars_compatible(m_char: str, s_char: str) -> bool:
    """
    Return True if motif character and sequence character are compatible
    under IUPAC ambiguity rules.
    """
    # Make them uniform by turning both into uppcase
    m_char = m_char.upper()
    s_char = s_char.upper()
    
   # If they are not present in my dictionbary return false
    if m_char not in IUPAC or s_char not in IUPAC:
        return False
    
    # else return the intersection
    return len(IUPAC[m_char].intersection(IUPAC[s_char])) > 0


# function for motif matching!!
def find_motif_positions(record: SequenceRecord, motif: str) -> list[tuple[int, int]]:
    """
    Find all positions where a motif matches a sequence.
    Returns a list of (start, end) tuples (0-based).
    """
    # Make them upper case to keep it uniform
    seq = record.seq.upper()
    motif = motif.upper()
    
    # Initialize empty list
    positions = []
    # variable for the length of the motif
    k = len(motif)

    # slide motif along the sequence
    for i in range(0, len(seq) - k + 1):

        match = True

        # compare motif to sequence at this position
        for j in range(k):
            if not chars_compatible(motif[j], seq[i + j]):
                match = False
                break

        if match:
            positions.append((i, i + k))
        
    return positions


