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

    def __init__(self, start: int, end: int, kind: str): # defining my constructor
        self.start = start # strta positions
        self.end = end #end position
        self.kind = kind # is it a exon or intron

    def length(self) -> int:
        return self.end - self.start # returns how long the segment is


class SequenceRecord:
    """ Represents one FATSA record with exon/intron segments"""

    def __init__(self, header: str, sequence: str): # contrcutor
        self.header = header
        self.seq = sequence
        self.length = len(sequence) # stores length once

        # List of Segment objects
        self.segments = []

        # Build exon/intron blocks immediately
        self._build_segments()

    def _build_segments(self):
        """Split sequence into exon/intron segments based on case. Uppercase = exon, lowercase = intron."""

        if self.seq[0].isupper():
            current_type = "exon"
        else:
            current_type = "intron"

        
        start = 0

        # Scan sequence
        for i in range(1, len(self.seq)):

            if self.seq[i].isupper():
                base_type = "exon"
            else:
                base_type = "intron"

            # If type changes, close old segment
            if base_type != current_type:

                seg = Segment(start, i, current_type)
                self.segments.append(seg)

                start = i
                current_type = base_type

        # Add final segment
        final_seg = Segment(start, len(self.seq), current_type)
        self.segments.append(final_seg)


def read_fasta(path: str) -> list[SequenceRecord]:
    """Read a FASTA file and return a list of SequenceRecord objects."""
    records = []
    header = None
    seq_lines = []

    with open(path, "r") as fh:

        for line in fh:

            line = line.strip()

            if line.startswith(">"):

                if header is not None:
                    sequence = "".join(seq_lines)
                    records.append(SequenceRecord(header, sequence))

                header = line
                seq_lines = []

            else:
                seq_lines.append(line)

        # Add last record
        if header is not None:
            sequence = "".join(seq_lines)
            records.append(SequenceRecord(header, sequence))

    return records

def read_motifs(path: str) -> list[str]:
    """Read a motifs text file (one motif per line) and return a list of motifs."""
    motifs = []

    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()

            motifs.append(line)

    return motifs



IUPAC = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},  # treat U as T
    "Y": {"C", "T"},
}

def chars_compatible(m_char: str, s_char: str) -> bool:
    """
    Return True if motif character and sequence character are compatible
    under IUPAC ambiguity rules.
    """
    m_char = m_char.upper()
    s_char = s_char.upper()

    if m_char not in IUPAC or s_char not in IUPAC:
        return False

    return len(IUPAC[m_char].intersection(IUPAC[s_char])) > 0


def find_motif_positions(record: SequenceRecord, motif: str) -> list[tuple[int, int]]:
    """
    Find all positions where a motif matches a sequence.
    Returns a list of (start, end) tuples (0-based, end-exclusive).
    """
    seq = record.seq.upper()
    motif = motif.upper()

    positions = []
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


