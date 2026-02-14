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


# Now here comes the drawing part
import cairo


def draw_gene_bases(records: list[SequenceRecord], motifs: list[str], out_png: str) -> None:
    """
    Draw genes with introns, exons, and motifs in one PNG.
    Spacing is computed automatically per gene.
    Motifs are semi-transparent so overlaps are visible.
    """

    # Set total image width in pixels
    width = 1400

    # Space on the left for gene labels
    left_margin = 220

    # Space on the right so lines don't touch edge
    right_margin = 40

    # Padding at the top of the image
    top_margin = 50

    # Padding at the bottom of the image
    bottom_margin = 40

    # Thickness of baseline line
    line_width = 2

    # Height of exon rectangles
    exon_height = 60

    # Height of motif rectangles
    motif_height = 60

    # Vertical padding above and below each gene
    pad_y = 10

    # Space reserved for gene name text
    label_space = 30

    # Transparency for motifs (0 = invisible, 1 = solid)
    motif_alpha = 0.60

    # Stop if no records were given
    if not records:
        raise ValueError("No records to draw")

    # Find longest gene in bases
    max_len = max(r.length for r in records)

    # Compute horizontal drawing space
    usable_w = width - left_margin - right_margin

    # Convert base pairs to pixels
    scale = usable_w / max_len

    # Convert base position to x coordinate
    def x(bp: int) -> float:
        return left_margin + bp * scale

    # Compute half exon height
    half_exon = exon_height / 2

    # Compute half motif height
    half_motif = motif_height / 2

    # Find the tallest feature
    half_feature = max(half_exon, half_motif)

    # Compute height needed for one gene row
    row_height = int(label_space + (2 * half_feature) + (2 * pad_y))

    # Compute full image height
    height = top_margin + len(records) * row_height + bottom_margin

    # Create blank PNG surface
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)

    # Create drawing context
    cr = cairo.Context(surface)

    # Set background color to white
    cr.set_source_rgb(1, 1, 1)

    # Paint entire background
    cr.paint()

    # Choose font family
    cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

    # Set font size
    cr.set_font_size(16)

    # Define colors for motifs
    palette = [
        (0.00, 0.80, 1.00),  # Neon Blue
        (0.00, 1.00, 0.40),  # Neon Green
        (1.00, 0.55, 0.00),  # Neon Orange
        (1.00, 0.10, 0.70),  # Neon Pink
    ]

    # Create dictionary mapping motif to color
    motif_to_color = {}

    # Assign each motif a color
    for i, motif in enumerate(motifs):
        motif_to_color[motif.upper()] = palette[i % len(palette)]

    # Loop over each gene
    for idx, rec in enumerate(records):

        # Compute top y position for this gene
        y_top = top_margin + idx * row_height

        # Compute baseline y position
        baseline_y = y_top + label_space + half_feature + pad_y

        # Get gene header
        gene_name = rec.header

        # Remove leading > if present
        if gene_name.startswith(">"):
            gene_name = gene_name[1:]

        # Set color for text
        cr.set_source_rgb(0, 0, 0)

        # Move pen to label position
        cr.move_to(20, baseline_y - half_feature - 15)

        # Draw gene name
        cr.show_text(gene_name)

        # Set line thickness
        cr.set_line_width(line_width)

        # Set line color to black
        cr.set_source_rgb(0, 0, 0)

        # Move pen to start of baseline
        cr.move_to(x(0), baseline_y)

        # Draw baseline to end of gene
        cr.line_to(x(rec.length), baseline_y)

        # Render baseline
        cr.stroke()

        # Loop through exon/intron segments
        for seg in rec.segments:

            # Skip introns
            if seg.kind != "exon":
                continue

            # Convert exon start to x pixel
            x0 = x(seg.start)

            # Convert exon end to x pixel
            x1 = x(seg.end)

            # Compute exon width
            w = x1 - x0

            # Ensure minimum width
            if w < 1:
                w = 1

            # Set color to black
            cr.set_source_rgb(0, 0, 0)

            # Draw exon rectangle
            cr.rectangle(x0, baseline_y - exon_height / 2, w, exon_height)

            # Fill exon box
            cr.fill()

        # Compute vertical position of motifs
        motif_y = baseline_y - motif_height / 2

        # Loop through each motif
        for motif in motifs:

            # Convert motif to uppercase
            motif_up = motif.upper()

            # Get motif color
            color = motif_to_color[motif_up]

            # Find all matches in this gene
            hits = find_motif_positions(rec, motif_up)

            # Draw each motif hit
            for start, end in hits:

                # Convert start position to x pixel
                x0 = x(start)

                # Convert end position to x pixel
                x1 = x(end)

                # Compute width in pixels
                w = x1 - x0

                # Ensure minimum width
                if w < 1:
                    w = 1

                # Set color with transparency
                cr.set_source_rgba(color[0], color[1], color[2], motif_alpha)

                # Draw motif rectangle
                cr.rectangle(x0, motif_y, w, motif_height)

                # Fill motif box
                cr.fill()

    # Save final image to file
    surface.write_to_png(out_png)



def fasta_to_png_name(fasta_path: str) -> str:
    """
    Convert FASTA filename to PNG filename.
    Example: Figure_1.fasta -> Figure_1.png
    """
    base = os.path.basename(fasta_path)
    prefix, _ = os.path.splitext(base)
    return f"{prefix}.png"


def main():
    args = get_args()
    records = read_fasta(args.fasta)
    motifs = read_motifs(args.motifs)

    out_png = fasta_to_png_name(args.fasta)
    draw_gene_bases(records, motifs, out_png)

    print("Saved:", out_png)

if __name__ == "__main__":
    main()
