#!/usr/bin/env python

# Import argparse for command line arguments
import argparse

# Import os to handle file name splitting
import os

# Import cairo for drawing
import cairo

# Import regex for reg exp
import re


# Define a function to collect command line arguments
def get_args():

    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Visualize sequence motifs on genomic regions using pycairo")
    # Add FASTA file argument
    parser.add_argument(    "-f",   "--fasta", help="Input fasta containing sequences", required=True)
    # Add motifs file argument
    parser.add_argument("-m", "--motifs",help="Input text file containing motifs (one per line)",required=True)

    # Parse and return the arguments
    return parser.parse_args()


# Create a class to represent a continuous exon or intron segment
class Segment:
    """Represents one exon or intron region"""

    # Constructor: runs when you create Segment(...)
    def __init__(self, start, end, kind):

        # Store the start coordinate
        self.start = start

        # Store the end coordinate (end-exclusive)
        self.end = end

        # Store whether this is exon or intron
        self.kind = kind

    # Function that returns the length of this segment
    def length(self):

        # Return end minus start
        return self.end - self.start


# Create a class that represents one FASTA record
class SequenceRecord:
    """Represents one FASTA record with exon/intron segments"""

    # Constructor: runs when you create SequenceRecord(...)
    def __init__(self, header, sequence):

        # Save header line
        self.header = header

        # Save raw sequence string
        self.seq = sequence

        # Save sequence length
        self.length = len(sequence)

        # Create an empty list of segments for exon/intron blocks
        self.segments = []

        # Build exon/intron segments from uppercase/lowercase pattern
        self._build_segments()

    # Internal helper to split exon/intron regions based on case
    def _build_segments(self):

        # If the first base is uppercase, we start in an exon
        if self.seq[0].isupper():
            current_type = "exon"
        else:
            current_type = "intron"

        # Start of the current segment
        start = 0

        # Loop through the sequence positions starting at index 1
        for i in range(1, len(self.seq)):

            # Decide what type this base is
            if self.seq[i].isupper():
                base_type = "exon"
            else:
                base_type = "intron"

            # If we changed from exon->intron or intron->exon
            if base_type != current_type:

                # Create a Segment for the previous run
                seg = Segment(start, i, current_type)

                # Add the segment to the segments list
                self.segments.append(seg)

                # Start a new segment at position i
                start = i

                # Update current_type to the new type
                current_type = base_type

        # After the loop, add the final segment
        final_seg = Segment(start, len(self.seq), current_type)

        # Store that final segment too
        self.segments.append(final_seg)


# New class: one motif match at one location in one gene
class MotifLocation:
    """Represents one motif match on one gene (start/end) with a color."""

    # Constructor to store motif match information
    def __init__(self, motif, start, end, color):

        # Store the motif string (example: "YGCY")
        self.motif = motif

        # Store match start (0-based)
        self.start = start

        # Store match end (end-exclusive)
        self.end = end

        # Store motif color (RGB tuple)
        self.color = color

    # Convenience function to get length of motif hit
    def length(self):

        # Return end minus start
        return self.end - self.start


# Function to read a FASTA file and return SequenceRecord objects
def read_fasta(path):

    # Make an empty list to store records
    records = []

    # Store current header (None before we see one)
    header = None

    # Store sequence lines for the current record
    seq_lines = []

    # Open the FASTA file
    with open(path, "r") as fh:

        # Loop through each line
        for line in fh:

            # Remove newline and spaces
            line = line.strip()

            # Skip blank lines
            if line == "":
                continue

            # If this is a header line
            if line.startswith(">"):

                # If we already have a header, finalize the previous record
                if header is not None:

                    # Combine stored sequence lines into one string
                    sequence = "".join(seq_lines)

                    # Create a SequenceRecord and append it
                    rec = SequenceRecord(header, sequence)
                    records.append(rec)

                # Update header to this new header line
                header = line

                # Reset sequence lines list
                seq_lines = []

            else:
                # Otherwise this is a sequence line, so store it
                seq_lines.append(line)

        # After loop ends, add the last record if header exists
        if header is not None:

            # Combine last sequence
            sequence = "".join(seq_lines)

            # Create final SequenceRecord
            rec = SequenceRecord(header, sequence)
            records.append(rec)

    # Return the list of SequenceRecord objects
    return records


# Function to read motifs file (one motif per line)
def read_motifs(path):

    # Make empty list for motifs
    motifs = []

    # Open motif file
    with open(path, "r") as fh:

        # Loop through each line
        for line in fh:

            # Strip newline
            line = line.strip()

            # Skip blank lines
            if line == "":
                continue

            # Add motif to list
            motifs.append(line)

    # Return list of motif strings
    return motifs


# Map IUPAC codes to regex character classes
IUPAC_REGEX = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",        # treat U as T
    "Y": "[CT]",     # pyrimidine
}

# Convert a motif like "YGCY" into a regex like "[CT]GC[CT]"
def motif_to_regex(motif):

    # Make motif uppercase
    motif = motif.upper()

    # Start with empty regex pattern
    pattern = ""

    # Build regex one character at a time
    for ch in motif:

        # If the char is known, add its regex
        if ch in IUPAC_REGEX:
            pattern = pattern + IUPAC_REGEX[ch]

        # If unknown, fail loudly (or you can return None)
        else:
            raise ValueError("Unsupported IUPAC code in motif: " + ch)

    # Return finished regex pattern
    return pattern


# Find motif hits using regex lookahead so overlaps are included
def find_motif_locations(record, motif, color):

    # Make sequence uppercase
    seq = record.seq.upper()

    # Convert motif to regex
    inner = motif_to_regex(motif)

    # Use lookahead so overlapping matches are found
    pattern = "(?=(" + inner + "))"

    # Store hits
    hits = []

    # Loop over every regex match
    for match in re.finditer(pattern, seq):

        # Start position is where lookahead begins
        start = match.start()

        # End is start + motif length
        end = start + len(motif)

        # Create MotifLocation object
        hit_obj = MotifLocation(motif.upper(), start, end, color)

        # Store it
        hits.append(hit_obj)

    # Return all hits
    return hits

# Convert FASTA filename to PNG name
def fasta_to_png_name(fasta_path):

    # Get base file name (no directories)
    base = os.path.basename(fasta_path)

    # Split name from extension
    prefix, ext = os.path.splitext(base)

    # Return new name with .png
    return prefix + ".png"


# Draw genes, motifs, and overlap lanes
def draw_gene_bases(records, motifs, out_png):

    # Width of the image
    width = 1400

    # Left margin for labels
    left_margin = 220

    # Right margin for spacing
    right_margin = 40

    # Top margin
    top_margin = 50

    # Bottom margin
    bottom_margin = 40

    # Baseline thickness
    line_width = 2

    # Exon rectangle height
    exon_height = 60

    # Main motif rectangle height
    motif_height = 60

    # Padding above/below features
    pad_y = 10

    # Space reserved for the label region
    label_space = 30

    # Motif transparency
    motif_alpha = 0.90

    # Lane bar height
    lane_height = 6

    # Lane gap (big gap)
    lane_gap = 10

    # Gap between gene and overlap strip
    strip_top_gap = 25

    # Extra strip padding at the bottom
    strip_bottom_pad = 10

    # If records is empty, stop
    if len(records) == 0:
        raise ValueError("No records to draw")

    # Find the maximum gene length
    max_len = 0
    for r in records:
        if r.length > max_len:
            max_len = r.length

    # Compute usable width
    usable_w = width - left_margin - right_margin

    # Compute scale
    scale = usable_w / max_len

    # Convert bp position to x coordinate
    def x(bp):
        return left_margin + bp * scale

    # Compute half exon height
    half_exon = exon_height / 2

    # Compute half motif height
    half_motif = motif_height / 2

    # Compute half feature height
    if half_exon > half_motif:
        half_feature = half_exon
    else:
        half_feature = half_motif

    # Create palette colors
    palette = [
        (0.00, 0.80, 1.00),
        (0.00, 1.00, 0.40),
        (1.00, 0.55, 0.00),
        (1.00, 0.10, 0.70),
    ]

    # Build motif to color mapping
    motif_to_color = {}

    # Loop over motif list
    i = 0
    while i < len(motifs):

        # Get motif
        motif = motifs[i]

        # Uppercase motif
        motif_up = motif.upper()

        # Pick color index
        color_index = i % len(palette)

        # Store mapping
        motif_to_color[motif_up] = palette[color_index]

        # Increment
        i += 1

    # Define helper to assign hits to lanes so overlaps don't share a lane
    def assign_lanes(hits):

        # Sort hits by start then end
        hits_sorted = sorted(hits, key=lambda h: (h.start, h.end))

        # Store lanes
        lanes = []

        # Store last end for each lane
        lane_last_end = []

        # Loop over hits
        for hit in hits_sorted:

            # Track if placed
            placed = False

            # Try each lane
            lane_index = 0
            while lane_index < len(lanes):

                # If it doesn't overlap the last one in this lane
                if hit.start >= lane_last_end[lane_index]:

                    # Add to this lane
                    lanes[lane_index].append(hit)

                    # Update last end
                    lane_last_end[lane_index] = hit.end

                    # Mark as placed
                    placed = True

                    # Stop searching lanes
                    break

                # Otherwise try next lane
                lane_index += 1

            # If not placed, create a new lane
            if placed is False:
                lanes.append([hit])
                lane_last_end.append(hit.end)

        # Return lanes
        return lanes

    # Precompute lanes for each gene
    per_gene_lanes = []

    # Precompute main hits for each gene (so we don't recalc twice)
    per_gene_hits = []

    # Loop over each record
    for rec in records:

        # Store all hits for this gene
        all_hits = []

        # Loop over motifs
        for motif in motifs:

            # Uppercase motif
            motif_up = motif.upper()

            # Get color
            color = motif_to_color[motif_up]

            # Find hits as MotifHit objects
            hits = find_motif_locations(rec, motif_up, color)

            # Append each hit into all_hits
            for h in hits:
                all_hits.append(h)

        # Assign lanes using hits list
        lanes = assign_lanes(all_hits)

        # Store lanes
        per_gene_lanes.append(lanes)

        # Store all hits for main drawing
        per_gene_hits.append(all_hits)

    # Find maximum lanes across genes
    max_lanes = 0
    for lanes in per_gene_lanes:
        if len(lanes) > max_lanes:
            max_lanes = len(lanes)

    # Compute strip height
    strip_height = 0
    if max_lanes > 0:
        strip_height = (max_lanes * lane_height) + ((max_lanes - 1) * lane_gap)

    # Compute row height
    row_height = int(
        label_space +
        (2 * half_feature) +
        (2 * pad_y) +
        strip_top_gap +
        strip_height +
        strip_bottom_pad
    )

    # Compute full image height
    height = top_margin + (len(records) * row_height) + bottom_margin

    # Create cairo surface
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)

    # Create cairo context
    cr = cairo.Context(surface)

    # Set background color white
    cr.set_source_rgb(1, 1, 1)

    # Paint background
    cr.paint()

    # Set font
    cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

    # Set font size
    cr.set_font_size(16)

    # Loop over genes for drawing
    idx = 0
    while idx < len(records):

        # Get this record
        rec = records[idx]

        # Compute y top
        y_top = top_margin + idx * row_height

        # Compute baseline y
        baseline_y = y_top + label_space + half_feature + pad_y

        # Get header
        gene_name = rec.header

        # Remove >
        if gene_name.startswith(">"):
            gene_name = gene_name[1:]

        # Set text color
        cr.set_source_rgb(0, 0, 0)

        # Move to text position
        cr.move_to(20, baseline_y - half_feature - 15)

        # Draw text
        cr.show_text(gene_name)

        # Set line width
        cr.set_line_width(line_width)

        # Set black for baseline
        cr.set_source_rgb(0, 0, 0)

        # Draw baseline
        cr.move_to(x(0), baseline_y)
        cr.line_to(x(rec.length), baseline_y)
        cr.stroke()

        # Draw exons
        for seg in rec.segments:
            if seg.kind != "exon":
                continue

            x0 = x(seg.start)
            x1 = x(seg.end)
            w = x1 - x0
            if w < 1:
                w = 1

            cr.set_source_rgb(0, 0, 0)
            cr.rectangle(x0, baseline_y - exon_height / 2, w, exon_height)
            cr.fill()

        # Draw motif hits on main track ON TOP of exons
        motif_y_main = baseline_y - motif_height / 2

        # Grab precomputed hits
        hits_for_gene = per_gene_hits[idx]

        # Draw each hit object
        for hit in hits_for_gene:

            x0 = x(hit.start)
            x1 = x(hit.end)
            w = x1 - x0
            if w < 1:
                w = 1

            cr.set_source_rgba(hit.color[0], hit.color[1], hit.color[2], motif_alpha)
            cr.rectangle(x0, motif_y_main, w, motif_height)
            cr.fill()

        # Compute overlap strip start
        strip_y0 = baseline_y + half_feature + strip_top_gap

        # Get lanes for this gene
        lanes = per_gene_lanes[idx]

        # Draw lanes
        li = 0
        while li < len(lanes):

            lane = lanes[li]

            lane_y = strip_y0 + li * (lane_height + lane_gap)

            for hit in lane:

                x0 = x(hit.start)
                x1 = x(hit.end)
                w = x1 - x0
                if w < 1:
                    w = 1

                cr.set_source_rgb(hit.color[0], hit.color[1], hit.color[2])
                cr.rectangle(x0, lane_y, w, lane_height)
                cr.fill()

            li += 1

        idx += 1

    # Save PNG
    surface.write_to_png(out_png)


# Main function
def main():

    # Parse arguments
    args = get_args()

    # Read FASTA records
    records = read_fasta(args.fasta)

    # Read motifs
    motifs = read_motifs(args.motifs)

    # Make output filename
    out_png = fasta_to_png_name(args.fasta)

    # Draw figure
    draw_gene_bases(records, motifs, out_png)

    # Print success
    print("Saved:", out_png)


# Run main if script executed directly
if __name__ == "__main__":
    main()
