""" Exon / transcript coordinate maths shared by the GTF/GFF3 parsers and the UTA converter

    Exons are tuples of (alt_start, alt_end, exon_id, cdna_start, cdna_end, gap): genomic coordinates
    0-based half-open, transcript (cdna) coordinates 1-based inclusive, gap a GFF3 Gap attribute string
    (eg 'M185 I3 M250') or None for a perfect alignment. See docs/coordinates_and_exons.md """


def create_perfect_exons(raw_exon_stranded_order):
    """ Perfectly matched exons are basically a no-gap case of cDNA match """
    exons = []
    cdna_start = 1
    exon_id = 0
    for exon_start, exon_end in raw_exon_stranded_order:
        exon_length = exon_end - exon_start
        cdna_end = cdna_start + exon_length - 1
        exons.append((exon_start, exon_end, exon_id, cdna_start, cdna_end, None))
        cdna_start = cdna_end + 1
        exon_id += 1
    return exons


def create_cdna_exons(cdna_matches_stranded_order):
    """ Adds on exon_id """
    exons = []
    exon_id = 0
    for (exon_start, exon_end, cdna_start, cdna_end, gap) in cdna_matches_stranded_order:
        exons.append((exon_start, exon_end, exon_id, cdna_start, cdna_end, gap))
        exon_id += 1
    return exons


def get_cdna_match_offset(cdna_match_gap, position: int, validate=True, end=False):
    """ cdna_match GAP attribute looks like: 'M185 I3 M250' which is code/length
        @see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md#the-gap-attribute
        codes operation
        M 	match
        I 	insert a gap into the reference sequence
        D 	insert a gap into the target (delete from reference)

        position is a 0-based genomic offset into the exon, in transcript direction.
        With end=True, position is the offset after the last base of a range (eg the stop codon),
        so it may sit right before a D. If you want the whole exon, then pass the exon's genomic length
    """

    if not cdna_match_gap:
        return 0

    genomic_consumed = 0  # M and D both use up genomic bases
    offset = 0
    for gap_op in cdna_match_gap.split():
        code = gap_op[0]
        length = int(gap_op[1:])
        if code == "M":
            if position < genomic_consumed + length:
                break
            genomic_consumed += length
        elif code == "I":
            offset += length
        elif code == "D":
            if end and position == genomic_consumed:
                break  # the range ends right before the D
            if validate and position < genomic_consumed + length:
                raise ValueError(
                    "Coordinate (%d) inside deletion (%s) - no mapping possible!" % (position + 1, gap_op))
            genomic_consumed += length
            offset -= length
        else:
            raise ValueError("Unknown code in cDNA GAP: %s" % gap_op)

    return offset


def get_transcript_position(transcript_strand, ordered_cdna_matches, genomic_coordinate, label=None, end=False):
    """ Returns a 0-based position along the whole transcript (issue #123)

        With end=True, genomic_coordinate is the end of a range in transcript direction (eg the stop
        codon), see get_cdna_match_offset

        The exon's own cdna_start is used as the offset, rather than the running sum of the
        preceding exon lengths. These are the same thing for the vast majority of transcripts,
        but a handful of RefSeq alignments leave a hole in the transcript coordinates - a run
        of transcript bases that aligns nowhere on the genome, so exon N+1 starts later than
        exon N ended. Summing exon lengths silently collapses those holes out, putting the
        codon positions in a different coordinate system to the exon cds_start/cds_end. """
    for (exon_start, exon_end, _exon_id, cdna_start, _cdna_end, cdna_match_gap) in ordered_cdna_matches:
        if exon_start <= genomic_coordinate <= exon_end:
            # We're inside this match
            if transcript_strand:
                position = genomic_coordinate - exon_start
            else:
                position = exon_end - genomic_coordinate
            # cdna_start is 1-based, so cdna_start - 1 is the exon's 0-based transcript start
            return (cdna_start - 1) + position + get_cdna_match_offset(cdna_match_gap, position, end=end)
    if label is None:
        label = "Genomic coordinate: %d" % genomic_coordinate
    raise ValueError('%s is not in any of the exons' % label)


def cdna_offset_to_genomic_offset(gap, cdna_offset):
    """ Convert an offset within an exon from cDNA to genomic, both 0-based in transcript
        direction. gap is a GFF3-style gap string (eg 'M185 I3 M250'): M consumes both,
        I cDNA only, D genomic only """
    if not gap:
        return cdna_offset
    cdna_consumed = 0
    genomic_consumed = 0
    for op_str in gap.split():
        code = op_str[0]
        length = int(op_str[1:])
        if code == "M":
            if cdna_offset < cdna_consumed + length:
                return genomic_consumed + (cdna_offset - cdna_consumed)
            cdna_consumed += length
            genomic_consumed += length
        elif code == "I":
            if cdna_offset < cdna_consumed + length:
                raise ValueError(f"cDNA offset {cdna_offset} is in gap '{op_str}' (unaligned transcript bases)")
            cdna_consumed += length
        elif code == "D":
            genomic_consumed += length
        else:
            raise ValueError(f"Unknown gap operation '{op_str}'")
    return genomic_consumed + (cdna_offset - cdna_consumed)


def transcript_position_to_genomic(strand, exons, transcript_position):
    """ The inverse of get_transcript_position: transcript_position is 0-based along the whole
        transcript, returns the 0-based genomic coordinate of that base """
    cdna_position = transcript_position + 1  # exon cdna_start/cdna_end are 1-based inclusive
    for (alt_start, alt_end, _exon_id, cdna_start, cdna_end, gap) in exons:
        if cdna_start <= cdna_position <= cdna_end:
            genomic_offset = cdna_offset_to_genomic_offset(gap, cdna_position - cdna_start)
            if strand == '+':
                return alt_start + genomic_offset
            else:
                return alt_end - 1 - genomic_offset
    raise ValueError(f"Transcript position {transcript_position} is not in any of the exons")
