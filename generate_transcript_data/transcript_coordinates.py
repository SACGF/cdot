""" Exon / transcript coordinate maths shared by the GTF/GFF3 parsers and the UTA converter

    Exons are tuples of (alt_start, alt_end, exon_id, cdna_start, cdna_end, gap): genomic coordinates
    0-based half-open, transcript (cdna) coordinates 1-based inclusive, gap a GFF3 Gap attribute string
    (eg 'M185 I3 M250') or None for a perfect alignment. See docs/coordinates_and_exons.md

    Coordinate conventions in this module:

    * Genomic coordinates are 0-based half-open (interbase) on both strands, as HTSeq reads them from
      the GTF/GFF3 and as cdot writes them to the exon array. A range is (start, end) with start < end
      in genomic order whatever the strand, so the CDS, a codon or a single base all arrive here as a
      genomic range, never as a bare coordinate.
    * Transcript positions are 0-based half-open along the whole transcript, the start_codon /
      stop_codon fields (biocommons cds_start_i / cds_end_i).
    * The only place a strand is used to turn a genomic range into a transcript position is
      _range_boundary_transcript_position. Callers never add or subtract 1 for the strand.
"""


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


def get_cdna_match_offset(cdna_match_gap, base_offset: int):
    """ cdna_match GAP attribute looks like: 'M185 I3 M250' which is code/length
        @see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md#the-gap-attribute
        codes operation
        M 	match
        I 	insert a gap into the reference sequence
        D 	insert a gap into the target (delete from reference)

        base_offset is the 0-based genomic offset of a base into the exon, in transcript direction.
        Returns what to add to it to get the base's cDNA offset into the exon. Raises ValueError if
        the base is in a D (genome only, so it has no cDNA position)
    """

    if not cdna_match_gap:
        return 0

    genomic_consumed = 0  # M and D both use up genomic bases
    offset = 0
    for gap_op in cdna_match_gap.split():
        code = gap_op[0]
        length = int(gap_op[1:])
        if code == "M":
            if base_offset < genomic_consumed + length:
                break
            genomic_consumed += length
        elif code == "I":
            offset += length
        elif code == "D":
            if base_offset < genomic_consumed + length:
                raise ValueError(
                    "Coordinate (%d) inside deletion (%s) - no mapping possible!" % (base_offset + 1, gap_op))
            genomic_consumed += length
            offset -= length
        else:
            raise ValueError("Unknown code in cDNA GAP: %s" % gap_op)

    return offset


def get_transcript_range_start(transcript_strand, ordered_cdna_matches, genomic_start, genomic_end):
    """ The transcript position (0-based) of the first base of a genomic range (0-based half-open),
        eg the CDS gives the start_codon. On the minus strand the first base in transcript direction
        is the one before genomic_end. Raises ValueError if that base is not in an exon (issue #123) """
    return _range_boundary_transcript_position(transcript_strand, ordered_cdna_matches,
                                               genomic_start, genomic_end, range_start=True)


def get_transcript_range_end(transcript_strand, ordered_cdna_matches, genomic_start, genomic_end):
    """ The transcript position (0-based, half-open) after the last base of a genomic range (0-based
        half-open), eg the CDS gives the stop_codon. On the minus strand the last base in transcript
        direction is genomic_start. Raises ValueError if that base is not in an exon """
    return _range_boundary_transcript_position(transcript_strand, ordered_cdna_matches,
                                               genomic_start, genomic_end, range_start=False)


def _range_boundary_transcript_position(transcript_strand, ordered_cdna_matches, genomic_start, genomic_end,
                                        range_start):
    """ A genomic boundary (an interbase coordinate, between two bases) can have two transcript
        positions: at a junction of abutting exons, and where transcript-only bases (Gap I) sit against
        it. The base of the range next to the boundary is unambiguous, so the boundary is placed by that
        base: the range's first base (in transcript direction) for its start, its last for its end.

        The exon's own cdna_start is used as the offset, rather than the running sum of the
        preceding exon lengths. These are the same thing for the vast majority of transcripts,
        but a handful of RefSeq alignments leave a hole in the transcript coordinates - a run
        of transcript bases that aligns nowhere on the genome, so exon N+1 starts later than
        exon N ended. Summing exon lengths silently collapses those holes out, putting the
        codon positions in a different coordinate system to the exon cds_start/cds_end. """
    if transcript_strand:
        base = genomic_start if range_start else genomic_end - 1
    else:
        base = genomic_end - 1 if range_start else genomic_start

    for (exon_start, exon_end, _exon_id, cdna_start, _cdna_end, cdna_match_gap) in ordered_cdna_matches:
        if exon_start <= base < exon_end:
            # The base's offset into the exon, in transcript direction
            if transcript_strand:
                base_offset = base - exon_start
            else:
                base_offset = exon_end - 1 - base
            # cdna_start is 1-based, so cdna_start - 1 is the exon's 0-based transcript start
            transcript_position = (cdna_start - 1) + base_offset + get_cdna_match_offset(cdna_match_gap, base_offset)
            if not range_start:
                transcript_position += 1  # Half-open end: after the last base
            return transcript_position
    which = "first" if range_start else "last"
    raise ValueError(f"Genomic range {genomic_start}-{genomic_end}: {which} base ({base}) is not in any of the exons")


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
    """ transcript_position is 0-based along the whole transcript, returns the 0-based genomic
        coordinate of that base on either strand. The inverse is get_transcript_range_start on the
        one base range (coordinate, coordinate + 1) """
    cdna_position = transcript_position + 1  # exon cdna_start/cdna_end are 1-based inclusive
    for (alt_start, alt_end, _exon_id, cdna_start, cdna_end, gap) in exons:
        if cdna_start <= cdna_position <= cdna_end:
            genomic_offset = cdna_offset_to_genomic_offset(gap, cdna_position - cdna_start)
            if strand == '+':
                return alt_start + genomic_offset
            else:
                return alt_end - 1 - genomic_offset
    raise ValueError(f"Transcript position {transcript_position} is not in any of the exons")
