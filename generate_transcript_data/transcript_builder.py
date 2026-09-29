""" Collects the genes, transcripts and their exon/CDS rows from a GTF/GFF3 (whatever the format or
    consortium) and turns them into cdot transcript dicts: exon arrays in genomic order, start/stop
    codon transcript positions, translation exceptions, and the per genome build layout.

    Issue #101 """
import importlib
import logging
import operator
from collections import defaultdict

from generate_transcript_data.transcript_coordinates import create_perfect_exons, create_cdna_exons, \
    get_transcript_position

CONTIG = "contig"
STRAND = "strand"
CODING_FEATURES = {"CDS", "start_codon", "stop_codon"}  # Use these to work out cds_start/cds_end

# These are fields that differ per transcript/genome build, anything NOT in here should be a property of
# the transcript across all builds
GENOME_BUILD_FIELDS = ["cds_start", "cds_end", "strand", "contig", "exons", "other_chroms", "source",
                       "tag", "note", "ccds", "transcript_support_level", "genome_mismatch", "warnings"]


class TranscriptBuilder:
    def __init__(self, genome_build, url, name_ac_map):
        self.genome_build = genome_build
        self.url = url
        self.name_ac_map = name_ac_map
        self.genes = {}
        self.transcripts = {}
        self.transcript_proteins = {}
        # Store features in separate dict as we don't need to write all as JSON
        self.transcript_features_by_type = defaultdict(lambda: defaultdict(list))
        self._warned_about_htseq_tag_attributes = False

    # Genes / transcripts

    def get_or_create_gene(self, feature, gene_accession):
        gene_data = self.genes.get(gene_accession)
        if gene_data is None:
            gene_name = feature.attr.get("gene_name") or feature.attr.get("Name")
            description = feature.attr.get("description")

            biotype_set = set()
            biotype = feature.attr.get("gene_biotype") or feature.attr.get("biotype")
            if biotype:
                biotype_set.add(biotype)

            source_set = set()
            if feature.source:
                source_set.add(feature.source)

            gene_data = {
                "gene_symbol": gene_name,
                "biotype": biotype_set,
                "source": source_set,
                "id": gene_accession,
                "description": description
            }
            self.genes[gene_accession] = gene_data
        return gene_data

    def create_transcript(self, feature, transcript_accession, gene_data):
        transcript_data = {
            "id": transcript_accession,
            "gene_name": gene_data.get("gene_symbol"),
            "gene_version": gene_data.get("id"),
            "exons": [],
            "biotype": set(),
            "source": set(),
            CONTIG: feature.iv.chrom,
            STRAND: feature.iv.strand,
        }
        if hgnc := gene_data.get("hgnc"):
            transcript_data["hgnc"] = hgnc
        self.transcripts[transcript_accession] = transcript_data
        return transcript_data

    @staticmethod
    def store_other_chrom(data, feature):
        other_chroms = data.get("other_chroms", set())
        other_chroms.add(feature.iv.chrom)
        data["other_chroms"] = other_chroms

    def set_protein(self, transcript_accession, protein_accession):
        self.transcript_proteins[transcript_accession] = protein_accession

    def add_tags(self, transcript_data, feature):
        # Ideally we only want to get this once per transcript
        # So we want to pick something like mRNA or transcript not CDS or exon
        if feature.type in ('mRNA', 'transcript'):
            attr_tuples = getattr(feature, "attr_tuples", None)
            if attr_tuples is None:
                attr_tuples = feature.attr.items()
                if not self._warned_about_htseq_tag_attributes:
                    self._warned_about_htseq_tag_attributes = True
                    htseq_version = importlib.metadata.version('HTSeq')
                    logging.warning("Your version of HTSeq (%s) can not handle duplicated tags. Some will be lost "
                                    "See https://github.com/htseq/htseq/issues/83", htseq_version)

            attr_list_vals = defaultdict(list)
            for tag, value in attr_tuples:
                attr_list_vals[tag].append(value)

            if tag_list := attr_list_vals.get("tag"):
                transcript_data["tag"] = ",".join(tag_list)

    # Exon / CDS rows

    def add_feature(self, transcript_accession, transcript, feature) -> bool:
        """ An exon/CDS/codon/cDNA_match row. Returns False if it is on another contig to the
            transcript, in which case only the contig is recorded """
        if feature.iv.chrom != transcript[CONTIG]:
            self.store_other_chrom(transcript, feature)
            return False

        if feature.type == "cDNA_match":
            target = feature.attr.get("Target")  # Target=NM_001304717.2 1 1110 +
            t_cols = target.split()
            cdna_start = int(t_cols[1])  # These are 1-based (as per above)
            cdna_end = int(t_cols[2])
            gap = feature.attr.get("Gap")
            feature_tuple = (feature.iv.start, feature.iv.end, cdna_start, cdna_end, gap)
        else:
            feature_tuple = (feature.iv.start, feature.iv.end)

        features_by_type = self.transcript_features_by_type[transcript_accession]
        features_by_type[feature.type].append(feature_tuple)
        if feature.type in CODING_FEATURES:
            features_by_type["coding_starts"].append(feature.iv.start)
            features_by_type["coding_ends"].append(feature.iv.end)

        if note := feature.attr.get("Note"):
            transcript["note"] = note
        return True

    def add_exon_coordinates(self, transcript_accession, feature):
        """ Use a row's coordinates as an exon, without the row itself being one """
        self.transcript_features_by_type[transcript_accession]["exon"].append((feature.iv.start, feature.iv.end))

    # Translation annotations, converted to codon numbers in finish() once the exons and CDS are known

    def add_transl_except(self, transcript_accession, start, end, amino_acid):
        """ A codon that codes for amino_acid rather than what the translation table says """
        self.transcript_features_by_type[transcript_accession]["transl_except"].append((start, end, amino_acid))

    def set_transl_table(self, transcript_accession, transl_table: int):
        self.transcript_features_by_type[transcript_accession]["transl_table"] = [transl_table]

    def add_ribosomal_slippage_cds(self, transcript_accession, feature):
        """ A CDS row marked as having a ribosomal frameshift """
        feature_tuple = (feature.iv.start, feature.iv.end)
        self.transcript_features_by_type[transcript_accession]["ribosomal_slippage_cds"].append(feature_tuple)

    def add_translation_exception(self, transcript_accession, exception):
        self.transcript_features_by_type[transcript_accession]["translation_exceptions"].append(exception)

    def add_genome_mismatch_exception(self, transcript_accession, exception):
        self.transcript_features_by_type[transcript_accession]["genome_mismatch_exceptions"].append(exception)

    # Finishing

    def finish(self):
        """ Returns (genes, transcripts) in cdot layout """
        self._build_transcripts()

        for gene_data in self.genes.values():
            # TODO: Can turn this back on if we want - just removing for diff
            gene_data.pop("id", None)
            gene_data["url"] = self.url

        # At the moment the transcript dict is flat - need to move it into "genome_builds" dict
        for transcript_accession, transcript_data in self.transcripts.items():
            if protein := self.transcript_proteins.get(transcript_accession):
                transcript_data["protein"] = protein

            genome_build_coordinates = {
                "url": self.url,
            }
            for field in GENOME_BUILD_FIELDS:
                value = transcript_data.pop(field, None)
                if value is not None:
                    genome_build_coordinates[field] = value

            # Make sure contig uses accession not chromosome names
            contig = genome_build_coordinates["contig"]
            genome_build_coordinates["contig"] = self.name_ac_map.get(contig, contig)
            transcript_data["genome_builds"] = {self.genome_build: genome_build_coordinates}

        return self.genes, self.transcripts

    def _build_transcripts(self):
        for transcript_accession, transcript_data in self.transcripts.items():
            features_by_type = self.transcript_features_by_type.get(transcript_accession, {})

            # Store coding start/stop transcript positions
            # For RefSeq, we need to deal with alignment gaps, so easiest is to convert exons w/o gaps
            # into cDNA match objects, so the same objects/algorithm can be used
            forward_strand = transcript_data[STRAND] == '+'
            cdna_matches = features_by_type.get("cDNA_match")
            if cdna_matches:
                cdna_matches_stranded_order = cdna_matches
                cdna_matches_stranded_order.sort(key=operator.itemgetter(0))
                if not forward_strand:
                    cdna_matches_stranded_order.reverse()
                # Need to add exon ID
                exons_stranded_order = create_cdna_exons(cdna_matches_stranded_order)

            else:
                raw_exon_stranded_order = features_by_type.get("exon", [])
                raw_exon_stranded_order.sort(key=operator.itemgetter(0))
                if not forward_strand:
                    raw_exon_stranded_order.reverse()
                exons_stranded_order = create_perfect_exons(raw_exon_stranded_order)

            if "coding_starts" in features_by_type:
                cds_min = min(features_by_type["coding_starts"])
                cds_max = max(features_by_type["coding_ends"])

                transcript_data["cds_start"] = cds_min
                transcript_data["cds_end"] = cds_max

                (coding_left, coding_right) = ("start_codon", "stop_codon")
                if not forward_strand:  # Switch
                    (coding_left, coding_right) = (coding_right, coding_left)

                try:
                    transcript_data[coding_left] = get_transcript_position(forward_strand, exons_stranded_order,
                                                                           cds_min)
                except ValueError as e:
                    logging.warning("Couldn't set %s transcript position from %s: %s", coding_left, cds_min, e)
                    self._add_warning(transcript_data, "codons_unplaced", coding_left)

                try:
                    transcript_data[coding_right] = get_transcript_position(forward_strand, exons_stranded_order,
                                                                            cds_max)
                except ValueError as e:
                    logging.warning("Couldn't set %s transcript positions from %s: %s", coding_right, cds_max, e)
                    self._add_warning(transcript_data, "codons_unplaced", coding_right)

            self._add_translation(transcript_accession, transcript_data, features_by_type, forward_strand,
                                  exons_stranded_order)

            exons_genomic_order = exons_stranded_order
            if not forward_strand:
                exons_genomic_order.reverse()
            transcript_data["exons"] = exons_genomic_order

            if "cds_start" in transcript_data:
                biotype = "mRNA"
            else:
                biotype = "ncRNA"
            transcript_data["biotype"].add(biotype)

            if gene_accession := transcript_data.get("gene_version"):
                gene_data = self.genes[gene_accession]
                gene_data["biotype"].update(transcript_data["biotype"])

    @staticmethod
    def _add_warning(transcript_data, warning, value=None):
        """ Problems converting this transcript (stored per genome build), keyed by warning type.
            Warnings with a value are a sorted list of them, otherwise True """
        warnings = transcript_data.setdefault("warnings", {})
        if value is None:
            warnings[warning] = True
        else:
            warnings[warning] = sorted(set(warnings.get(warning, [])) | {value})

    @staticmethod
    def _is_transcript_transl_except(amino_acid, codon):
        """ Sec, stops (TERM, eg completed by the poly(A) tail), 'Other' (eg stop codon readthrough) and
            non-AUG starts are properties of the transcript. The rest are RefSeq correcting a codon where this
            genome differs from the transcript, which depends on the genome build """
        return amino_acid in {"Sec", "TERM", "Other"} or codon == 1

    def _add_translation(self, transcript_accession, transcript_data, features_by_type, forward_strand,
                         exons_stranded_order):
        translation = {}
        genome_mismatch = {}
        if transl_table := features_by_type.get("transl_table"):
            translation["transl_table"] = transl_table[0]

        has_codons = "start_codon" in transcript_data and "stop_codon" in transcript_data
        if transl_except := features_by_type.get("transl_except"):
            if has_codons:
                codons, unplaced = self._get_transl_except_codons(transcript_accession, transcript_data,
                                                                  forward_strand, exons_stranded_order,
                                                                  transl_except)
                for amino_acid, codon_numbers in codons.items():
                    for codon in codon_numbers:
                        if self._is_transcript_transl_except(amino_acid, codon):
                            d = translation
                        else:
                            d = genome_mismatch
                        d.setdefault("transl_except", {}).setdefault(amino_acid, []).append(codon)
            else:
                unplaced = {amino_acid for _, _, amino_acid in transl_except}
            # Consumers can't tell a missing transl_except from none, so record it
            for amino_acid in unplaced:
                self._add_warning(transcript_data, "transl_except_unplaced", amino_acid)

        if slippage_cds := features_by_type.get("ribosomal_slippage_cds"):
            unplaced = True
            if has_codons:
                slippage, unplaced = self._get_ribosomal_slippage(transcript_accession, transcript_data,
                                                                  forward_strand, exons_stranded_order,
                                                                  slippage_cds)
                if slippage:
                    translation["ribosomal_slippage"] = slippage
            if unplaced:
                self._add_warning(transcript_data, "ribosomal_slippage_unplaced")

        if exceptions := features_by_type.get("translation_exceptions"):
            translation["exceptions"] = sorted(set(exceptions))
        if exceptions := features_by_type.get("genome_mismatch_exceptions"):
            genome_mismatch["exceptions"] = sorted(set(exceptions))

        if translation:
            transcript_data["translation"] = translation
        if genome_mismatch:
            transcript_data["genome_mismatch"] = genome_mismatch

    @staticmethod
    def _get_transl_except_codons(transcript_accession, transcript_data, forward_strand, exons_stranded_order,
                                  transl_except):
        """ Returns ({amino_acid: [codon numbers]}, {unplaced amino acids}), the codon numbers 1-based
            within the CDS (ie the amino acid positions in the protein) """
        codons_by_amino_acid = defaultdict(set)
        cds_length = transcript_data["stop_codon"] - transcript_data["start_codon"]
        not_codon_starts = []  # (amino_acid, cds_position or None, start, end)
        for start, end, amino_acid in transl_except:
            # First base of the codon, in transcript direction
            genomic_coordinate = start if forward_strand else end
            try:
                transcript_position = get_transcript_position(forward_strand, exons_stranded_order,
                                                              genomic_coordinate)
            except ValueError as e:
                logging.warning("%s: couldn't place transl_except %s at %d-%d: %s",
                                transcript_accession, amino_acid, start + 1, end, e)
                not_codon_starts.append((amino_acid, None, start, end))
                continue
            cds_position = transcript_position - transcript_data["start_codon"]
            if 0 <= cds_position < cds_length and cds_position % 3 == 0:
                codons_by_amino_acid[amino_acid].add(cds_position // 3 + 1)
            else:
                not_codon_starts.append((amino_acid, cds_position, start, end))

        unplaced = set()
        for amino_acid, cds_position, start, end in not_codon_starts:
            if cds_position is not None:
                if 0 <= cds_position < cds_length and cds_position // 3 + 1 in codons_by_amino_acid[amino_acid]:
                    continue  # Rest of a codon split across exons (Ensembl writes a row per exon)
                logging.warning("%s: transl_except %s at %d-%d is not a codon of the CDS",
                                transcript_accession, amino_acid, start + 1, end)
            unplaced.add(amino_acid)
        codons = {amino_acid: sorted(codons) for amino_acid, codons in codons_by_amino_acid.items() if codons}
        return codons, unplaced

    @staticmethod
    def _get_ribosomal_slippage(transcript_accession, transcript_data, forward_strand, exons_stranded_order,
                                slippage_cds):
        """ RefSeq marks a ribosomal frameshift with 'exception=ribosomal slippage' on the CDS rows, and
            splits the CDS there: rows overlapping by 1 base is a -1 frameshift (that base is read twice),
            a 1 base gap between rows is a +1 frameshift (that base is skipped).
            Returns ([{"cds_position": 1-based position in the CDS (as c. numbering), "shift": -1 or 1}],
                     any_unplaced) """
        slippage = []
        unplaced = False
        slippage_cds = sorted(set(slippage_cds))
        for (_, prev_end), (next_start, _) in zip(slippage_cds, slippage_cds[1:]):
            if next_start == prev_end - 1:
                shift, base = -1, prev_end - 1
            elif next_start == prev_end + 1:
                shift, base = 1, prev_end
            else:
                continue  # an intron
            try:
                genomic_coordinate = base if forward_strand else base + 1
                transcript_position = get_transcript_position(forward_strand, exons_stranded_order,
                                                              genomic_coordinate)
            except ValueError as e:
                logging.warning("%s: couldn't place ribosomal slippage at %d: %s", transcript_accession, base + 1, e)
                unplaced = True
                continue
            cds_position = transcript_position - transcript_data["start_codon"] + 1
            slippage.append({"cds_position": cds_position, "shift": shift})
        if not slippage:
            unplaced = True
            logging.warning("%s: CDS has a ribosomal slippage exception, but couldn't find where",
                            transcript_accession)
        return sorted(slippage, key=lambda s: s["cds_position"]), unplaced
