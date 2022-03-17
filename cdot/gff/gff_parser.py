import abc
import logging
import operator
import os
import re
from collections import Counter, defaultdict
from hashlib import md5

import HTSeq

# Keys used in dictionary (serialized to JSON)
CONTIG = "contig"
START = "start"
END = "stop"
STRAND = "strand"


def file_md5sum(filename):
    m = md5()
    with open(filename, "rb") as f:
        m.update(f.read())
    return m.hexdigest()


class GFFParser(abc.ABC):
    CODING_FEATURES = {"CDS", "start_codon", "stop_codon"}  # Use these to work out cds_start/cds_end
    FEATURE_ALLOW_LIST = {}
    FEATURE_IGNORE_LIST = {"biological_region", "chromosome", "region", "scaffold", "supercontig"}

    def __init__(self, filename, discard_contigs_with_underscores=True):
        self.filename = filename
        self.discard_contigs_with_underscores = discard_contigs_with_underscores

        self.discarded_contigs = Counter()
        self.genes_by_id = {}
        self.transcripts_by_id = {}
        self.gene_id_by_name = {}
        # Store features in separate dict as we don't need to write all as JSON
        self.transcript_features_by_type = defaultdict(lambda: defaultdict(list))

    @abc.abstractmethod
    def handle_feature(self, feature):
        pass

    def parse(self):
        for feature in HTSeq.GFF_Reader(self.filename):
            if self.FEATURE_ALLOW_LIST and feature.type not in self.FEATURE_ALLOW_LIST:
                continue
            if feature.type in self.FEATURE_IGNORE_LIST:
                continue

            try:
                contig = feature.iv.chrom
                if self.discard_contigs_with_underscores and not contig.startswith("NC_") and "_" in contig:
                    self.discarded_contigs[contig] += 1
                    continue
                self.handle_feature(feature)
            except Exception as e:
                print("Could not parse '%s': %s" % (feature.get_gff_line(), e))
                raise e

    def finish(self):
        self._process_coding_features()

        if self.discarded_contigs:
            print("Discarded contigs: %s" % self.discarded_contigs)

    @staticmethod
    def _create_gene(gene_name, feature):
        biotypes = set()

        gene = {
            "name": gene_name,
            "transcripts": set(),
            "biotype": biotypes,
            CONTIG: feature.iv.chrom,
            START: feature.iv.start,
            END: feature.iv.end,
            STRAND: feature.iv.strand
        }

        # Attempt to get some biotypes in there if available
        if feature.type == "gene":
            gene_version = feature.attr.get("version")
            biotype = feature.attr.get("biotype")
            description = feature.attr.get("description")
            if description:
                gene["description"] = description
        else:
            gene_version = feature.attr.get("gene_version")
            biotype = feature.attr.get("gene_biotype")

        if biotype:
            biotypes.add(biotype)

        if gene_version:
            gene["version"] = int(gene_version)
        return gene

    @staticmethod
    def _create_transcript(feature):
        return {
            "exons": [],
            "biotype": set(),
            CONTIG: feature.iv.chrom,
            START: feature.iv.start,
            END: feature.iv.end,
            STRAND: feature.iv.strand,
        }

    @staticmethod
    def _store_other_chrom(data, feature):
        other_chroms = data.get("other_chroms", set())
        other_chroms.add(feature.iv.chrom)
        data["other_chroms"] = other_chroms

    @staticmethod
    def _get_biotype_from_transcript_id(transcript_id):
        biotypes_by_transcript_id_start = {"NM_": "protein_coding", "NR_": "non_coding"}
        for (start, biotype) in biotypes_by_transcript_id_start.items():
            if transcript_id.startswith(start):
                return biotype

        if "tRNA" in transcript_id:
            return "tRNA"
        return None

    def _add_transcript_data(self, transcript_id, transcript, feature):
        if feature.iv.chrom != transcript[CONTIG]:
            self._store_other_chrom(transcript, feature)
            return

        if feature.type == "cDNA_match":
            target = feature.attr.get("Target")
            t_cols = target.split()
            cdna_start = int(t_cols[1])
            cdna_end = int(t_cols[2])
            gap = feature.attr.get("Gap")
            feature_tuple = (feature.iv.start, feature.iv.end, cdna_start, cdna_end, gap)
        else:
            feature_tuple = (feature.iv.start, feature.iv.end)

        features_by_type = self.transcript_features_by_type[transcript_id]
        features_by_type[feature.type].append(feature_tuple)
        if feature.type in self.CODING_FEATURES:
            features_by_type["coding_starts"].append(feature.iv.start)
            features_by_type["coding_ends"].append(feature.iv.end)

    def _process_coding_features(self):
        for transcript_id, transcript in self.transcripts_by_id.items():
            features_by_type = self.transcript_features_by_type.get(transcript_id)

            # Store coding start/stop transcript positions
            # For RefSeq, we need to deal with alignment gaps, so easiest is to convert exons w/o gaps
            # into cDNA match objects, so the same objects/algorithm can be used
            forward_strand = transcript[STRAND] == '+'
            cdna_matches = features_by_type.get("cDNA_match")
            if cdna_matches:
                cdna_matches_stranded_order = cdna_matches
                cdna_matches_stranded_order.sort(key=operator.itemgetter(0))
                if not forward_strand:
                    cdna_matches_stranded_order.reverse()
                # Need to add exon ID
                exons_stranded_order = self._create_cdna_exons(cdna_matches_stranded_order)

            else:
                raw_exon_stranded_order = features_by_type["exon"]
                raw_exon_stranded_order.sort(key=operator.itemgetter(0))
                if not forward_strand:
                    raw_exon_stranded_order.reverse()
                exons_stranded_order = self._create_perfect_exons(raw_exon_stranded_order)

            if "coding_starts" in features_by_type:
                cds_min = min(features_by_type["coding_starts"])
                cds_max = max(features_by_type["coding_ends"])

                transcript["cds_start"] = cds_min
                transcript["cds_end"] = cds_max

                try:
                    (coding_left, coding_right) = ("start_codon", "stop_codon")
                    if not forward_strand:  # Switch
                        (coding_left, coding_right) = (coding_right, coding_left)
                    transcript[coding_left] = GFFParser._get_transcript_position(forward_strand, exons_stranded_order,
                                                                                 cds_min)
                    transcript[coding_right] = GFFParser._get_transcript_position(forward_strand, exons_stranded_order,
                                                                                  cds_max)
                except Exception as e:
                    logging.warning("Couldn't set coding start/end transcript positions: %s", e)

            exons_genomic_order = exons_stranded_order
            if not forward_strand:
                exons_genomic_order.reverse()
            transcript["exons"] = exons_genomic_order

    @staticmethod
    def _create_perfect_exons(raw_exon_stranded_order):
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

    @staticmethod
    def _create_cdna_exons(cdna_matches_stranded_order):
        """ Adds on exon_id """
        exons = []
        exon_id = 0
        for (exon_start, exon_end, cdna_start, cdna_end, gap) in cdna_matches_stranded_order:
            exons.append((exon_start, exon_end, exon_id, cdna_start, cdna_end, gap))
            exon_id += 1
        return exons

    @staticmethod
    def get_cdna_match_offset(cdna_match_gap, position: int, validate=True):
        """ cdna_match GAP attribute looks like: 'M185 I3 M250' which is code/length
            @see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md#the-gap-attribute
            codes operation
            M 	match
            I 	insert a gap into the reference sequence
            D 	insert a gap into the target (delete from reference)

            If you want the whole exon, then pass the end
        """

        if not cdna_match_gap:
            return 0

        position_1_based = position + 1
        cdna_match_index = 1
        offset = 0
        for gap_op in cdna_match_gap.split():
            code = gap_op[0]
            length = int(gap_op[1:])
            if code == "M":
                cdna_match_index += length
            elif code == "I":
                if validate and position_1_based < cdna_match_index + length:
                    raise ValueError(
                        "Coordinate (%d) inside insertion (%s) - no mapping possible!" % (position_1_based, gap_op))
                offset += length
            elif code == "D":
                if validate and position < cdna_match_index + length:
                    raise ValueError(
                        "Coordinate (%d) inside deletion (%s) - no mapping possible!" % (position_1_based, gap_op))
                offset -= length
            else:
                raise ValueError("Unknown code in cDNA GAP: %s" % gap_op)

            if cdna_match_index > position_1_based:
                break

        return offset

    @staticmethod
    def _get_transcript_position(transcript_strand, ordered_cdna_matches, genomic_coordinate, label=None):
        cdna_offset = 0
        for (exon_start, exon_end, _exon_id, cdna_start, cdna_end, cdna_match_gap) in ordered_cdna_matches:
            if exon_start <= genomic_coordinate <= exon_end:
                # We're inside this match
                if transcript_strand:
                    position = genomic_coordinate - exon_start
                else:
                    position = exon_end - genomic_coordinate
                return cdna_offset + position + GFFParser.get_cdna_match_offset(cdna_match_gap, position)
            else:
                length = cdna_end - cdna_start + 1
                cdna_offset += length
        if label is None:
            label = "Genomic coordinate: %d" % genomic_coordinate
        raise ValueError('%s is not in any of the exons' % label)

    def get_data(self):
        self.parse()
        self.finish()

        gene_ids_by_biotype = defaultdict(set)
        for gene_id, gene in self.genes_by_id.items():
            for biotype in gene["biotype"]:
                gene_ids_by_biotype[biotype].add(gene_id)

        return {
            "reference_gtf": {"path": os.path.abspath(self.filename),
                              "md5sum": file_md5sum(self.filename)},
            "genes_by_id": self.genes_by_id,
            "transcripts_by_id": self.transcripts_by_id,
            "gene_id_by_name": self.gene_id_by_name,
            "gene_ids_by_biotype": gene_ids_by_biotype,
        }


class GTFParser(GFFParser):
    """ GTF (GFF2) - used by Ensembl, @see http://gmod.org/wiki/GFF2

        GFF2 only has 2 levels of feature hierarchy, so we have to build or 3 levels of gene/transcript/exons ourselves
    """
    GTF_TRANSCRIPTS_DATA = GFFParser.CODING_FEATURES | {"exon"}
    FEATURE_ALLOW_LIST = GTF_TRANSCRIPTS_DATA | {"gene"}

    def __init__(self, *args, **kwargs):
        super(GTFParser, self).__init__(*args, **kwargs)

    def handle_feature(self, feature):
        gene_id = feature.attr["gene_id"]
        # Non mandatory - Ensembl doesn't have on some RNAs
        gene_name = None
        if feature.type == "gene":
            gene_name = feature.attr.get("Name")
        else:
            gene_name = feature.attr.get("gene_name")
        if gene_name:
            self.gene_id_by_name[gene_name] = gene_id  # Shouldn't be dupes per file

        gene = self.genes_by_id.get(gene_id)
        if gene is None:
            gene = self._create_gene(gene_name, feature)
            self.genes_by_id[gene_id] = gene
        else:
            self._update_extents(gene, feature)

        transcript_id = feature.attr.get("transcript_id")
        transcript_version = feature.attr.get("transcript_version")
        if transcript_version:
            transcript_id += "." + transcript_version

        if transcript_id:
            gene["transcripts"].add(transcript_id)
            transcript = self.transcripts_by_id.get(transcript_id)
            if transcript is None:
                transcript = self._create_transcript(feature)
                self.transcripts_by_id[transcript_id] = transcript
            else:
                self._update_extents(transcript, feature)

            # No need to store chrom/strand for each feature, will use transcript
            if feature.type in self.GTF_TRANSCRIPTS_DATA:
                self._add_transcript_data(transcript_id, transcript, feature)

            biotype = feature.attr.get("gene_biotype")
            if biotype is None:
                # Ensembl GTFs store biotype info under gene_type or transcript_type
                biotype = feature.attr.get("gene_type")
                if biotype is None:
                    biotype = self._get_biotype_from_transcript_id(transcript_id)

            if biotype:
                gene["biotype"].add(biotype)
                transcript["biotype"].add(biotype)

    @staticmethod
    def _update_extents(genomic_region_dict, feature):
        if feature.iv.chrom == genomic_region_dict[CONTIG]:
            start = genomic_region_dict[START]
            if feature.iv.start < start:
                genomic_region_dict[START] = feature.iv.start

            end = genomic_region_dict[END]
            if feature.iv.end > end:
                genomic_region_dict[END] = feature.iv.end
        else:
            GFFParser._store_other_chrom(genomic_region_dict, feature)


class GFF3Parser(GFFParser):
    """ GFF3 - @see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        Used by RefSeq and later Ensembl (82 onwards)

        GFF3 support arbitrary hierarchy

    """

    GFF3_GENES = {"gene", "pseudogene"}
    GFF3_TRANSCRIPTS_DATA = {"exon", "CDS", "cDNA_match", "five_prime_UTR", "three_prime_UTR"}

    def __init__(self, *args, **kwargs):
        super(GFF3Parser, self).__init__(*args, **kwargs)
        self.gene_id_by_feature_id = defaultdict()
        self.transcript_id_by_feature_id = defaultdict()
        self.hgnc_pattern = re.compile(r".*\[Source:HGNC.*Acc:(HGNC:)?(\d+)\]")

    def handle_feature(self, feature):
        parent_id = feature.attr.get("Parent")
        # Genes never have parents
        # RefSeq genes are always one of GFF3_GENES, Ensembl has lots of different types (lincRNA_gene etc)
        # Ensembl treats pseudogene as a transcript (has parent)
        if parent_id is None and (feature.type in self.GFF3_GENES or "gene_id" in feature.attr):
            gene_id = feature.attr.get("gene_id")
            dbxref = self._get_dbxref(feature)
            if not gene_id:
                gene_id = dbxref.get("GeneID")
            if not gene_id:
                raise ValueError("Could not obtain 'gene_id', tried 'gene_id' and 'Dbxref[GeneID]'")

            gene_name = feature.attr.get("Name")
            # Gene can have multiple loci, thus entries in GFF, keep original so all transcripts are added
            gene = self.genes_by_id.get(gene_id)
            if gene is None:
                gene = self._create_gene(gene_name, feature)
                # If a gene already exists - then need to merge it...
                self.genes_by_id[gene_id] = gene

            # RefSeq stores HGNC in dbxref
            if hgnc := dbxref.get("HGNC"):
                # Might have HGNC: (5 characters) at start of it
                if hgnc.startswith("HGNC"):
                    hgnc = hgnc[5:]
                gene["hgnc"] = hgnc
            elif description := gene.get("description"):
                # Ensembl stores HGNC in description comment. Sometimes has "HGNC:" prefix eg:
                # Can be either "[Source:HGNC Symbol%3BAcc:HGNC:8907]" or "[Source:HGNC Symbol%3BAcc:37102]"
                if m := self.hgnc_pattern.match(description):
                    gene["hgnc"] = m.group(2)

            if gene_name:
                self.gene_id_by_name[gene_name] = gene_id
            self.gene_id_by_feature_id[feature.attr["ID"]] = gene_id
        else:
            if feature.type in self.GFF3_TRANSCRIPTS_DATA:
                if feature.type == 'cDNA_match':
                    target = feature.attr["Target"]
                    transcript_id = target.split()[0]
                else:
                    # Some exons etc may be for miRNAs that have no transcript ID, so skip those (won't have parent)
                    if parent_id:
                        transcript_id = self.transcript_id_by_feature_id.get(parent_id)
                    else:
                        logging.warning("Transcript data has no parent: %s" % feature.get_gff_line())
                        transcript_id = None

                if transcript_id:
                    transcript = self.transcripts_by_id[transcript_id]
                    self._handle_transcript_data(transcript_id, transcript, feature)
            else:
                # There are so many different transcript ontology terms just taking everything that
                # has a transcript_id and is child of gene (ie skip miRNA etc that is child of primary_transcript)
                transcript_id = feature.attr.get("transcript_id")
                if transcript_id:
                    transcript_version = feature.attr.get("version")
                    if transcript_version:
                        transcript_id += "." + transcript_version
                    assert parent_id is not None
                    gene_id = self.gene_id_by_feature_id.get(parent_id)
                    if not gene_id:
                        raise ValueError("Don't know how to handle feature type %s (not child of gene)" % feature.type)
                    gene = self.genes_by_id[gene_id]
                    self._handle_transcript(gene, transcript_id, feature)

    @staticmethod
    def _get_dbxref(feature):
        """ RefSeq stores attribute with more keys, eg: 'Dbxref=GeneID:7840,HGNC:HGNC:428,MIM:606844' """
        dbxref = {}
        dbxref_str = feature.attr.get("Dbxref")
        if dbxref_str:
            dbxref = dict(d.split(":", 1) for d in dbxref_str.split(","))
        return dbxref

    def _handle_transcript(self, gene, transcript_id, feature):
        """ Sometimes we can get multiple transcripts in the same file - just taking 1st """
        if transcript_id not in self.transcripts_by_id:
            # print("_handle_transcript(%s, %s)" % (gene, feature))
            gene["transcripts"].add(transcript_id)
            transcript = self._create_transcript(feature)
            biotype = self._get_biotype_from_transcript_id(transcript_id)
            if biotype:
                gene["biotype"].add(biotype)
                transcript["biotype"].add(biotype)
            partial = feature.attr.get("partial")
            if partial:
                transcript["partial"] = 1
            self.transcripts_by_id[transcript_id] = transcript
        self.transcript_id_by_feature_id[feature.attr["ID"]] = transcript_id

    def _handle_transcript_data(self, transcript_id, transcript, feature):
        self._add_transcript_data(transcript_id, transcript, feature)
