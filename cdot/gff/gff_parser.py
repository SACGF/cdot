import abc
import logging
import operator
import re
from collections import Counter, defaultdict
from typing import Optional

import HTSeq

# Keys used in dictionary (serialized to JSON)
from bioutils.assemblies import make_name_ac_map

CONTIG = "contig"
STRAND = "strand"


class GFFParser(abc.ABC):
    CODING_FEATURES = {"CDS", "start_codon", "stop_codon"}  # Use these to work out cds_start/cds_end
    FEATURE_ALLOW_LIST = {}
    FEATURE_IGNORE_LIST = {"biological_region", "chromosome", "region", "scaffold", "supercontig"}

    def __init__(self, filename, genome_build, url, discard_contigs_with_underscores=True):
        self.filename = filename
        self.genome_build = genome_build
        self.url = url
        self.discard_contigs_with_underscores = discard_contigs_with_underscores

        self.discarded_contigs = Counter()
        self.gene_data_by_accession = {}
        self.transcript_data_by_accession = {}
        # Store features in separate dict as we don't need to write all as JSON
        self.transcript_features_by_type = defaultdict(lambda: defaultdict(list))
        self.name_ac_map = make_name_ac_map(genome_build)

    @abc.abstractmethod
    def handle_feature(self, feature):
        pass

    def _parse(self):
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

    def _finish(self):
        self._process_coding_features()

        for gene_data in self.gene_data_by_accession.values():
            # TODO: Can turn this back on if we want - just removing for diff
            gene_data.pop("id", None)
            #############
            # Turn set into comma sep string
            gene_data["biotype"] = ",".join(sorted(gene_data["biotype"]))
            gene_data["url"] = self.url

        # At the moment the transcript dict is flat - need to move it into "genome_builds" dict
        GENOME_BUILD_FIELDS = ["cds_start", "cds_end", "strand", "contig", "exons", "other_chroms"]
        for transcript_data in self.transcript_data_by_accession.values():
            genome_build_coordinates = {
                "url": self.url,
            }
            for field in GENOME_BUILD_FIELDS:
                if value := transcript_data.pop(field, None):
                    genome_build_coordinates[field] = value

            # Make sure contig uses accession not chromosome names
            contig = genome_build_coordinates["contig"]
            genome_build_coordinates["contig"] = self.name_ac_map.get(contig, contig)
            transcript_data["genome_builds"] = {self.genome_build: genome_build_coordinates}

        if self.discarded_contigs:
            print("Discarded contigs: %s" % self.discarded_contigs)

    @staticmethod
    def _create_gene(feature, gene_accession):
        biotype_set = set()
        description = None

        # Non mandatory - Ensembl doesn't have some stuff on some RNAs
        if feature.type in {"gene", "pseudogene"}:
            gene_name = feature.attr.get("Name")
            description = feature.attr.get("description")
            biotype = feature.attr.get("biotype")
        else:
            gene_name = feature.attr.get("gene_name")
            biotype = feature.attr.get("gene_biotype")

        if biotype:
            biotype_set.add(biotype)

        return {
            "gene_symbol": gene_name,
            "biotype": biotype_set,
            "id": gene_accession,
            "description": description
        }

    @staticmethod
    def _create_transcript(feature, transcript_accession, gene_data):
        transcript_data = {
            "id": transcript_accession,
            "gene_name": gene_data["gene_symbol"],
            "gene_version": gene_data["id"],
            "exons": [],
            "biotype": set(),
            CONTIG: feature.iv.chrom,
            STRAND: feature.iv.strand,
        }
        if hgnc := gene_data.get("hgnc"):
            transcript_data["hgnc"] = hgnc
        return transcript_data

    @staticmethod
    def _store_other_chrom(data, feature):
        other_chroms = data.get("other_chroms", set())
        other_chroms.add(feature.iv.chrom)
        data["other_chroms"] = other_chroms

    @staticmethod
    def _get_biotype_from_transcript_accession(transcript_accession):
        biotypes_by_transcript_id_start = {"NM_": "protein_coding", "NR_": "non_coding"}
        for (start, biotype) in biotypes_by_transcript_id_start.items():
            if transcript_accession.startswith(start):
                return biotype

        if "tRNA" in transcript_accession:
            return "tRNA"
        return None

    def _add_transcript_data(self, transcript_accession, transcript, feature):
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

        features_by_type = self.transcript_features_by_type[transcript_accession]
        features_by_type[feature.type].append(feature_tuple)
        if feature.type in self.CODING_FEATURES:
            features_by_type["coding_starts"].append(feature.iv.start)
            features_by_type["coding_ends"].append(feature.iv.end)

    def _process_coding_features(self):
        for transcript_accession, transcript_data in self.transcript_data_by_accession.items():
            features_by_type = self.transcript_features_by_type.get(transcript_accession)

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

                transcript_data["cds_start"] = cds_min
                transcript_data["cds_end"] = cds_max

                try:
                    (coding_left, coding_right) = ("start_codon", "stop_codon")
                    if not forward_strand:  # Switch
                        (coding_left, coding_right) = (coding_right, coding_left)
                    transcript_data[coding_left] = GFFParser._get_transcript_position(forward_strand,
                                                                                      exons_stranded_order,
                                                                                      cds_min)
                    transcript_data[coding_right] = GFFParser._get_transcript_position(forward_strand,
                                                                                       exons_stranded_order,
                                                                                       cds_max)
                except Exception as e:
                    logging.warning("Couldn't set coding start/end transcript positions: %s", e)

            exons_genomic_order = exons_stranded_order
            if not forward_strand:
                exons_genomic_order.reverse()
            transcript_data["exons"] = exons_genomic_order

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

    @staticmethod
    def _get_transcript_accession(feature, version_key) -> Optional[str]:
        transcript_accession = None
        if transcript_id := feature.attr.get("transcript_id"):
            if transcript_version := feature.attr.get(version_key):
                transcript_accession = f"{transcript_id}.{transcript_version}"
            else:
                # print(f"warning: Couldn't get out {version_key} from {feature.type=} {feature.attr=}")
                transcript_accession = transcript_id
        return transcript_accession

    @staticmethod
    def _get_gene_accession(feature) -> Optional[str]:
        """ This can sometimes fail, in which case RefSeq will use dbxRef """
        gene_accession = None
        if gene_id := feature.attr.get("gene_id"):
            if feature.type == "gene":
                gene_version = feature.attr.get("version")
            else:
                gene_version = feature.attr.get("gene_version")

            if gene_version:
                gene_accession = f"{gene_id}.{gene_version}"
            else:
                gene_accession = gene_id
        return gene_accession


    def get_genes_and_transcripts(self):
        self._parse()
        self._finish()

        return self.gene_data_by_accession, self.transcript_data_by_accession


class GTFParser(GFFParser):
    """ GTF (GFF2) - used by Ensembl, @see http://gmod.org/wiki/GFF2

        GFF2 only has 2 levels of feature hierarchy, so we have to build or 3 levels of gene/transcript/exons ourselves
    """
    GTF_TRANSCRIPTS_DATA = GFFParser.CODING_FEATURES | {"exon"}
    FEATURE_ALLOW_LIST = GTF_TRANSCRIPTS_DATA | {"gene"}

    def __init__(self, *args, **kwargs):
        super(GTFParser, self).__init__(*args, **kwargs)

    def handle_feature(self, feature):
        gene_accession = self._get_gene_accession(feature)
        gene_data = self.gene_data_by_accession.get(gene_accession)
        if gene_data is None:
            gene_data = self._create_gene(feature, gene_accession)
            self.gene_data_by_accession[gene_accession] = gene_data

        if transcript_accession := self._get_transcript_accession(feature, version_key="transcript_version"):
            transcript = self.transcript_data_by_accession.get(transcript_accession)
            if transcript is None:
                transcript = self._create_transcript(feature, transcript_accession, gene_data)
                self.transcript_data_by_accession[transcript_accession] = transcript
            else:
                if feature.iv.chrom != transcript[CONTIG]:
                    GFFParser._store_other_chrom(transcript, feature)

            # No need to store chrom/strand for each feature, will use transcript
            if feature.type in self.GTF_TRANSCRIPTS_DATA:
                self._add_transcript_data(transcript_accession, transcript, feature)

            biotype = feature.attr.get("gene_biotype")
            if biotype is None:
                # Ensembl GTFs store biotype info under gene_type or transcript_type
                biotype = feature.attr.get("gene_type")
                if biotype is None:
                    biotype = self._get_biotype_from_transcript_accession(transcript_accession)

            if biotype:
                gene_data["biotype"].add(biotype)
                transcript["biotype"].add(biotype)


class GFF3Parser(GFFParser):
    """ GFF3 - @see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        Used by RefSeq and later Ensembl (82 onwards)

        GFF3 support arbitrary hierarchy

    """

    GFF3_GENES = {"gene", "pseudogene"}
    GFF3_TRANSCRIPTS_DATA = {"exon", "CDS", "cDNA_match", "five_prime_UTR", "three_prime_UTR"}

    def __init__(self, *args, **kwargs):
        super(GFF3Parser, self).__init__(*args, **kwargs)
        self.gene_accession_by_feature_id = defaultdict()
        self.transcript_accession_by_feature_id = defaultdict()
        self.hgnc_pattern = re.compile(r".*\[Source:HGNC.*Acc:(HGNC:)?(\d+)\]")

    def handle_feature(self, feature):
        parent_id = feature.attr.get("Parent")
        # Genes never have parents
        # RefSeq genes are always one of GFF3_GENES, Ensembl has lots of different types (lincRNA_gene etc)
        # Ensembl treats pseudogene as a transcript (has parent)
        if parent_id is None and (feature.type in self.GFF3_GENES or "gene_id" in feature.attr):
            gene_accession = self._get_gene_accession(feature)
            dbxref = self._get_dbxref(feature)
            if not gene_accession:
                gene_accession = dbxref.get("GeneID")  # RefSeq no versions
                if not gene_accession:
                    raise ValueError("Could not obtain 'gene_id', even using 'Dbxref[GeneID]'")

            # Gene can have multiple loci, thus entries in GFF, keep original so all transcripts are added
            gene_data = self.gene_data_by_accession.get(gene_accession)
            if gene_data is None:
                gene_data = self._create_gene(feature, gene_accession)
                # If a gene already exists - then need to merge it...
                self.gene_data_by_accession[gene_accession] = gene_data

            # RefSeq stores HGNC in dbxref
            if hgnc := dbxref.get("HGNC"):
                # Might have HGNC: (5 characters) at start of it
                if hgnc.startswith("HGNC"):
                    hgnc = hgnc[5:]
                gene_data["hgnc"] = hgnc
            elif description := gene_data.get("description"):
                # Ensembl stores HGNC in description comment. Sometimes has "HGNC:" prefix eg:
                # Can be either "[Source:HGNC Symbol%3BAcc:HGNC:8907]" or "[Source:HGNC Symbol%3BAcc:37102]"
                if m := self.hgnc_pattern.match(description):
                    gene_data["hgnc"] = m.group(2)

            self.gene_accession_by_feature_id[feature.attr["ID"]] = gene_accession
        else:
            if feature.type in self.GFF3_TRANSCRIPTS_DATA:
                if feature.type == 'cDNA_match':
                    target = feature.attr["Target"]
                    transcript_accession = target.split()[0]
                else:
                    # Some exons etc may be for miRNAs that have no transcript ID, so skip those (won't have parent)
                    if parent_id:
                        transcript_accession = self.transcript_accession_by_feature_id.get(parent_id)
                    else:
                        logging.warning("Transcript data has no parent: %s" % feature.get_gff_line())
                        transcript_accession = None

                if transcript_accession:
                    transcript = self.transcript_data_by_accession[transcript_accession]
                    self._handle_transcript_data(transcript_accession, transcript, feature)
            else:
                # There are so many different transcript ontology terms just taking everything that
                # has a transcript_id and is child of gene (ie skip miRNA etc that is child of primary_transcript)
                if transcript_accession := self._get_transcript_accession(feature, version_key="version"):
                    assert parent_id is not None
                    gene_accession = self.gene_accession_by_feature_id.get(parent_id)
                    if not gene_accession:
                        raise ValueError("Don't know how to handle feature type %s (not child of gene)" % feature.type)
                    gene_data = self.gene_data_by_accession[gene_accession]
                    self._handle_transcript(gene_data, transcript_accession, feature)

    @staticmethod
    def _get_dbxref(feature):
        """ RefSeq stores attribute with more keys, eg: 'Dbxref=GeneID:7840,HGNC:HGNC:428,MIM:606844' """
        dbxref = {}
        dbxref_str = feature.attr.get("Dbxref")
        if dbxref_str:
            dbxref = dict(d.split(":", 1) for d in dbxref_str.split(","))
        return dbxref

    def _handle_transcript(self, gene_data, transcript_accession, feature):
        """ Sometimes we can get multiple transcripts in the same file - just taking 1st """
        if transcript_accession not in self.transcript_data_by_accession:
            # print("_handle_transcript(%s, %s)" % (gene, feature))
            transcript_data = self._create_transcript(feature, transcript_accession, gene_data)
            if biotype := self._get_biotype_from_transcript_accession(transcript_accession):
                gene_data["biotype"].add(biotype)
                transcript_data["biotype"].add(biotype)

            if feature.attr.get("partial"):
                transcript_data["partial"] = 1
            self.transcript_data_by_accession[transcript_accession] = transcript_data
        self.transcript_accession_by_feature_id[feature.attr["ID"]] = transcript_accession

    def _handle_transcript_data(self, transcript_accession, transcript, feature):
        self._add_transcript_data(transcript_accession, transcript, feature)
