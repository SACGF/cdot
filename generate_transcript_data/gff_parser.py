""" GTF / GFF3 readers: work out the gene -> transcript -> exon/CDS hierarchy from each file format
    and feed it to a TranscriptBuilder, asking the AnnotationConsortium (RefSeq/Ensembl) for anything
    the two consortia write differently. See issue #101 """
import abc
import logging
from bioutils.assemblies import make_name_ac_map
from collections import Counter
from typing import Optional

import HTSeq

from generate_transcript_data.annotation_consortium import get_annotation_consortium, detect_annotation_consortium
from generate_transcript_data.transcript_builder import TranscriptBuilder, CONTIG, CODING_FEATURES

EXCLUDE_BIOTYPES = {"transcript"}  # feature.type we won't put into biotype


def get_name_ac_map(assembly_name):
    if assembly_name == "GRCh37":
        assembly_name = 'GRCh37.p13'  # Original build didn't have MT
    return make_name_ac_map(assembly_name)


class GFFParser(abc.ABC):
    FEATURE_ALLOW_LIST = {}
    FEATURE_IGNORE_LIST = {"biological_region", "chromosome", "region", "scaffold", "supercontig"}

    def __init__(self, filename, genome_build, url,
                 discard_contigs_with_underscores=True, no_contig_conversion=False,
                 skip_missing_parents=False, annotation_consortium: Optional[str] = None):
        """ annotation_consortium: 'refseq', 'ensembl' or 'generic', detected from the file if not given """
        self.filename = filename
        self.genome_build = genome_build
        self.url = url
        self.discard_contigs_with_underscores = discard_contigs_with_underscores
        self.skip_missing_parents = skip_missing_parents

        self.discarded_contigs = Counter()
        self.skipped_features_no_parents = Counter()  # if skip_missing_parents

        name_ac_map = {}
        if not no_contig_conversion:
            try:
                name_ac_map = get_name_ac_map(genome_build)
            except FileNotFoundError as e:
                raise FileNotFoundError(f"Your genome build '{genome_build}' doesn't have an assembly in biocommons "
                                        f"bioutils. You need this for Biocommons HGVS conversion but if you are using "
                                        f"this for another purpose, you can try adding --no-contig-conversion") from e
        self.name_ac_map = name_ac_map

        if annotation_consortium is None:
            annotation_consortium = detect_annotation_consortium(filename)
        self.consortium = get_annotation_consortium(annotation_consortium)
        self.builder = TranscriptBuilder(genome_build, url, name_ac_map)

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

    def get_genes_and_transcripts(self):
        self._parse()
        genes, transcripts = self.builder.finish()

        if self.discarded_contigs:
            print("Discarded contigs: %s" % self.discarded_contigs)

        if self.skipped_features_no_parents:
            print("Skipped features w/o parents: %s" % self.skipped_features_no_parents)

        return genes, transcripts

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
            # GFF3 gene rows put it in 'version', GTF gene rows (and everything else) in 'gene_version'
            gene_version = feature.attr.get("gene_version")
            if feature.type == "gene":
                gene_version = feature.attr.get("version") or gene_version

            if gene_version:
                gene_accession = f"{gene_id}.{gene_version}"
            else:
                gene_accession = gene_id
        return gene_accession

    def _add_transcript_data(self, transcript_accession, transcript, feature):
        """ An exon/CDS/codon/cDNA_match row of a transcript """
        if self.builder.add_feature(transcript_accession, transcript, feature):
            self.consortium.add_feature_annotations(transcript_accession, feature, self.builder)

    def _handle_protein_version(self, transcript_accession, feature):
        if protein_accession := self.consortium.get_protein_accession(feature):
            self.builder.set_protein(transcript_accession, protein_accession)


class GTFParser(GFFParser):
    """ GTF (GFF2) - used by Ensembl, @see https://gmod.org/wiki/GFF2

        GFF2 only has 2 levels of feature hierarchy, so we have to build or 3 levels of gene/transcript/exons ourselves

        We *have* to use GTF as Ensembl GFF3s don't include the protein version (just the ID)

    """
    GTF_TRANSCRIPTS_DATA = CODING_FEATURES | {"exon"}
    FEATURE_ALLOW_LIST = GTF_TRANSCRIPTS_DATA | {"gene", "transcript", "Selenocysteine"}

    def handle_feature(self, feature):
        gene_accession = self._get_gene_accession(feature)
        if gene_accession is None:
            gene_data = {}  # Empty
            # logging.warning("Read gene accession = None for %s", feature)  # Think this may not happen now with GTFs
        else:
            gene_data = self.builder.get_or_create_gene(feature, gene_accession)

        if transcript_accession := self._get_transcript_accession(feature, version_key="transcript_version"):
            transcript = self.builder.transcripts.get(transcript_accession)
            if transcript is None:
                transcript = self.builder.create_transcript(feature, transcript_accession, gene_data)
            else:
                if feature.iv.chrom != transcript[CONTIG]:
                    self.builder.store_other_chrom(transcript, feature)

            # No need to store chrom/strand for each feature, will use transcript
            if feature.type in self.GTF_TRANSCRIPTS_DATA:
                self._add_transcript_data(transcript_accession, transcript, feature)
                if transcript_support_level := self.consortium.get_transcript_support_level(feature):
                    transcript["transcript_support_level"] = transcript_support_level
                if ccds := self.consortium.get_ccds(feature):
                    transcript["ccds"] = ccds
            else:
                self.consortium.handle_transcript_feature(transcript_accession, feature, self.builder)

            biotype = feature.attr.get("gene_biotype")
            if biotype is None:
                # Ensembl GTFs store biotype info under gene_type or transcript_type
                biotype = feature.attr.get("gene_type")

            if biotype:
                gene_data["biotype"].add(biotype)
                transcript["biotype"].add(biotype)

            if feature.source:
                gene_data["source"].add(feature.source)
                transcript["source"].add(feature.source)

            self._handle_protein_version(transcript_accession, feature)
            self.builder.add_tags(transcript, feature)


class GFF3Parser(GFFParser):
    """ GFF3 - @see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        Used by RefSeq and later Ensembl (82 onwards)

        GFF3 support arbitrary hierarchy

    """

    GFF3_GENES = {"gene", "pseudogene", "ncRNA_gene"}
    GFF3_TRANSCRIPTS_DATA = {"exon", "CDS", "cDNA_match", "five_prime_UTR", "three_prime_UTR"}

    def __init__(self, *args, **kwargs):
        super(GFF3Parser, self).__init__(*args, **kwargs)
        self.gene_accession_by_feature_id = {}
        self.transcript_accession_by_feature_id = {}

    def handle_feature(self, feature):
        parent_id = feature.attr.get("Parent")
        # Genes never have parents
        # RefSeq genes are always one of GFF3_GENES, Ensembl has lots of different types (lincRNA_gene etc)
        # Ensembl treats pseudogene as a transcript (has parent)
        if parent_id is None and (feature.type in self.GFF3_GENES or "gene_id" in feature.attr):
            # Gene
            gene_accession = self._get_gene_accession(feature)
            if not gene_accession:
                gene_accession = self.consortium.get_gene_accession_fallback(feature)
                if not gene_accession:
                    raise ValueError("Could not obtain 'gene_id', even using 'Dbxref[GeneID]'")

            # Gene can have multiple loci, thus entries in GFF, keep original so all transcripts are added
            gene_data = self.builder.get_or_create_gene(feature, gene_accession)

            if hgnc := self.consortium.get_hgnc(feature, gene_data):
                gene_data["hgnc"] = hgnc

            if feature.source:
                gene_data["source"].add(feature.source)

            self.gene_accession_by_feature_id[feature.attr["ID"]] = gene_accession
        else:
            # Transcripts
            transcript_accession = None
            if feature.type in self.GFF3_TRANSCRIPTS_DATA:
                if feature.type == 'cDNA_match':
                    target = feature.attr["Target"]
                    transcript_accession = target.split()[0]
                else:
                    # Some exons etc may be for miRNAs that have no transcript ID, so skip those (won't have parent)
                    if parent_id:
                        transcript_accession = self.transcript_accession_by_feature_id.get(parent_id)
                        if transcript_accession:
                            self._handle_protein_version(transcript_accession, feature)
                    else:
                        logging.warning("Transcript data has no parent: %s" % feature.get_gff_line())

                if transcript_accession:
                    transcript = self.builder.transcripts.get(transcript_accession)
                    if not transcript:
                        msg = f"Couldn't find transcript data for accession '{transcript_accession}'"
                        if self.skip_missing_parents:
                            logging.warning(msg)
                            self.skipped_features_no_parents[feature.type] += 1
                            return
                        raise ValueError(msg)
                    self._gff_handle_transcript_data(transcript_accession, transcript, feature)

            if transcript_accession is None:
                # There are so many different transcript ontology terms just taking everything that
                # has a transcript_id and is child of gene (ie skip miRNA etc that is child of primary_transcript)
                if transcript_and_parent := self._get_transcript_accession_and_parent(feature):
                    transcript_accession, replaced_parent_id = transcript_and_parent
                    if replaced_parent_id:
                        parent_id = replaced_parent_id
                    gene_accession = self.gene_accession_by_feature_id.get(parent_id)
                    if not gene_accession:
                        if self.skip_missing_parents:
                            self.skipped_features_no_parents[feature.type] += 1
                            return
                        msg = f"Don't know how to handle feature type {feature.type} (couldn't find parent gene {parent_id})"
                        raise ValueError(msg)
                    gene_data = self.builder.genes[gene_accession]
                    self._handle_transcript(gene_data, transcript_accession, feature)
                    transcript = self.builder.transcripts[transcript_accession]

                    if feature.type == 'CDS':
                        # We should only be here if transcript_id wasn't on it (chrM)
                        self._gff_handle_transcript_data(transcript_accession, transcript, feature)
                    elif feature.type not in EXCLUDE_BIOTYPES:
                        transcript["biotype"].add(feature.type)

                    if feature.source:
                        transcript["source"].add(feature.source)

    def _get_transcript_accession_and_parent(self, feature) -> Optional[tuple[str, Optional[str]]]:
        if transcript_accession := self._get_transcript_accession(feature, version_key="version"):
            return transcript_accession, None
        # eg RefSeq mitochondrial CDS, which have no transcript
        return self.consortium.get_fake_transcript_accession_and_parent(feature, self.builder)

    def _handle_transcript(self, gene_data, transcript_accession, feature):
        """ Sometimes we can get multiple transcripts in the same file - just taking 1st """
        if transcript_accession not in self.builder.transcripts:
            transcript_data = self.builder.create_transcript(feature, transcript_accession, gene_data)

            if feature.attr.get("partial"):
                transcript_data["partial"] = 1

            self.builder.add_tags(transcript_data, feature)
        self.transcript_accession_by_feature_id[feature.attr["ID"]] = transcript_accession

    def _gff_handle_transcript_data(self, transcript_accession, transcript, feature):
        self._add_transcript_data(transcript_accession, transcript, feature)
        if ccds := self.consortium.get_ccds(feature):
            transcript["ccds"] = ccds
