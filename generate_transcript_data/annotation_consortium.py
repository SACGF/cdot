""" Where RefSeq and Ensembl differ in how they annotate, independent of whether the file is GTF or GFF3

    The GTF/GFF3 parsers (gff_parser.py) work out the gene/transcript/exon hierarchy, then ask the
    consortium for anything written differently: HGNC IDs, CCDS, protein versions, translation
    exceptions, and RefSeq's transcript-less mitochondrial genes. The plain AnnotationConsortium
    is used for files from neither (eg UCSC GTFs), which only have the standard columns.

    Issue #101 """
import gzip
import logging
import re
from typing import Optional
from urllib.parse import unquote

from generate_transcript_data.transcript_builder import TranscriptBuilder


class AnnotationConsortium:
    name = "generic"

    def get_gene_accession_fallback(self, feature) -> Optional[str]:
        """ Gene accession when the row has no gene_id """
        return None

    def get_hgnc(self, feature, gene_data) -> Optional[str]:
        """ HGNC ID (without 'HGNC:' prefix) from a gene row """
        return None

    def get_ccds(self, feature) -> Optional[str]:
        return None

    def get_transcript_support_level(self, feature) -> Optional[str]:
        return None

    def get_protein_accession(self, feature) -> Optional[str]:
        """ Versioned protein accession from a CDS row, or None for other rows """
        if feature.type == "CDS":
            if protein_accession := feature.attr.get("protein_id"):
                if protein_version := feature.attr.get("protein_version"):
                    protein_accession = f"{protein_accession}.{protein_version}"
                if "." not in protein_accession:
                    raise ValueError(f"Protein '{protein_accession}' missing version")
                return protein_accession
        return None

    def add_feature_annotations(self, transcript_accession: str, feature, builder: TranscriptBuilder):
        """ Called for each exon/CDS etc row of a transcript (on the transcript's contig) """
        pass

    def add_transcript_annotations(self, transcript_accession: str, feature, builder: TranscriptBuilder):
        """ Called for the row that creates a GFF3 transcript (eg mRNA) """
        pass

    def handle_transcript_feature(self, transcript_accession: str, feature, builder: TranscriptBuilder) -> bool:
        """ A row of a transcript that isn't a standard exon/CDS/codon. Returns whether it was handled """
        return False

    def get_fake_transcript_accession_and_parent(self, feature, builder: TranscriptBuilder) \
            -> Optional[tuple[str, Optional[str]]]:
        """ For rows that belong to a transcript but have no transcript_id """
        return None


class RefSeq(AnnotationConsortium):
    """ RefSeq GFF3: IDs in Dbxref, exceptions and transl_except on CDS rows, and the mitochondrial
        proteins (YP_) have no NM_ transcript so we make fake ones """
    name = "refseq"

    EXCEPTION_GENOME_MISMATCH = "annotated by transcript or proteomic data"
    EXCEPTION_RIBOSOMAL_SLIPPAGE = "ribosomal slippage"
    MITO_CONTIG = "NC_012920.1"
    # Note=The RefSeq transcript has 1 substitution%2C 1 non-frameshifting indel compared to this genomic sequence
    GENOME_MISMATCH_NOTE = re.compile(r"The RefSeq (transcript|protein) has (.+?) compared to this genomic sequence")
    GENOME_MISMATCH_COVERAGE = re.compile(r"aligns at (\d+(?:\.\d+)?)% coverage")
    # cDNA_match attributes describing how well the transcript aligns to the genome
    ALIGNMENT_STATS = ["num_mismatch", "gap_count", "pct_identity_gap", "pct_coverage"]

    @staticmethod
    def _get_dbxref(feature):
        """ RefSeq stores attribute with more keys, eg: 'Dbxref=GeneID:7840,HGNC:HGNC:428,MIM:606844' """
        dbxref = {}
        dbxref_str = feature.attr.get("Dbxref")
        if dbxref_str:
            dbxref = dict(d.split(":", 1) for d in dbxref_str.split(","))
        return dbxref

    def get_gene_accession_fallback(self, feature) -> Optional[str]:
        return self._get_dbxref(feature).get("GeneID")  # RefSeq no versions

    def get_hgnc(self, feature, gene_data) -> Optional[str]:
        if hgnc := self._get_dbxref(feature).get("HGNC"):
            # Might have HGNC: (5 characters) at start of it
            if hgnc.startswith("HGNC:"):
                hgnc = hgnc[5:]
        return hgnc

    def get_ccds(self, feature) -> Optional[str]:
        return self._get_dbxref(feature).get("CCDS")

    @staticmethod
    def _parse_exception(exception):
        """ RefSeq 'exception' attribute, eg 'alternative start codon%2C annotated by transcript or proteomic data' """
        if not exception:
            return []
        return [e.strip() for e in exception.replace("%2C", ",").split(",") if e.strip()]

    @staticmethod
    def _parse_transl_except(transl_except):
        """ RefSeq, eg '(pos:complement(31105951..31105953)%2Caa:Sec),(pos:4261..4262%2Caa:TERM)'
            Yields (start, end, amino_acid) with 0-based half-open genomic coordinates """
        for pos, amino_acid in re.findall(r"pos:(.+?)(?:%2C|,)aa:(\w+)", transl_except):
            coordinates = [int(c) for c in re.findall(r"\d+", pos)]
            yield min(coordinates) - 1, max(coordinates), amino_acid

    @staticmethod
    def _number(value: str):
        """ '100' -> 100, '99.962' -> 99.962 """
        number = float(value)
        return int(number) if number.is_integer() else number

    @classmethod
    def _parse_genome_mismatch_note(cls, note) -> dict[str, dict]:
        """ RefSeq Note, eg 'The RefSeq transcript has 5 substitutions%2C 5 non-frameshifting indels compared to
            this genomic sequence' -> {"transcript": {"substitutions": 5, "non_frameshifting_indels": 5}}
            A Note can have other sentences, and one for the protein as well """
        mismatches = {}
        for sequence_type, differences in cls.GENOME_MISMATCH_NOTE.findall(unquote(note)):
            counts = {}
            for item in re.split(r",\s*|\s+and\s+", differences):
                if m := re.fullmatch(r"(\d+) (.+)", item.strip()):
                    kind = re.sub(r"[^a-z0-9]+", "_", m.group(2).lower()).strip("_")
                    if not kind.endswith("s"):
                        kind += "s"
                    counts[kind] = int(m.group(1))
                elif m := cls.GENOME_MISMATCH_COVERAGE.fullmatch(item.strip()):
                    counts["pct_coverage"] = cls._number(m.group(1))
                else:
                    logging.warning("Unknown genome mismatch '%s' in Note: %s", item, note)
            if counts:
                mismatches[sequence_type] = counts
        return mismatches

    def _add_genome_mismatch_note(self, transcript_accession, feature, builder):
        if note := feature.attr.get("Note"):
            for sequence_type, counts in self._parse_genome_mismatch_note(note).items():
                builder.set_genome_mismatch(transcript_accession, sequence_type, counts)

    def _add_alignment_stats(self, transcript_accession, feature, builder):
        """ cDNA_match rows repeat the stats for the whole alignment. Only kept if it isn't perfect """
        stats = {}
        for key in self.ALIGNMENT_STATS:
            if (value := feature.attr.get(key)) is not None:
                stats[key] = self._number(value)
        imperfect = stats.get("num_mismatch", 0) > 0 or stats.get("gap_count", 0) > 0 \
            or stats.get("pct_identity_gap", 100) < 100 or stats.get("pct_coverage", 100) < 100
        if imperfect:
            builder.set_genome_mismatch(transcript_accession, "alignment", stats)

    def add_transcript_annotations(self, transcript_accession, feature, builder):
        self._add_genome_mismatch_note(transcript_accession, feature, builder)

    def add_feature_annotations(self, transcript_accession, feature, builder):
        self._add_genome_mismatch_note(transcript_accession, feature, builder)
        if feature.type == "cDNA_match":
            self._add_alignment_stats(transcript_accession, feature, builder)
        exceptions = self._parse_exception(feature.attr.get("exception"))
        if feature.type == "CDS":
            # RefSeq repeats these on every CDS row
            if transl_except := feature.attr.get("transl_except"):
                for start, end, amino_acid in self._parse_transl_except(transl_except):
                    builder.add_transl_except(transcript_accession, start, end, amino_acid)
            if transl_table := feature.attr.get("transl_table"):
                builder.set_transl_table(transcript_accession, int(transl_table))
            if self.EXCEPTION_RIBOSOMAL_SLIPPAGE in exceptions:
                builder.add_ribosomal_slippage_cds(transcript_accession, feature)
            for exception in exceptions:
                if exception == self.EXCEPTION_GENOME_MISMATCH:
                    builder.add_genome_mismatch_exception(transcript_accession, exception)
                else:
                    builder.add_translation_exception(transcript_accession, exception)
        else:
            # RefSeq exon rows: the transcript differs from this genome
            for exception in exceptions:
                builder.add_genome_mismatch_exception(transcript_accession, exception)

    def get_fake_transcript_accession_and_parent(self, feature, builder):
        # In RefSeq there are no transcript_ids for MT genes/mRNAs
        # The proteins have "YP_" prefix (no corresp. NM_ transcript, so we will create fake ones
        # Some GFFs have 'mRNA' features, however some only have gene and CDS
        if feature.iv.chrom == self.MITO_CONTIG:
            if feature.type == "mRNA":
                transcript_accession = "fake-" + feature.attr.get("ID")
                return transcript_accession, None
            elif feature.type == "CDS":
                GENE_PREFIX = "gene-"
                RNA_PREFIX = "rna-"
                if parent := feature.attr.get("Parent"):
                    if parent.startswith(GENE_PREFIX):
                        rna_parent = parent.replace(GENE_PREFIX, RNA_PREFIX)
                        transcript_accession = "fake-" + rna_parent
                        # Chuck a fake exon in there
                        builder.add_exon_coordinates(transcript_accession, feature)
                        return transcript_accession, parent
        return None


class Ensembl(AnnotationConsortium):
    """ Ensembl GTF: protein versions in their own attribute, selenocysteines as their own rows,
        HGNC only in the gene description (GENCODE metadata is used instead, see cdot_json.py) """
    name = "ensembl"

    MITO_CONTIG = "MT"
    # Ensembl keeps the genetic code in its database (seq_region attribute 'codon_table'), not in the GTF/GFF3
    MITO_TRANSL_TABLE = 2  # Vertebrate mitochondrial

    # Can be either "[Source:HGNC Symbol%3BAcc:HGNC:8907]" or "[Source:HGNC Symbol%3BAcc:37102]"
    HGNC_PATTERN = re.compile(r".*\[Source:HGNC.*Acc:(HGNC:)?(\d+)]")

    def get_hgnc(self, feature, gene_data) -> Optional[str]:
        if description := gene_data.get("description"):
            if m := self.HGNC_PATTERN.match(description):
                return m.group(2)
        return None

    def get_ccds(self, feature) -> Optional[str]:
        return feature.attr.get("ccds_id")

    def get_transcript_support_level(self, feature) -> Optional[str]:
        return feature.attr.get("transcript_support_level")

    def get_protein_accession(self, feature) -> Optional[str]:
        # Ensembl GTF: CDS protein_id = ENSP00000477624 protein_version = 1
        # Ensembl GFF3 (release 114+): CDS protein_id = ENSP00000477624 version = 1
        # Ensembl GFF3 (release 113 and earlier): CDS protein_id = ENSP00000477624  <----- no version, can't use
        if feature.type == "CDS":
            if protein_id := feature.attr.get("protein_id"):
                protein_version = feature.attr.get("protein_version") or feature.attr.get("version")
                if not protein_version:
                    raise ValueError(f"Protein '{protein_id}' missing version. Ensembl GFF3 files before release "
                                     f"114 do not carry protein versions, use the Ensembl GTF instead "
                                     f"(https://ftp.ensembl.org/pub/release-<N>/gtf/)")
                return f"{protein_id}.{protein_version}"
        return None

    def add_feature_annotations(self, transcript_accession, feature, builder):
        if feature.type == "CDS" and feature.iv.chrom == self.MITO_CONTIG:
            builder.set_transl_table(transcript_accession, self.MITO_TRANSL_TABLE)

    def handle_transcript_feature(self, transcript_accession, feature, builder) -> bool:
        if feature.type == "Selenocysteine":
            builder.add_transl_except(transcript_accession, feature.iv.start, feature.iv.end, "Sec")
            return True
        return False


CONSORTIA = {c.name: c for c in [AnnotationConsortium, RefSeq, Ensembl]}


def get_annotation_consortium(name: str) -> AnnotationConsortium:
    try:
        return CONSORTIA[name.lower()]()
    except KeyError:
        raise ValueError(f"Unknown annotation consortium '{name}', expected one of {sorted(CONSORTIA)}")


_ENSEMBL_ACCESSION = re.compile(r"\bENS[A-Z]*[GT]\d{11}")
_REFSEQ_GENE_ID = re.compile(r"GeneID:\d+")
_DETECT_MAX_LINES = 1000


def detect_annotation_consortium(filename: str) -> str:
    """ Sniff the first rows of a GTF/GFF3: Ensembl uses ENSG/ENST accessions, RefSeq puts
        'GeneID:' in Dbxref (GFF3) or db_xref (GTF). Anything else is generic """
    opener = gzip.open if filename.endswith(".gz") else open
    name = AnnotationConsortium.name
    with opener(filename, "rt") as f:
        for i, line in enumerate(f):
            if i >= _DETECT_MAX_LINES:
                break
            if line.startswith("#"):
                continue
            if _ENSEMBL_ACCESSION.search(line):
                name = Ensembl.name
                break
            if _REFSEQ_GENE_ID.search(line):
                name = RefSeq.name
                break
    logging.info("Detected '%s' annotation conventions in '%s'", name, filename)
    return name
