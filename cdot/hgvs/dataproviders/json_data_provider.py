import abc
import gzip
import requests

from collections import defaultdict
from lazy import lazy
from hgvs.dataproviders.interface import Interface
from hgvs.dataproviders.seqfetcher import SeqFetcher
from hgvs.exceptions import HGVSDataNotAvailableError
from intervaltree import IntervalTree
from typing import List

from bioutils.assemblies import make_ac_name_map, make_name_ac_map

from cdot import get_data_schema_int, __version__
from cdot import models

def get_ac_name_map(assembly_name):
    if assembly_name == "GRCh37":
        assembly_name = 'GRCh37.p13'  # Original build didn't have MT
    return make_ac_name_map(assembly_name)

def get_name_ac_map(assembly_name):
    if assembly_name == "GRCh37":
        assembly_name = 'GRCh37.p13'  # Original build didn't have MT
    return make_name_ac_map(assembly_name)



class AbstractJSONDataProvider(Interface):
    # All cdot data is 'splign', it's the method used in NCBI/Ensembl GTFs, and we also only pull out 'splign' from UTA
    NCBI_ALN_METHOD = "splign"
    required_version = "1.1"

    def __init__(self, assemblies: List[str] = None, mode=None, cache=None, seqfetcher=None):
        """ assemblies: defaults to ["GRCh37", "GRCh38"]
            seqfetcher defaults to biocommons SeqFetcher()
        """
        if assemblies is None:
            assemblies = ["GRCh37", "GRCh38"]

        super().__init__(mode=mode, cache=cache)
        if seqfetcher:
            try:
                seqfetcher.set_data_provider(self)
            except AttributeError:
                pass
        else:
            seqfetcher = SeqFetcher()
        self.seqfetcher = seqfetcher
        self.assembly_maps = {}
        for assembly_name in assemblies:
            self.assembly_maps[assembly_name] = get_ac_name_map(assembly_name)
        self.assembly_by_contig = {}
        for assembly_name, contig_map in self.assembly_maps.items():
            self.assembly_by_contig.update({contig: assembly_name for contig in contig_map.keys()})


    @abc.abstractmethod
    def _get_transcript(self, tx_ac):
        pass

    @abc.abstractmethod
    def _get_gene(self, gene):
        pass

    def get_tx_versions(self, accession: str) -> list[int]:
        """ Return every stored integer version for a versionless transcript accession,
            ascending (eg "NM_000059" -> [3, 4]). Empty list if the accession is unknown.

            Used by cdot.hgvs.gene_hgvs.resolve_transcript_version() to substitute an
            adjacent version when the requested one isn't available. Subclasses that can
            enumerate their transcripts implement this; the default raises NotImplementedError
            so providers that can't (and callers who don't ask for version fallback) are
            unaffected.
        """
        raise NotImplementedError("This data provider cannot enumerate transcript versions")

    def _get_transcript_coordinates_for_contig(self, transcript, alt_ac):
        assembly = self.assembly_by_contig.get(alt_ac)
        if assembly is None:
            supported_assemblies = ", ".join(self.assembly_maps.keys())
            raise ValueError(f"Contig '{alt_ac}' not supported. Supported assemblies: {supported_assemblies}")

        return transcript["genome_builds"].get(assembly)

    @staticmethod
    def _get_contig_start_end_strand(build_data):
        contig = build_data["contig"]
        strand = 1 if build_data["strand"] == "+" else -1
        start = build_data["exons"][0][0]
        end = build_data["exons"][-1][1]
        return contig, start, end, strand

    def data_version(self):
        return self.required_version

    def schema_version(self):
        return self.required_version

    def get_assembly_map(self, assembly_name):
        """return a list of accessions for the specified assembly name (e.g., GRCh38.p5) """
        assembly_map = self.assembly_maps.get(assembly_name)
        if assembly_map is None:
            supported_assemblies = ", ".join(self.assembly_maps.keys())
            raise ValueError(f"Assembly '{assembly_name}' not supported. Supported assemblies: {supported_assemblies}")

        return assembly_map

    def sequence_source(self):
        return self.seqfetcher.source

    def get_seq(self, ac, start_i=None, end_i=None):
        return self.seqfetcher.fetch_seq(ac, start_i, end_i)

    @staticmethod
    def _convert_gap_to_cigar(gap):
        """
                gap = 'M196 I1 M61 I1 M181'
                CIGAR = '194=1D60=1D184='
        """

        # This has to/from sequences inverted, so insertion is a deletion
        OP_CONVERSION = {
            "M": "=",
            "I": "D",
            "D": "I",
        }

        cigar_ops = []
        for gap_op in gap.split():
            gap_code = gap_op[0]
            length = int(gap_op[1:])

            cigar_ops.append(str(length) + OP_CONVERSION[gap_code])

        return "".join(cigar_ops)

    def get_tx_exons(self, tx_ac, alt_ac, alt_aln_method):
        self._check_alt_aln_method(alt_aln_method)
        transcript = self._get_transcript(tx_ac)
        if not transcript:
            return None

        assembly_coordinates = self._get_transcript_coordinates_for_contig(transcript, alt_ac)
        if assembly_coordinates is None:
            return None

        tx_exons = []  # Genomic order
        exons = assembly_coordinates["exons"]
        alt_strand = 1 if assembly_coordinates["strand"] == "+" else -1

        for (alt_start_i, alt_end_i, exon_id, cds_start_i, cds_end_i, gap) in exons:
            # cdot tx_start_i/tx_end_i is 1 based while UTA is 0 based
            tx_start_i = cds_start_i - 1
            tx_end_i = cds_end_i
            if gap is not None:
                cigar = self._convert_gap_to_cigar(gap)
            else:
                length = alt_end_i - alt_start_i  # Will be same for both transcript/genomic
                cigar = str(length) + "="

            exon_data = {
                'tx_ac': tx_ac,
                'alt_ac': alt_ac,
                'alt_strand': alt_strand,
                'alt_aln_method': alt_aln_method,
                'ord': exon_id,
                'tx_start_i': tx_start_i,
                'tx_end_i': tx_end_i,
                'alt_start_i': alt_start_i,
                'alt_end_i': alt_end_i,
                'cigar': cigar,
            }
            tx_exons.append(exon_data)

        return tx_exons

    def get_tx_identity_info(self, tx_ac):
        # Get any transcript as it's assembly independent
        transcript = self._get_transcript(tx_ac)
        if not transcript:
            return None

        tx_info = self._get_transcript_info(transcript)

        # Only using lengths (same in each build) not coordinates so grab anything
        exons = []
        for build_coordinates in transcript["genome_builds"].values():
            exons = build_coordinates["exons"]
            break

        stranded_order_exons = sorted(exons, key=lambda e: e[2])  # sort by exon_id
        tx_info["lengths"] = [ex[4] + 1 - ex[3] for ex in stranded_order_exons]
        tx_info["tx_ac"] = tx_ac
        tx_info["alt_ac"] = tx_ac  # Same again
        tx_info["alt_aln_method"] = "transcript"
        return tx_info

    @staticmethod
    def _get_transcript_info(transcript):
        gene_name = transcript["gene_name"]
        cds_start_i = transcript.get("start_codon")
        cds_end_i = transcript.get("stop_codon")
        return {
            "hgnc": gene_name,
            "cds_start_i": cds_start_i,
            "cds_end_i": cds_end_i,
        }

    def get_tx_info(self, tx_ac, alt_ac, alt_aln_method):
        self._check_alt_aln_method(alt_aln_method)

        if transcript := self._get_transcript(tx_ac):
            for build_data in transcript["genome_builds"].values():
                if alt_ac == build_data["contig"]:
                    tx_info = self._get_transcript_info(transcript)
                    tx_info["tx_ac"] = tx_ac
                    tx_info["alt_ac"] = alt_ac
                    tx_info["alt_aln_method"] = self.NCBI_ALN_METHOD
                    return tx_info

        raise HGVSDataNotAvailableError(
            f"No tx_info for (tx_ac={tx_ac},alt_ac={alt_ac},alt_aln_method={alt_aln_method})"
        )

    def get_tx_mapping_options(self, tx_ac):
        mapping_options = []
        if transcript := self._get_transcript(tx_ac):
            for build_coordinates in transcript["genome_builds"].values():
                mo = {
                    "tx_ac": tx_ac,
                    "alt_ac": build_coordinates["contig"],
                    "alt_aln_method": self.NCBI_ALN_METHOD,
                }
                mapping_options.append(mo)
        return mapping_options

    def get_acs_for_protein_seq(self, seq):
        """
            This is not implemented. The only caller has comment: 'TODO: drop get_acs_for_protein_seq'
            And is only ever called as a backup when get_pro_ac_for_tx_ac fails
        """
        return []

    def get_gene_info(self, gene):
        gene_info = None
        if g := self._get_gene(gene):
            # UTA produces aliases that look like '{DCML,IMD21,MONOMAC,NFE1B}'
            uta_style_aliases = '{' + g["aliases"].replace(" ", "") + '}'
            gene_info = {
                "hgnc": g["gene_symbol"],
                "maploc": g["map_location"],
                "descr": g["description"],
                "summary": g["summary"],
                "aliases": uta_style_aliases,
                "added": None,  # Don't know where this is stored/comes from (hgnc?)
            }
        return gene_info

    def get_pro_ac_for_tx_ac(self, tx_ac):
        pro_ac = None
        if transcript := self._get_transcript(tx_ac):
            pro_ac = transcript.get("protein")
        return pro_ac

    def get_similar_transcripts(self, tx_ac):
        """ UTA specific functionality that uses tx_similarity_v table
            This is not used by the HGVS library """
        raise NotImplementedError()

    def get_alignments_for_region(self, alt_ac, start_i, end_i, alt_aln_method=None):
        """ Prefer to use get_tx_for_region as this may be removed/deprecated
            This is never called externally, only used to implement get_tx_for_region in UTA data provider. """
        if alt_aln_method is None:
            alt_aln_method = self.NCBI_ALN_METHOD
        return self.get_tx_for_region(alt_ac, alt_aln_method, start_i, end_i)

    def _check_alt_aln_method(self, alt_aln_method):
        if alt_aln_method != self.NCBI_ALN_METHOD:
            raise HGVSDataNotAvailableError(f"cdot only supports alt_aln_method={self.NCBI_ALN_METHOD}")

    @abc.abstractmethod
    def get_tx_for_gene(self, gene):
        pass

    @abc.abstractmethod
    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        pass

    def _validate_schema_compatibility(self, json_schema_version: str):
        """ Raise an error if versions out of sync """
        cdot_client_data_schema_int = get_data_schema_int(__version__)
        cdot_data_schema_version = get_data_schema_int(json_schema_version)
        if cdot_client_data_schema_int < cdot_data_schema_version:
            raise ValueError(f"This cdot client ({__version__}) cannot read {json_schema_version=} - please upgrade.")

class LocalDataProvider(AbstractJSONDataProvider):
    """ For JSON and Redis providers (implemented in cdot_rest)
        https://github.com/SACGF/cdot_rest - cdot_rest.redis_data_provider.RedisDataProvider """

    @abc.abstractmethod
    def _get_transcript_ids_for_gene(self, gene):
        pass

    @abc.abstractmethod
    def _get_contig_interval_tree(self, alt_ac):
        pass

    def get_tx_for_gene(self, gene):
        """ return transcript info records for supplied gene, in order of decreasing length """

        tx_list = []  # Store in tuples with length, so we can sort before returning
        for transcript_id in self._get_transcript_ids_for_gene(gene):
            transcript_data = self._get_transcript(transcript_id)
            cds_start_i = transcript_data.get("start_codon")
            cds_end_i = transcript_data.get("stop_codon")
            for build_data in transcript_data["genome_builds"].values():
                contig, tx_start, tx_end, _ = self._get_contig_start_end_strand(build_data)
                length = tx_end - tx_start
                tx_data = {
                    "hgnc": gene,
                    "cds_start_i": cds_start_i,
                    "cds_end_i": cds_end_i,
                    "tx_ac": transcript_id,
                    "alt_ac": contig,
                    "alt_aln_method": self.NCBI_ALN_METHOD,
                }
                tx_list.append((length, tx_data))

        return [x[1] for x in sorted(tx_list, key=lambda x: x[0], reverse=True)]

    def _get_transcript_tags(self, tx_ac: str, transcript_data: dict, genome_build: str) -> list[str]:
        """
        Return the tag list for a single transcript in the given genome build.

        The default implementation reads the comma-separated ``tag`` field from
        the cdot JSON (e.g. ``"MANE_Select,basic"`` → ``["MANE_Select", "basic"]``).

        Override this method in a subclass to supplement or replace the cdot tag
        data with information from an external source (e.g. a database table that
        stores MANE/canonical status separately from the transcript JSON).  The
        override receives the accession, the full transcript dict, and the build
        name, so it can fall back to the cdot data when the external source has no
        entry. (To answer for many accessions in one query, override
        :meth:`_get_tags_by_tx_ac` instead.)

        Example override (Django)::

            def _get_transcript_tags(self, tx_ac, transcript_data, genome_build):
                tags = super()._get_transcript_tags(tx_ac, transcript_data, genome_build)
                if not tags:
                    # supplement from DB MANE table
                    tags = MyMANEModel.tags_for(tx_ac, genome_build)
                return tags
        """
        build_data = transcript_data["genome_builds"].get(genome_build, {})
        tag_str = build_data.get("tag", "")
        return [t.strip() for t in tag_str.split(",") if t.strip()] if tag_str else []

    def _get_tags_by_tx_ac(self, tx_acs: list[str], genome_build: str) -> dict[str, list[str]]:
        """
        Return {tx_ac: tags} for the given accessions in the given genome build.

        The default loops :meth:`_get_transcript` + :meth:`_get_transcript_tags`
        per accession. Override in a subclass that can answer for all accessions
        in a single query (e.g. a Django provider with a MANE table) to avoid the
        N+1 query pattern on the gene-symbol resolution path.
        """
        result = {}
        for tx_ac in tx_acs:
            transcript_data = self._get_transcript(tx_ac)
            result[tx_ac] = self._get_transcript_tags(tx_ac, transcript_data, genome_build)
        return result

    def get_tx_ac_tags_for_gene(self, gene: str, genome_build: str) -> list[tuple[str, list[str]]]:
        """
        Return [(tx_ac, tags), ...] for all transcripts of the given gene in the
        given genome build, sorted by decreasing transcript length.

        tags is a list of tag strings (e.g. ['MANE_Select', 'basic']).
        Empty list if the transcript has no tags for that build.

        Tag data comes from :meth:`_get_tags_by_tx_ac` (which by default defers to
        the per-transcript :meth:`_get_transcript_tags`); either can be overridden
        in subclasses to supplement the cdot JSON with an external source.
        """
        lengths = []
        for transcript_id in self._get_transcript_ids_for_gene(gene):
            transcript_data = self._get_transcript(transcript_id)
            build_data = transcript_data["genome_builds"].get(genome_build)
            if build_data is None:
                continue
            # Spliced (exonic) transcript length - the sum of exon lengths, NOT
            # the genomic span (which would include introns). Exons are
            # [alt_start, alt_end, ...] so each exon's length is alt_end - alt_start.
            length = sum(exon[1] - exon[0] for exon in build_data["exons"])
            lengths.append((length, transcript_id))

        tags_by_ac = self._get_tags_by_tx_ac([tx_ac for _, tx_ac in lengths], genome_build)
        lengths.sort(key=lambda x: x[0], reverse=True)
        return [(tx_ac, tags_by_ac.get(tx_ac, [])) for _, tx_ac in lengths]

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        """ return transcripts that overlap given region """

        self._check_alt_aln_method(alt_aln_method)

        tx_list = []
        contig_iv_tree = self._get_contig_interval_tree(alt_ac)
        for interval in contig_iv_tree[start_i:end_i+1]:
            transcript_id = interval.data
            transcript_data = self._get_transcript(transcript_id)
            build_data = self._get_transcript_coordinates_for_contig(transcript_data, alt_ac)
            contig, tx_start, tx_end, strand = self._get_contig_start_end_strand(build_data)
            if contig == alt_ac:
                tx_list.append({
                    "alt_ac": alt_ac,
                    "alt_aln_method": self.NCBI_ALN_METHOD,
                    "alt_strand": strand,
                    "start_i": tx_start,
                    "end_i": tx_end,
                    "tx_ac": transcript_id,
                })
        return tx_list

    @staticmethod
    def _get_tx_by_gene_and_intervals(transcript_iter_items):
        # The region query works on exons, but storing all of these makes the interval tree huge
        # So we just store the start/end of each transcript ID, then look up the exons at retrieval time

        tx_by_gene = defaultdict(set)
        tx_intervals = defaultdict(IntervalTree)
        for transcript_id, transcript_data in transcript_iter_items:
            if gene_name := transcript_data['gene_name']:
                tx_by_gene[gene_name].add(transcript_id)

            for build_data in transcript_data["genome_builds"].values():
                contig = build_data["contig"]
                tx_start = build_data["exons"][0][0]
                tx_end = build_data["exons"][-1][1]
                tx_intervals[contig][tx_start:tx_end] = transcript_id
        return tx_by_gene, tx_intervals


class JSONDataProvider(LocalDataProvider):
    """ Local JSON file """
    def __init__(self, file_or_filename_list, mode=None, cache=None, seqfetcher=None):
        assemblies = set()
        self.transcripts = {}
        self.genes = {}
        for file_or_filename in file_or_filename_list:
            if isinstance(file_or_filename, str):
                open_func = gzip.open if file_or_filename.endswith(".gz") else open
                f = open_func(file_or_filename)
            else:
                f = file_or_filename
            # msgspec decodes ~2.5x faster than json.load and into more compact typed
            # structs (see cdot/models.py). The structs are dict-read-compatible, so the
            # rest of the provider is unchanged and consumers see identical output.
            data = models.loads(f.read())
            assemblies.update(data.genome_builds)
            self.transcripts.update(data.transcripts)
            if genes := data.genes:
                for g in genes.values():
                    if gene_symbol := g.gene_symbol:
                        self.genes[gene_symbol] = g
            cdot_data_version_str = data.cdot_version
            self._validate_schema_compatibility(cdot_data_version_str)
            self.cdot_data_version = tuple(int(v) for v in cdot_data_version_str.split("."))

        super().__init__(assemblies=assemblies, mode=mode, cache=cache, seqfetcher=seqfetcher)

    def _get_transcript(self, tx_ac):
        return self.transcripts.get(tx_ac)

    def get_tx_versions(self, accession: str) -> list[int]:
        prefix = accession + "."  # so "NM_000059" doesn't also match "NM_0000591"
        versions = []
        for tx_ac in self.transcripts:
            if tx_ac.startswith(prefix):
                _base, _dot, version = tx_ac.rpartition(".")
                if version.isdigit():
                    versions.append(int(version))
        return sorted(versions)

    def _get_gene(self, gene):
        return self.genes.get(gene)

    def _get_transcript_ids_for_gene(self, gene):
        tx_by_gene, _ = self._tx_by_gene_and_intervals
        return tx_by_gene[gene]

    def _get_contig_interval_tree(self, alt_ac):
        _, tx_intervals = self._tx_by_gene_and_intervals
        return tx_intervals[alt_ac]

    def get_pro_ac_for_tx_ac(self, tx_ac):
        if self.cdot_data_version < (0, 2, 8):
            cdot_version = '.'.join(str(v) for v in self.cdot_data_version)
            msg = f"ProteinID not in your JSON data version '{cdot_version}'. " \
                  "Please use data generated from cdot >= 0.2.8"
            raise NotImplementedError(msg)
        return super().get_pro_ac_for_tx_ac(tx_ac)

    @lazy
    def _tx_by_gene_and_intervals(self):
        return self._get_tx_by_gene_and_intervals(self.transcripts.items())

    def get_gene_info(self, gene):
        if self.cdot_data_version < (0, 2, 10):
            cdot_version = '.'.join(str(v) for v in self.cdot_data_version)
            msg = f"Gene Info not in your JSON data version '{cdot_version}'. " \
                  "Please use data generated from cdot >= 0.2.10"
            raise NotImplementedError(msg)
        return super().get_gene_info(gene)


class RESTDataProvider(AbstractJSONDataProvider):

    def __init__(self, url=None, secure=True, mode=None, cache=None, seqfetcher=None):
        assemblies = ["GRCh37", "GRCh38"]
        super().__init__(assemblies=assemblies, mode=mode, cache=cache, seqfetcher=seqfetcher)
        if url is None:
            if secure:
                url = "https://cdotlib.org"
            else:
                url = "http://cdotlib.org"
        self.url = url
        self.transcripts = {}
        self.genes = {}

    def _get_from_url(self, url):
        data = None
        response = requests.get(url)
        if response.ok:
            if 'application/json' in response.headers.get('Content-Type'):
                data = response.json()
            else:
                raise ValueError("Non-json response received for '%s' - are you behind a firewall?" % url)
        return data

    def _post_to_url(self, url, json_body):
        """ POST helper mirroring _get_from_url.

            Returns the parsed JSON response, or None if the endpoint is absent (404/405) so
            the caller can fall back to an older code path. Any other non-ok status raises.
        """
        response = requests.post(url, json=json_body)
        if response.status_code in (404, 405):
            return None  # Server predates this endpoint - caller falls back
        if not response.ok:
            response.raise_for_status()
        if 'application/json' in response.headers.get('Content-Type', ''):
            return response.json()
        raise ValueError("Non-json response received for '%s' - are you behind a firewall?" % url)

    def _get_transcript(self, tx_ac):
        # We store None for 404 on REST
        if tx_ac in self.transcripts:
            return self.transcripts[tx_ac]

        transcript = self._get_from_url(self.url + "/transcript/" + tx_ac)
        self.transcripts[tx_ac] = transcript
        return transcript

    def get_tx_versions(self, accession: str) -> list[int]:
        """ For a versionless accession (eg "NM_000059") the /transcript/<ac> endpoint returns
            an object containing every stored version, keyed by full accession. We parse the
            versions out and warm the cache with each versioned transcript as a bonus (so a
            following _get_transcript() for the chosen version is a hit).
            See https://cdotlib.org/static/api-docs.html
        """
        versions = []
        if data := self._get_from_url(self.url + "/transcript/" + accession):
            for full_ac, transcript in data.items():
                self.transcripts[full_ac] = transcript  # warm cache
                _base, _dot, version = full_ac.rpartition(".")
                if version.isdigit():
                    versions.append(int(version))
        return sorted(versions)

    def prefetch(self, tx_acs, batch=True, max_workers=10):
        """ Read-ahead cache warming: populate self.transcripts up front.

            The biocommons HGVS Interface is one-transcript-at-a-time, and each c_to_g makes a
            single /transcript/<ac> call that's then reused for get_tx_info/get_tx_exons/etc.
            When you know the transcripts up front (eg bulk HGVS processing) calling this first
            means every later _get_transcript() is a cache hit - no Interface change required.

            By default this POSTs the whole accession list to the batch /transcripts endpoint in
            a single round-trip. Accessions may be versionless (eg "NM_000059"), in which case
            the server expands them to every available version - useful for warming the
            version-bump path. Servers without the batch endpoint fall back to a concurrent
            thread-pool of single /transcript/<ac> requests (set batch=False to force this).

            Already-cached accessions are skipped. Missing/404 transcripts are stored as None,
            matching _get_transcript. Returns the number of transcripts added to the cache.
        """
        to_fetch = {ac for ac in tx_acs if ac not in self.transcripts}
        if not to_fetch:
            return 0

        if batch:
            data = self._post_to_url(self.url + "/transcripts", {"ids": sorted(to_fetch)})
            if data is not None:
                # Keyed by full accession: versionless ids expand to all their versions,
                # versioned misses come back as null. Drop straight into the cache.
                for tx_ac, transcript in data.items():
                    self.transcripts[tx_ac] = transcript
                return len(data)
            # else: no batch endpoint on this server - fall through to concurrent singles

        return self._prefetch_concurrent(to_fetch, max_workers=max_workers)

    def _prefetch_concurrent(self, to_fetch, max_workers=10):
        """ Fallback pre-warm: fire one /transcript/<ac> per accession concurrently. Turns N
            sequential round-trips into ceil(N/max_workers) waves. Used when the server has no
            batch endpoint. Versionless accessions aren't expanded here - pass exact versions.
        """
        from concurrent.futures import ThreadPoolExecutor

        def _fetch(tx_ac):
            try:
                return tx_ac, self._get_from_url(self.url + "/transcript/" + tx_ac)
            except Exception:
                # Match _get_transcript's "store None for missing" behaviour; a failed
                # prefetch must not be fatal - a later _get_transcript can retry/raise.
                return tx_ac, None

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for tx_ac, transcript in executor.map(_fetch, to_fetch):
                self.transcripts[tx_ac] = transcript
        return len(to_fetch)

    def prefetch_from_hgvs(self, hgvs_strings, clean=True, **kwargs):
        """ Convenience wrapper: extract transcript accessions from HGVS strings and prefetch.

            Runs clean_hgvs() first (pure string, no provider) so messy input still yields a
            usable accession, then warms the cache for every distinct transcript accession.
            Non-transcript references (gene-only HGVS, genomic NC_ contigs) are skipped.

            Versions are kept as-is, so "NM_000059:c.1A>G" (no version) prefetches every
            available version in one round-trip - the input the version-bump path needs.

            Extra kwargs (batch, max_workers) are passed through to prefetch().
            Returns the number of transcripts added to the cache.
        """
        from cdot.hgvs.clean import clean_hgvs

        accessions = set()
        for hgvs_string in hgvs_strings:
            if clean:
                hgvs_string, _fixes = clean_hgvs(hgvs_string)
            if accession := self._accession_from_hgvs(hgvs_string):
                accessions.add(accession)
        return self.prefetch(accessions, **kwargs)

    @staticmethod
    def _accession_from_hgvs(hgvs_string):
        """ Pull the transcript accession (before the ':') from an HGVS string, or None.

            Strips a gene-in-parens suffix ("NM_000059.4(BRCA2)" -> "NM_000059.4") and returns
            None for anything that isn't a transcript (gene symbols, genomic NC_ refs).
        """
        from cdot.hgvs.clean import _looks_like_transcript

        accession, sep, _allele = hgvs_string.partition(":")
        if not sep:
            return None
        accession = accession.strip().split("(", 1)[0]  # drop "(GENE)" suffix if present
        if _looks_like_transcript(accession):
            return accession
        return None

    def _get_gene(self, gene_name):
        # We store None for 404 on REST
        if gene_name in self.genes:
            return self.genes[gene_name]

        gene = self._get_from_url(self.url + "/gene/" + gene_name)
        self.genes[gene_name] = gene
        return gene

    def get_tx_for_gene(self, gene_name):
        tx_list = []
        if data := self._get_from_url(f"{self.url}/transcripts/gene/{gene_name}"):
            tx_list = data["results"]
        return tx_list

    def get_tx_ac_tags_for_gene(self, gene: str, genome_build: str) -> list[tuple[str, list[str]]]:
        """ Return [(tx_ac, tags), ...] for a gene/build, sorted longest-first.

            Used by cdot.hgvs.gene_hgvs.resolve_gene_hgvs() to pick the best
            (MANE/canonical) transcript. The server returns the same data as
            LocalDataProvider.get_tx_ac_tags_for_gene(); ranking stays client-side.
        """
        tx_list = []
        url = f"{self.url}/transcripts/gene/{gene}/tags/{genome_build}"
        if data := self._get_from_url(url):
            tx_list = [(tx_ac, tags) for tx_ac, tags in data["results"]]
        return tx_list

    def get_tx_for_region(self, alt_ac, alt_aln_method, start_i, end_i):
        self._check_alt_aln_method(alt_aln_method)

        tx_list = []
        url = f"{self.url}/transcripts/region/{alt_ac}/{alt_aln_method}/{start_i}/{end_i}"
        if data := self._get_from_url(url):
            tx_list = data["results"]
        return tx_list
