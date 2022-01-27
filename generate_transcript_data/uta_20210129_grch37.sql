\copy (SELECT transcript.ac, string_agg(distinct transcript.hgnc, ',') as hgnc, string_agg(distinct origin.url, ',') as origin_url,
string_agg(distinct es.alt_ac::varchar, ',') as contig,
string_agg(distinct es.alt_strand::varchar, ',') as strand,
transcript.cds_start_i,
transcript.cds_end_i,
string_agg(exon.start_i::varchar, ',' order by exon.ord) as exon_starts,
string_agg(exon.end_i::varchar, ',' order by exon.ord) as exon_ends,
string_agg(exon_aln.cigar, ',' order by exon.ord) as cigars
from uta_20210129.transcript transcript
inner join uta_20210129.exon_set es on (transcript.ac = es.tx_ac AND alt_aln_method = 'splign')
inner join uta_20210129.origin origin on (transcript.origin_id = origin.origin_id)
Inner join uta_20210129.exon as exon on (es.exon_set_id = exon.exon_set_id)
inner join uta_20210129.exon_aln exon_aln on (exon_aln.alt_exon_id = exon.exon_id)
WHERE es.alt_ac in
('NC_000001.10', 'NC_000002.11', 'NC_000003.11', 'NC_000004.11', 'NC_000005.9', 'NC_000006.11', 'NC_000007.13', 'NC_000008.10', 'NC_000009.11', 'NC_000010.10', 'NC_000011.9', 'NC_000012.11', 'NC_000013.10', 'NC_000014.8', 'NC_000015.9', 'NC_000016.9', 'NC_000017.10', 'NC_000018.9', 'NC_000019.9', 'NC_000020.10', 'NC_000021.8', 'NC_000022.10', 'NC_000023.10', 'NC_000024.9') and origin.origin_id not in (10, 11)
group by transcript.ac) TO 'uta_20210129_grch37.csv' CSV HEADER;