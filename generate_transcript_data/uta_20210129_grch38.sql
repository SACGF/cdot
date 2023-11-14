\copy (SELECT transcript.ac, string_agg(distinct transcript.hgnc, ',') as hgnc, 'http://www.ncbi.nlm.nih.gov/refseq/' as origin_url,
string_agg(distinct aln_v.alt_ac::varchar, ',') as contig,
string_agg(distinct aln_v.alt_strand::varchar, ',') as strand,
transcript.cds_start_i,
transcript.cds_end_i,
string_agg(aln_v.alt_start_i::varchar, ',' order by aln_v.alt_exon_id) as exon_starts,
string_agg(aln_v.alt_end_i::varchar, ',' order by aln_v.alt_exon_id) as exon_ends,
string_agg(aln_v.cigar, ',' order by aln_v.alt_exon_id) as cigars,
string_agg(distinct aa.pro_ac, ',' order by aa.pro_ac) as protein
from uta_20210129.transcript transcript
inner join uta_20210129.tx_exon_aln_v aln_v on (transcript.ac = aln_v.tx_ac AND alt_aln_method = 'splign')
left outer join uta_20210129.associated_accessions aa on (transcript.ac = aa.tx_ac)
WHERE aln_v.alt_ac in
('NC_000001.11', 'NC_000002.12', 'NC_000003.12', 'NC_000004.12', 'NC_000005.10', 'NC_000006.12', 'NC_000007.14', 'NC_000008.11', 'NC_000009.12', 'NC_000010.11', 'NC_000011.10', 'NC_000012.12', 'NC_000013.11', 'NC_000014.9', 'NC_000015.10', 'NC_000016.10', 'NC_000017.11', 'NC_000018.10', 'NC_000019.10', 'NC_000020.11', 'NC_000021.9', 'NC_000022.11', 'NC_000023.11', 'NC_000024.10') and origin.origin_id not in (10, 11)
group by transcript.ac) TO 'uta_20210129_grch38.csv' CSV HEADER;