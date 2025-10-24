# After 0.2.22 we split version into separate code (pip) and data schema versions
# The cdot client will use its own major/minor to determine whether it can read these data files

# DATA changelog (TODO: Move this to changelog)

# 0.2.29 - Ensembl now has HGNC added from outside GTFs
# 0.2.30 - Ensembl GRCh37 has canonical transcripts added from outside GTFs
# 0.2.31 - Add 'metadata' - method/urls
# 0.2.32 - Add 'source' (GTF column #2) to build data
# 0.2.33 - Store "ccds", "transcript_support_level"
JSON_SCHEMA_VERSION = "0.2.33"
