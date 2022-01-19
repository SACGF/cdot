
# if PyHGVS below a certain version - do it the old way

# If after - can use different exons, and start/stop codons to make it faster


# Changes from old loading:

# See dot has no cds_start/end if non-coding
#  PyHGVS expects cds_start/cds_end be equal to end/end for non-coding transcripts (so coding length ie end-start = 0)
#     cds_start = transcript_data.get("cds_start", end)
#     cds_end = transcript_data.get("cds_end", end)



# VG loader also expects biotype to be comma sep, now is list