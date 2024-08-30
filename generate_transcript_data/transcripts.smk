include: "cdot_utils.smk"

annotation_consortium = config["annotation_consortium"]
genome_build = config["genome_build"]
urls = config["urls"]
cdot_output_dir = os.path.join(annotation_consortium, genome_build)

cdot_file_template = "cdot-" + cdot_data_version + "-{name}.json.gz"

def get_cdot_command(wildcards):
    url = urls[wildcards.name]
    cdot_command = "gff_to_json" if url.endswith(".gff.gz") else "gtf_to_json"
    return cdot_command


rule all_transcripts:
    input:
        expand(os.path.join(cdot_output_dir, cdot_file_template), name=urls.keys())

rule download_files:
    output:
        # Don't re-download if snakemake script changes
        protected("downloads/{name}.gz")
    params:
        url=lambda wildcards: urls[wildcards.name]
    shell:
        "curl -o {output} {params.url}"

rule cdot_json:
    input:
        gene_info_json=gene_info_json,
        gff_file="downloads/{name}.gz"
    output:
        protected(os.path.join(cdot_output_dir, cdot_file_template))
    params:
        url=lambda wildcards: urls[wildcards.name],
        cdot_command=get_cdot_command
    shell:
        """
            PYTHONPATH={cdot_dir} \
            {cdot_json} \
                {params.cdot_command} \
                "{input.gff_file}" \
                --url "{params.url}" \
                --genome-build="{genome_build}" \
                --output "{output}" \
                --gene-info-json="{input.gene_info_json}"
        """
