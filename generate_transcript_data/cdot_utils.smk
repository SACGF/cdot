import subprocess
from datetime import datetime

cdot_json = os.path.join(workflow.basedir, "cdot_json.py")
cdot_dir = os.path.dirname(workflow.basedir)
cdot_output_raw = subprocess.check_output(f"{cdot_json} --version", shell=True, env={"PYTHONPATH": cdot_dir})
cdot_data_version = cdot_output_raw.decode().strip()

# Name it based on date as it may vary
today = datetime.now().date().isoformat()
gene_info_download_filename = f"Homo_sapiens.gene_info.{today}.gz"
gene_info_json = f"Homo_sapiens.gene-info-{cdot_data_version}.json.gz"
