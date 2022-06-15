from snakemake.utils import validate
import pandas as pd

##### load config and sample sheets #####
configfile: "config/config.yaml"

samplesheet = pd.read_table(config["samples"], sep="\t").set_index("submitted_donor_id", drop=False)
samplesheet.index.names = ["sample_id"]

#Predefine output folders
OUT_DIR=config["out_dir"]

BIN_VALS = config["bins"]
