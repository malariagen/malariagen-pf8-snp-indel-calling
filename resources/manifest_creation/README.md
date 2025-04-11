# Manifest Creation

## Install dependencies in a virtual environment

Using virtualenv:
```bash
virualenv -p python3 venv
source venv/bin/activate
pip install numpy pandas petlx mysqlclient
```

Using conda:
```bash
conda create -n create_manifest python numpy pandas petl mysqlclient
conda activate create_manifest
```

## Preparation

### Create a directory containing manifests used for previous builds
This script will only include samples that have not been included in previous manifests/data releases, or samples for which new lanes have been added since a previous release. To ensure that the correct samples are included, you need to point it to a directory containing manifest files for previous releases. To do this, invoke the script with the `-p <dir>` or `--previous_imports_dir <dir>` (default: current working directory).

Manifests for previous builds can be found at
```
<INSERT PATH HERE>/parasite-ops/work/80_Pf_6_3_initial_irods_manifest/pf_60_import.txt
<INSERT PATH HERE>/parasite-ops/work/80_Pf_6_3_initial_irods_manifest/pf_61_import.txt
<INSERT PATH HERE>/parasite-ops/work/80_Pf_6_3_initial_irods_manifest/pf_62_import.txt
<INSERT PATH HERE>/parasite-ops/work/80_Pf_6_3_initial_irods_manifest/pf_63_irods_manifest_20191015.txt
``` 

Manifest files stored in the given `<dir>` should be named pf_6x_import.txt, where x is substituted for a build number (0,1,2,3).

If you need to include another manifest (because, for example, new lanes have been added since you last created the manifest), you will also need to amend the script by adding another item to the dictionary `previous_imports`.

### Supplying database login credentials
Without any other modifications, when running the script, you will prompted for username/password login credentials for databases that are required. The script uses two mySQL databases, the `subtrack` database on host `shap` (ENA submission database) and the `mlwarehouse` on host `mlwh-db-ro` (iRODS database). You can find these credentials in [pathdev-passwords](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pathdev_passwords) repository before running, but any credentials with appropriate access will work.

To avoid typing in credentials interactively, the script will also read them from environment variables `SUBTRACK_USER`, `SUBTRACK_PASSWORD`, `MLWH_USER`, `MLWH_PASSWORD`. If you set these, remember to set them in such a way that they are not recorded in the logs/process history (e.g. set them in a wrapper script, or a file that you source before running the `create_manifest.py` script, and restrict permissions for these files with `chmod 700 <file_with_credentials>` or similar).

## Running the script
To run the script, log into the farm, activate the relevant virtual environment, then run:
```bash
./create_manifest.py -p <manifests_dir> -o <path_to_output_manifest>
```

Where \<manifests_dir\> is substituted as described in [Preparation](#preparation). By default, the output file `manifest.txt` will be created in the current working directory.

It's unclear what role `write_gvcf_entry_manifest.py` played in Pf8, but will leave as is. 