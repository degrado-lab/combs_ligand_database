# combs_ligand_database
Parse a mirror of the PDB for protein-ligand contacts.

### using the code
This code is intended to generate the necessary pickled objects for the generation of van der Mer databases 
using the Combs2 software package for the design of ligand-binding proteins.  A mirror of the PDB and 
validation information for each structure is a necessary prerequisite.  These can be downloaded from the 
PDB FTP server as follows:

```bash
> rsync -rlpt -v -z --delete --port=33444
  rsync.rcsb.org::ftp_data/structures/divided/pdb/ $LOCAL_PDB_MIRROR_PATH
  
> rsync -rlpt -v -z --delete --include="*/" --include="*.xml.gz" --exclude="*"
  --port=33444 rsync.rcsb.org::ftp/validation_reports/ $LOCAL_VALIDATION_PATH
```

Binaries for the probe, rotalyze, and reduce software from the Richardson Lab are also necessary.  These 
can be downloaded as part of the MolProbity software package (https://github.com/rlabduke/MolProbity) and 
paths to the binaries can be passed as command-line arguments to ligand_database.py

Lastly, the Combs2 conda environment has all the necessary prerequisite Python packages for 
ligand_database.py, the script that provides the primary user interface of this package, to run.

ligand_database.py was designed for compatibility with the directory structure of PDB mirrors from the 
PDB FTP server, which groups structures based on the middle two characters of their accession code. As 
such, the code should be run independently on each of the 1060 directories in the mirror. Output from 
ligand_database.py, parallelized in batch jobs on the UCSF Wynton cluster, is available at 
https://ucsf.box.com/s/5xtg7j4ekhrtqba61h3hzbprwp8u9qd7
