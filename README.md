# NBIS Support \#5568-cray

This is the home of the NBIS support project with the Redmine issue \#5568,
*"Single cell-RNA sequencing for identifying hemocyte subpopulations"*.

A Conda environment specification is provided in the environment.yml file and
a Nextflow workflow allows reproducibility of the analyses. Data is not stored
in this repository.

## Reproducibility

The transcriptome FASTA and BLAST XML files used in these analyses come directly
from the Söderhall group's previous work, while the mitochondrion genome can be
found at the [NCBI website](https://www.ncbi.nlm.nih.gov/nuccore/NC_033509.1/).
The analyses can be reproduced using either Docker or Conda.

Nextflow is run from it's own conda environment whether using docker, singularity, or conda packages. The conda environment is provided
as a shared environment for all members of the project.
```bash
# conda create --prefix "/proj/snic20XX-YY-ZZ/NBIS_<id>/conda_nextflow-env" -f environment.yml
conda activate /proj/snic20XX-YY-ZZ/NBIS_<id>/conda_nextflow-env
nextflow run ...
conda deactivate
```

### Using Docker

Public images are used where possible

```bash
# Run the pre-processing workflow
nextflow run . -profile uppmax,docker -entry pre_processing

# Run the analysis workflow
nextflow run . -profile uppmax,docker
```

### Using Conda

```bash
# Run the pre-processing workflow
nextflow run . -profile uppmax,conda -entry pre_processing

# Run the analysis workflow
nextflow run . -profile uppmax,conda
```
