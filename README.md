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

### Using Docker

```bash
# Build the Docker image
docker build -t nbis-5568 .

# Run the workflow
nextflow run . -use-docker nbis-5568
```

### Using Conda

```bash
# Create and activate the Conda environment
conda env create -p 5568-env -f environment.yml
conda activate 5568-env/

# Run the workflow
nextflow run .
```
