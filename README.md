# pymolint
A python package to characterise molecular interactions at atomic level in PDB structures and MD simulations
![](https://github.com/intbio/pymolint/workflows/Testing/badge.svg)
## Installation via conda
```
conda install -c intbio -c conda-forge -c bioconda pymolint
```


## Usage
### In Jupyter
For Jupyter use see [example.ipynb](example.ipynb) - THIS IS RECOMMENDED




## Debuging and common problems


## For developers 
---- In development
- Testing can be done as
```
git clone https://github.com/intbio/pymolint.git
docker run   --workdir "/wd" -v "$PWD/pymolint:/wd" intbio/pymolint_test pytest
```
See test_results folder for results.
