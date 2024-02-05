# FabA33DataProcessing

The data processing for the simulation data of Fab A33



## Feature calculation of Fab A33

## Content

SASA

Non-polar SASA

Sum of APR ΔSASA of the 7 APRs 

Native contacts

RMSF

Net charges

Salt bridges

Hydrogen bonds



### The rationale behind the molecular features selected for model building


| Feature                     | Rationale                                                                                                                                                            | Action                                                                                                                                          |
|-----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| SASA (global)               | SASA is a measure of how much of the area of a molecule is available to the solvent. The atomic SASA can be used as a measure of the steric availability of an atom. | Calculate the polar SASA for all the atoms inside one residue and then plus them together to get a final polar SASA at a residue level. GROMACS |
| Non-polar SASA              | SASA-NP is the accessible surface area of all non-polar atoms.                                                                                                       | GROMACS + index file                                                                                                                            |
| APR                         | APR score tells you which residue is most solvent accessible                                                                                                         | Comput Struct Biotechnol J. 2021; 19: 2726–2741. Cheng Zhang.                                                                                   |
| Fraction of native contacts | The percentage of the contacts that exist within a native state in every frame                                                                                       | MDanalysis package                                                                                                                              |
| RMSF                        | It measures the fluctuation of an atom or e.g., a protein residue along the course of a simulation, and is useful to identify the most mobile regions                | GROMACS                                                                                                                                         |
| Net charges                 | It indicates the net charge of the protein at different pHs                                                                                                          | PropKa                                                                                                                                          |
| Salt bridges                | It is the average of the salt bridge occurrence in the last 50ns                                                                                                     | MDanalysis                                                                                                                                      |
| Hydrogen bonds              | It measures the number of hydrogen bonds over time                                                                                                                   | GROMACS                                                                                                                                         |


## Methods

#### Regression analysis

* [PCA](http://thegrantlab.org/bio3d/reference/pca.xyz.html)

* [DCCM](http://thegrantlab.org/bio3d/reference/dccm.html)

* Pearson correlation coefficient


#### Model building

* [Lazy predict](https://github.com/shankarpandala/lazypredict)

* [SHAP](https://shap.readthedocs.io/en/latest/)



## Programming languages
### Python
python packages
* pandas
* Dtale
* re
* fnmatch
* lazy predict
* seaborn
* matplotlib.pyplot
  

### R
R packages
* [bio3d](http://thegrantlab.org/bio3d/)

### TCL

### Bash

