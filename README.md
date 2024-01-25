# FabA33DataProcessing

The data processing for the simulation data of Fab A33

[![Open in VS Code](https://img.shields.io/badge/Visual_Studio_Code-0078D4?style=flat&logo=visual%20studio%20code&logoColor=white)](https://open.vscode.dev/man-group/dtale)

## Feature calculation of Fab A33

## Content
SASA

Native contacts

RMSF, RMSD, Gyrate

Distribution

Net charges

Salt bridges

Hydrogen bonds

### The rationale behind the molecular features selected for model building


<img width="853" alt="image" src="https://github.com/wintermute221/FabA33DataProcessing/assets/57851709/f3529767-740a-4577-8763-751fe879606c">

| Feature                     | Rationale                                                                                                                                                            | Action                                                                                                                                            |
|-----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| SASA (global)               | SASA is a measure of how much of the area of a molecule is available to the solvent. The atomic SASA can be used as a measure of the steric availability of an atom. | Calculate the polar SASA for all the atoms inside one residue and then plus them together to get a final polar SASA at a residue level. GROMACS   |
| SASA-NP                     | SASA-NP is the accessible surface area of all non-polar atoms.                                                                                                       | GROMACS + index file                                                                                                                              |
| SASA-P                      | SASA-P is the accessible surface area of all polar atoms.                                                                                                            |                                                                                                                                                   |
| APR                         | APR score tells you which residue is most solvent accessible                                                                                                         | APR tools (extracted from Nesrine)                                                                                                                |
| Fraction of native contacts | The percentage of the contacts that exist within a native state in every frame                                                                                       | MDanalysis package                                                                                                                                |
| RMSF                        | It measures the fluctuation of an atom or e.g., a protein residue along the course of a simulation, and is useful to identify the most mobile regions                | GROMACS                                                                                                                                           |
| Hydropathy                  | The hydropathy index of an amino acid is a number representing the hydrophobic or hydrophilic properties of its sidechain.                                           | - Maybe discuss this Paul. I always think there should be a better, more insightful, deeper explanation of the reason why people use this feature |

## Methods

Regression analysis

* [PCA](http://thegrantlab.org/bio3d/reference/pca.xyz.html)

* [DCCM](http://thegrantlab.org/bio3d/reference/dccm.html)

* 

## Programming languages
### Python
* python packages
* [pandas]()
* [Dtale]()

### R
* R packages
* [bio3d](http://thegrantlab.org/bio3d/)

### TCL

### Bash

