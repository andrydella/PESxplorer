# PESxplorer
Pesxplorer is a tool that creates a synergy between existing software from the literature to achieve the exploration of Potential Energy Surfaces (PES) by identifying all unique conformers, possible pathways among them that involve the breakage and formation of up to one bond and finally attempting at finding saddle points that connect couples of isomers. 
The workflow is structured as follows:
1. Generate one or more lists of conformers by setting up multiple nanoreactor simualtions with CREST, or the recently added (CREST version 3.0) MS-REACT keyword.
2. Run <explorer.py> to
  * collect the cartesian coordinates of all conformers from the output/s of CREST in a unique text file
  * Sort it by energy
  * Select unique isomers by computing a connectivity matrix for each conformer
  * Exclude bimolecular products if present in the original output file
  * Compute all the possible permutations of lines containing the same atom for each isomer and perform a further check through the connectivity matrix of each possible permutated isomer
  * Create a list of all possible reaction pathways and eliminate the doubles
  * Set up subfolders for each species and preapare the necessary input for each unique isomerization channel identified
  * Atuomatically run the Growing String Method (GSM) softwre to reliably find possible transition state (TS) structures connecting the minima. The GSM code relies on Gaussian for the computation of gradients.
3. Post process the results by checking the reaction paths identified by the GSM. This step is manual, as many different cases may arise, such as: successful run and TS identified, multiple TSs found connecting reactant and product, the reactant converged to a different isomer at the new level of theory used for Gaussian calculations of the gradient of the energy.
## Requirements
PESxplorer exploits existing software from the literature to achieve the identification of minima on a potential energy surface, characterized by the same brute formula, e.g. C<sub>x</sub>O<sub>y</sub>N<sub>z</sub>H<sub>i</sub>, but the code is not limited to these atoms only. Any atom type that can be properly described by the theoretical mehtodologies employed by CREST and Gaussian can be included in the molecule. 
* CREST
* GSM
* Gaussian
* python v.3.7 or higher
