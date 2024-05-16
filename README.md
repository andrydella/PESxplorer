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
* CREST, short for Conformer-Rotamer Ensemble Sampling Tool, is a utility and driver program for the semiempirical quantum chemistry package. Its functionalities include but are not limited to: conformational sampling, ensemble sorting, costrained sampling, protonation site calculation and conforomational entropy calculation. For more details, the documentation can be found [here](https://crest-lab.github.io/crest-docs/page/examples). The [installation guide](https://crest-lab.github.io/crest-docs/page/installation/install_compile.html) can be followed to compile from source code.
* GSM, the Growng String Method is a reaction path (RP) and transition state (TS) finding method that develops a RP by iteratively adding new nodes (each node is a chemical structure along the RP) and optimizing them until a complete RP with a TS and a stable intermediate on each side of the string are present. More details and instructions on installing and running the code can be found [here](https://github.com/ZimmermanGroup/molecularGSM/wiki)
* [Gaussian](https://gaussian.com/). It is used by thr GSM to compute energy gradients. It can be replaced with another electronic structure calculations software, as the GSM allows to interface with MOLPRO, QCHEM, MOPAC, ASE or GAUSSIAN.
* Python v.3.7 or higher. That's the python version I used, but I limited the use of python libraries so there shouldn'r be any compatibility problems with Python version > 3.0

## WORKFLOW
The main working directory will be referred to as CWD.
### User input
* In CWD, generate a subfolder in which all CREST calculations and output files will be collected: CWD/MD_PATH
* Go into MD_PATH and create the necessary input files for CREST:
  * `starting_geometry.xyz` : it contains the cartesian coordinates of the initial structure (it could be found online or generated with RDKIT starting from InChI, SMILES or MOLBLOCK file of a molecule.
  * `.UHF` and `.CHRG` : contain spin and charge of the molecule. necessary only if not default (0 for both spin and charge)
* Consider running an initial optimization of your structure
```bash
crest starting_geometry.xyz --opt --gfn2 -T 36 > crst_opt.out &
```
* Run the `MSREACT` protocol with either of the following commands taken from the official documentation of CREST (the second one could be used for more complex systems, but trying them both is possible, as the outputs can be manually combined in a unique file. --T is the number of cores used by the software.
```bash
crest starting_geometry_opt.xyz --msreact --msmolbar -T 36 > crst_samp.out &
crest starting_geometry_opt.xyz --msreact --msnshifts 5000 --msnshifts2 50 --msmolbar --msiso -T 36 --ewin 500.0 > crst_samp.out &&
```
* The keyword `--msiso` tells CREST to collect the unique isomers that it identifies in the `isomers.xyz` trajectory file.
* If using an older version of CREST (< 3.0), you can set up different metadynamics simulation in subfolders CWD/MD_PATH/metadyn1, CWD/MD_PATH/metadyn2 etc. The prefix of the metadynamics subfolders must be always the same (e.g. "metadyn" in this case. A few parameters need to be tweaked in `explorer.py`:
  * in the function `unite_xyz`, the name of the CREST output file is hard coded and is currently set to crest_products.xyz. Make sure it is correct.
  * in the `main` function, edit the variables `crest_version`, `md_path` (the relative path to MD_PATH is sufficient) and `crest_out_name`. This last argument is set to "isomers.xyz" for the newer CREST version (> 3.0), but needs to be set to the prefix of the metadynamics subfolders, in this case "metadyn".
* In the `main` function, the calls to the functions can be commented to run the script step by step, or left as is to perform all the steps at the same time.
* `explorer.py` is called from CWD, and here all the intermediate and output files and subfolders will be generated. It is suggested to comment the last function that calls for the GSM runs initially, to make sure that the previous results are satisfactory.
* Currently, `explorer.py` runs the GSM sfotware through a bash script placed in the `bin` that sends the job to a queue management system. "rungsm" can be replaced with any other command or script that executes the GSM. The counter that follows refers to the index of the GSM input files, without the preceeding zeros (GSM input files are in the form "name_of_file_0001.xyz".
* the variable model_data_gsm must be the path to the model data folder that contains the necessary files for setting up the GSM calculation.
