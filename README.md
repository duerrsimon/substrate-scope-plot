# Radial Scope plot 
### The matplotlib/rdkit version

#### Orginal idea
This figure made its rounds on twitter based on this tweet https://twitter.com/NarayanLab/status/1206652692321292288

Original credit goes to Attabey Rodríguez Benítez [@ScienceBey](https://twitter.com/ScienceBay) and Alisson Narayan 

This is a completely automated matplotlib version of this plot using RDKit for generating the molecules based on SMILES strings.

### How does it work

1.Get a SMILES string from ChemDraw  or any other software able to output SMILES strings (such as [PubChem](https://pubchem.ncbi.nlm.nih.gov/edit3/index.html).
For every $R_x$ group put a methyl group.
2. Setup the python dictionary for the settings
3. Setup a dictionary for each $R_x$ you want to replace with a radial scope plot. Three examples can be found below
4. Execute all cells the Run button or Ctrl+Enter

### Installation

Use this within a Python Anaconda distribution, e.g [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
Anaconda is available on Windows, Mac and Linux and this should work on all platforms.

To install RDKIT

    conda install -c conda-forge rdkit

If you do not have matplotlib etc. do

    conda install matplotlib svgutils

Clone this repo using git

    git clone https://github.com/duerrsimon/substrate-scope-plot.git

Then fire up a jupyter-notebook from the terminal like so

    cd substrate-scope-plot
    jupyter-notebook radial_scope_plot.ipynb
    
### How does it look like
![radial scope plot example](https://raw.githubusercontent.com/duerrsimon/substrate-scope-plot/master/substrate_scope_replaced.svg "radial scope plot example")


#### Changelog 
v0.2| 23.12.2019
- put all the plotting in a class that is called from the notebook
- multiple radial scope plots possible 
- correct automatic positioning based on the index of the atom from the SMILES
- addition of rendering of small organic rests (hacky, the user will need to reposition the images using a vector software)
- added some options such as
    - rounding options
    - bw atom theme
    - bold font
    - rounding based on threshold value
- text becomes white at 80% of max value

v0.1 | 17.12.2019
- first rough reproduction of the figure using matplotlib
