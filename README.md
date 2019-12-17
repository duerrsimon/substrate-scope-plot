# Substrate-scope-plot

Inspired by https://twitter.com/NarayanLab/status/1206652692321292288
I made a matplotlib version using RDKit. 

Simply get a SMILES string from ChemDraw or save a svg in ChemDraw, name it rdkit.svg and put it in the root directory of the of this document. 

Use this within a Python Anaconda distribution, e.g Miniconda https://docs.conda.io/en/latest/miniconda.html


To install RDKIT

    conda install -c conda-forge rdkit

If you do not have matplotlib etc. do

    conda install matplotlib svgutils

Clone this repo using git

    git clone https://github.com/duerrsimon/substrate-scope-plot.git

Then fire up a jupyter-notebook from the terminal like so

    cd substrate-scope-matplotlib
    jupyter-notebook substrate_scope.ipynb
    
Execute cells (Ctrl+Enter) after entering all variables
