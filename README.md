# masstodon

Your Python3 module for investigating the Electron Transfer Dissociation in Mass Spectrometry, and more generally finding stuff in the mass spectrum.

# Prerequisites
Python3

# Installation

The package can be installed directly from the Python Package Index:
```
pip install masstodon
```

Otherwise, download this github repo, use the default branch, and install the software with:

```
pip install .
```
from the folder containing **setup.py**.
This will let you import `masstodon` simply by

```{python}
import masstodon
``` 

# Running

## command line interface

It is possible to run `masstodon` from command line after installing it from PIP.
Given that it would be not comfortable to input multiple compounds with their possible modifications and charges from the terminal, we have restricted (for now) the possibility to run `masstodon` in the traditional case of studying precisely one protein and its c/z fragmentation patterns.

To see a detailed description of the possible arguments, open the terminal (or the Anaconda prompt on Windows, but seriously, Windows? You should feel guilty...) and write
```{bash}
    masstodon -h
```

To run `masstodon`, you do need to supply at least:
* the file with the spectrum
* the tolerance for the search in the m/z axis
* the amino acid sequence (fasta)
* the charge of the substance

All other parameters can be skipped, but it might be silly.
For instance, if you know that there is a particular PTM on a given amino acid,
you should supply the `-m` parameter, and so on.


## Python scripting

**Importing data.** You can analyze individual mass spectra with `masstodon`.
This might very well be individual scans of an *Orbitrap*, or a general mass spectrum.
The easiest way to import the mass spectrum is to present a plain **ASCII** file with m/z and intensity values in each row, separated by tab or some other whitespace sign.

**Running masstodon**.
To call `masstodon` with a single compound, use the `masstodon_single` function.

For instance:

```{python}
from masstodon.data.substanceP_wh15_wv400 import mz, intensity, fasta, modifications, z # input data
from masstodon.masstodon import masstodon_single # masstodon-proper
from pprint import pprint # nicer output in the terminal

print(mz)
# array([  61.01 ,   64.152,   66.515, ..., 1403.9  , 1485.578, 1497.226])
print(intensity)
# array([844.4 ,  25.35, 190.1 , ...,  15.38,  55.62,  21.  ])
print(fasta)
# 'RPKPQQFFGLM'
print(modifications)
# {'11': {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}}
# (modify 11th amino-acid [starting from 1 for user-friendliness],
# by adding the chemical diff unto the group of atoms that include the C_carbo.)
print(z)
# 3

m = masstodon_single(mz, intensity, fasta, z,
                     min_mz_diff        = 1.1, # this is actually the default
                     modifications      = modifications,
                     orbitrap           = False, # this ain't an orbitrap profile spectrum
                     threshold          = "0.05 Da", # how far away to search for
                     isotopic_coverage  = .999, # IsoSpec isotopic envelope coverage
                     min_prob           = .8, # minimal acceptance probability of a candidate theoretical envelope.
                     std_cnt            = 3 # m/z standard deviation count filter)
# for min_prob and std_cnt check out the publication under the filtering procedures.

print("Save spectrum to where you started your python from.")
m.plotly("spectrum.html", shape='rectangles', show=True)

print("Getting stats on inintial and final number of nodes and edges in the deconvolution graph.")
pprint(m.ome.G_stats)

print("Errors: absolute deviations and relative distance between spectra.")
pprint(m.imperator.errors())

print("Estimates of intensities of reactions and fragmentations.")
pprint(m.cz_simple.intensities)
pprint(m.cz_simple.probabilities)

print("Save all things under the given path.")
m.write(".")
```

<!-- # Running *masstodon*

There are now two ways to run the program:

1. In terminal.
2. As part of another Python script.

### Terminal Call

To run MassTodonPy in terminal, simply type

```
masstodon <spectrum> <fasta> <charge> <mz_tol> -o <output_path>
```
where: 
* <spectrum> is the path to the file containing the spectrum: *mzXml*, *mzml*, or a raw txt file with two columns (recorded m/z and intensities):

```bash
191.932 17.36
271.183 98.33
415.8 17.23
425.948 15.21
444.232 208.4
444.359 6.41
444.568 117.6
445.236 44.26
449.284 19.72
... ...
```

* <fasta> is the fasta amino acidic sequence, e.g. AAA of the molecule you look for in the spectrum.
* <charge> is the charge of the observed precursor
* <mz_tol> is the distance from the theoretical m/z that MassTodon will look around for signals.
* <output_path> is the path to where you want to write the output of the software.


For a full list of possibilities type:
```bash
masstodon -h
```

Say you want to add a modification to your fasta specified protein.
Suppose that the original peptide was AAAGGGVVAGV, had 2 charges, and included a C-terminal amidation,
i.e. a replacement of -COOH with -CONH2. 
The modification diff consists of a change we can symbolically write as O=-1, N=1, H=-1.
The C-terminal is on the right valine, which is the eleventh amino acid counting from the N-terminal.
We label the atoms included in the backbone of any amino acid as N, C\_alpha, and C\_carbo. 
The amidation only modifies the C_carbo atom, i.e. the one present in the carboxyl group -COOH.
Calling masstodon would look like:

```bash
masstodon spectrum.mzXml AAAGGGVVAGV 2 .01 -modifications '11 C_carbo H=-1 N=1 O=-1'.
```

It is crucial to specify the positionment of modification w.r.t. correct backbone atom,
as the fragmentations severes bonds between them.
If you specify it wrongly, then the diff will modify part of the resulting fragments.

MassTodonPy supports the insertion of multiple modifications:
```bash
masstodon spectrum.mzXml AAAGGGVVAGV 2 .01 -modifications '11 C_carbo H=-1 N=1 O=-1 | 2 N H=-1 Li=1'.
```

The above silly modification replaces a hydrogen atom with a lithium atom (pardon my chem-lish).


### Visualizing Spectra
MassTodonPy includes a seperate submodule to plot raw spectra using Bokeh library.
By default, it is added to the *bin* folder, and can be called:

```bash
plot_mass_spectrum <spectrum>
```

For more options (including properties of the browser plot), type 
```bash
plot_mass_spectrum -h
```

### Python Scripting

The simplest way to use the MassTodon in your Python script is to import the **MassTodon** function from the MassTodonPy module. 
A simple script used to run the previous example peptide, would look like this:

```python
from MassTodonPy import MassTodon

res = MassTodonize( fasta           = AAAGGGVVAGV,
                    precursor_charge= 2,
                    mz_prec         = .01,
                    spectrum        = "spectrum.mzXml",
                    modifications   = {11: {'C_carbo': {'H': 1, 'O': -1, 'N': 1}}} )

res.write('output_path')
```

This will run the software and save the results to 'output_path'.

To visualize the outputs, add:

```python
from MassTodonPy.Plot import bokeh_spectrum
from MassTodonPy.Plot import bokeh_aggregated_precursors
from MassTodonPy.Plot import bokeh_aggregated_fragments
from MassTodonPy.Plot import bokeh_estimated_aggregated_fragments

bokeh_spectrum(res)
bokeh_aggregated_precursors(res)
bokeh_aggregated_fragments(res)
bokeh_estimated_aggregated_fragments(res)
```
 -->