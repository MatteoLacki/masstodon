# MassTodonPy

![Alt](/docs2/source/figs/masstodon_cracowian.jpg)

A Python module investigating the Electron Transfer Dissociation in Mass Spectrometry

# Prerequisites

In order to play with MassTodonPy on your computer shall need:
1. UNIX based operating system, such as Linux (eg. Ubuntu, Fedora, Gentoo, or similar) or macOS (eg. macOS Sierra)
2. python2.7 interpreter

To check if you have the correct interpreter, open the terminal and simply type

```{bash}
python2.7
```
and you should get:

![Alt](/Webpage/python_terminal.png)

To leave the terminal press cntr+d

To install globally MassTodonPy, run in terminal
```{bash}
pip2 install MassTodonPy
```

This will install MassTodonPy and the required dependencies.
To check the installation, simply write now:

```{bash}
masstodon_example_call
```

This will run an example session of the program to check that you are able to get any output. By no means is it necessary to run it more than once, just after installation.

Great, you seem to have all the tools in place for some massive calculations!

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
