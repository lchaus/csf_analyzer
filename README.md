# CSF Analyzer

This is a command-line tool to analyze wavefunctions from ORCA (CASSCF, ICE-CI...).

## Capabilities

* Parse the weights and occupations of state-averaged calculations in ORCA
* Classify different types of excitations after reading user-defined orbital blocks to evaluate charge transfer
* Classify degrees of excitation with respect to a given reference (coming soon)

## Usage
First, you need to run the analysis. You will need a text file "orbitals.dict", storing orbitals names and there corresponding label in ORCA (Warning: 0 based), and the printed wavefunction from ORCA, "orca.out".
Here is an example of how 


We suggest to print the wavefunction using "TPrint 1e-6" in the CI block and trim it with a text editor or suitable script.

Then run:
```python $path_csf_analyzer/main.py analyze orca.out orbitals.dict```
