# CSF Analyzer

This is a command-line tool to analyze wavefunctions from ORCA (CASSCF, ICE-CI...).

## Capabilities

* Parse the weights and occupations of state-averaged calculations in ORCA.
* Classify different types of excitations after reading user-defined orbital blocks to evaluate charge transfer.
* Classify degrees of excitation with respect to a given reference.

## Usage

### Charge-transfer evaluation
This program is designed to evaluate charge-transfer contributions to a localized multiconfigurational wavefunction. Preferably, one should include in the active space and localize full orbital shells (e.g. 3d, 4d...). The program will evaluate, with respect to the main configuration of each root, local contributions (electron count is the same in every block), and charge-transfer contributions (sum of the block to block electron transfer). We can later plot, in the form of pie charts or bar plots:

* The main CSF (Tn).
* Local contributions (all summed together).
* All block to block contributions above a given threshold (e.g. 3d -> 4d, 3d -> 3p).

Orbital blocks and their labels in the ORCA wavefunction are defined by the user in a file ```orbitals.dict``` that you should keep in this format where you will run the program:

```
{
    r"Cu}^{3s3p}": [0,1,2,3], 
    r"Cu}^{3d}": [4,5,6,7,8],
    r"Cu}^{4s4d}": [13,21,23,24,25,33],
    r"O}_2^{2p}": [9,10,11,12,14,15],
    r"O}_2^{3p}": [26,27,28,29,30,31],
    r"N}_{bonding}": [16,17,18,19],
    r"N}_{antibonding}": [20,22,32,34],
}
```
You also have to print the wavefunction with ORCA (for example using "TPrint 1e-6" in the CI block) and copy the output file where you will run the program.

Then run:

```python $path_csf_analyzer/main.py analyze_ct cuo2.out orbitals.dict```

This creates a file ```excitation_classes.txt``` containing all charge-transfer information for each root in the wavefunction.

To visualize the results you can run:

```python $path_csf_analyzer/main.py visualize excitation_classes.txt```

Which will create a directory ```./plots``` containing bar plots and pie charts of the data.
You can control the threshold for the plots with (TODO)



### Excitation degree evaluation


