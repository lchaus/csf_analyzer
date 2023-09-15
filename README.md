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

Orbital blocks and their labels in the ORCA wavefunction are defined by the user in a file ```orbitals.dict``` that you should keep in this format where you will run the program (to preserve matplotlib LaTeX formatting):

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
You also have to print the wavefunction with ORCA (for example using "TPrint 1e-6" or lower in the CI block) and copy the output file where you will run the program.

Then run:

```python $path_csf_analyzer/main.py analyze_ct orca.out orbitals.dict```

This creates a file called ```excitation_classes.txt``` containing all charge-transfer information for each root in the wavefunction.

To visualize the results you can run:

```python $path_csf_analyzer/main.py visualize excitation_classes.txt```

Which will create a directory ```./plots``` containing bar plots and pie charts of the data.
You can also control the threshold for plotting a given contribution in either bar plots of pie charts using respectively options ```-tb 0.01``` or ```-tp 0.01```:

```python $path_csf_analyzer/main.py visualize -tp 0.01 excitation_classes.txt```

### Excitation degree evaluation

Alternately, you can run an analysis of the wavefunction in terms of degrees of excitation with respect to a given reference (e.g. ROHF, LMCT...).
The program will classify all configurations as such:

* ROHF (user defined)
* LMCT (user defined)
* Reference configurations (In case of multireference states, all configurations above a certain user-defined threshold, with important contributions, are printed separately (e.g. d-d excitations)
* Single, Double, Triple, Quadruple...  (with respect to ROHF confifuration)
* LMCT + Single are a special class of double excitations that are printed separetly (TODO: add option to deactivate)

To run the program, you will need, in this order:

* ORCA output file with the printed wavefunction (for example using "TPrint 1e-6" or lower in the CI block)
* A file containing orbital names and labels for retrieving the reference configurations. For later use in matplotlib (with LaTeX), please use this formatting:

```
{
    r"Cu}_{3s}": 0, 
    r"Cu}_{3p_z}": 1, 
    r"Cu}_{3p_x}": 2, 
    r"Cu}_{3p_y}": 3, 
    r"d}_{xz}": 4, 
    r"d}_{yz}": 5, 
    r"d}_{x^2-y^2}": 6, 
    r"d}_{xy}": 7, 
    r"d}_{z^2}": 17, 
    r"4s}": 18,
}
```

* You will also need a file ```ref_csf_threshold.txt``` containing for each root in the calculation, the threshold that should be used to consider a configuration as a "Reference configuration" (and not classify it as Single, Double...)

```
0 0.03
1 0.03
2 0.03
3 0.03
4 0.03
5 0.03
6 0.03
7 0.03
8 0.03
9 0.03
10 0.03
```

* Finally, you need to manually specify in a text file ```ref_def.txt``` the index of the ROHF CSF in ORCA's output and the LMCT loss and gain of electron indices as such:

```
GS 0,0,1
LMCT 15,17
```
Here, the GS is found in the 0th root, as the 0th CSF of that root. The LMCT consists in removing electron from orbital 15 and putting it in orbital 17.

Finally, you can run the program as:

```python $path_csf_analyzer/main.py analyze_degree orca.out orbitals.dict ref_csf_threshold.txt ref_def.txt```

This creates a file called ```excitation_classes.txt``` containing excitation degree information for each root in the wavefunction.

To visualize the results you can run:

```python $path_csf_analyzer/main.py -da visualize excitation_classes.txt```

Which will create a directory ```./plots``` containing bar plots and pie charts of the data.
Remember to add option -da when running ```visualize``` on this type of analysis.

You can also control the threshold for plotting a given contribution in either bar plots of pie charts using respectively options ```-tb 0.01``` or ```-tp 0.01```:

```python $path_csf_analyzer/main.py visualize -tp 0.01 -da excitation_classes.txt```

