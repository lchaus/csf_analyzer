import os
import ast 
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

class Wavefunction:
    def __init__(self, orca_output_file=None, orbitals_file=None, thresh_bar=0.01, thresh_pie=0.04, ref_csf_threshold=None, ref_define=None, degree_analysis_vis=False):
        self.thresh_bar = thresh_bar
        self.thresh_pie = thresh_pie
        self.ref_csf_threshold = ref_csf_threshold
        self.ref_csf_threshold_dict = {}
        self.ref_define = ref_define
        self.ground_state_csf = None
        self.lmct_gl = None 
        self.orbitals_file = orbitals_file
        self.orca_output_file = orca_output_file
        self.wavefunction = {}
        self.orbitals = None
        self.degree_analysis_vis = degree_analysis_vis

    ### DATA EXTRACTION ###
    #Read orbitals names and labels from a predefined file 
    def read_orbitals(self):
        # Read the content of the orbitals.dict file (contains the orbitals names and labels)
        with open(self.orbitals_file, 'r') as file:
            file_content = file.read()

        # Use ast.literal_eval to parse the content of the file into a dictionary
        self.orbitals = ast.literal_eval(file_content)

    #Parse orca output file and extract CSFs and weights
    def extract_data(self):
        def parse_csf_string(csf_string):
            parsed_csf = []
            i = 0
            while i < len(csf_string):
                if csf_string[i] == '2' or csf_string[i] == '0':
                    parsed_csf.append(csf_string[i])
                    i += 2
                elif csf_string[i] == '1':
                    parsed_csf.append(csf_string[i:i+2])
                    i += 2
            return parsed_csf

        with open(self.orca_output_file, "r") as file:
            data = file.read()

        # We extract wavefunction written after last iteration in the output file
        pattern_csf = r"                 <<<<<<<<<<<<<<<<<<INITIAL CI STATE CHECK>>>>>>>>>>>>>>>>>>(.*?)                 <<<<<<<<<<<<<<<<<<INITIAL CI STATE CHECK>>>>>>>>>>>>>>>>>>"
        data = re.search(pattern_csf, data, re.DOTALL)

        pattern = re.compile(r"STATE\s+(\d+).+?\n((?:\s+\d+\.\d+.+?\n)+)")
        matches = pattern.finditer(data.group(0))

        self.wavefunction = {}

        for match in matches:
            state_num = int(match.group(1))
            state_lines = match.group(2).strip().split("\n")
    
            self.wavefunction[state_num] = []
            
            for line in state_lines:
                line_data = line.strip().split()
                weight = float(line_data[0])
                csf = " ".join(line_data[4:]).lstrip("=").strip()
                parsed_csf = parse_csf_string(csf) 
                self.wavefunction[state_num].append((weight, parsed_csf))

    ### ANALYSIS ###
    def local_analysis(self):
        excitation_classes = defaultdict(float)

        number_to_string = {1: 'Single', 2: 'Double', 3: 'Triple', 4: 'Quadruple'}

        for root, csfs in self.wavefunction.items():
            ref_weight, ref_csf = csfs[0]
            excitation_classes[(root, f'T{root} ')] = ref_weight

            ref_block_occ = {}
            for block_name, block_indices in self.orbitals.items():
                ref_block_occ[block_name] = sum(int(ref_csf[i][0]) for i in block_indices)

            for weight, curr_csf in csfs:
                if curr_csf == ref_csf:
                    continue

                curr_block_occ = {}
                for block_name, block_indices in self.orbitals.items():
                    curr_block_occ[block_name] = sum(int(curr_csf[i][0]) for i in block_indices)

                source_block = target_block = None
                excitation_degree = 0
                total_difference = 0
                for block_name in self.orbitals.keys():
                    difference = ref_block_occ[block_name] - curr_block_occ[block_name]
                    total_difference += abs(difference)
                    if difference > 0:
                        source_block = block_name
                        excitation_degree += difference
                    elif difference < 0:
                        target_block = block_name

                if source_block is not None and target_block is not None:
                    excitation_classes[(root, r'$\mathrm{' + source_block + r'\rightarrow\mathrm{' + target_block + r'$ ' +  number_to_string.get(excitation_degree, 'Unsupported') )] += weight
                elif total_difference == 0:
                    excitation_classes[(root, 'Local')] += weight

            #Write all contributions to a text file for later use (visualize)
            with open("excitation_classes.txt", "w") as file:
                for excitation_class, weight in excitation_classes.items():
                    file.write(f"{excitation_class}: {weight}\n")
    
    def reference_csf(self):
        # Parse self.ref_csf_threshold (user defined) each root and store it in a dictionary
        self.ref_csf_threshold_dict = {}
        
        with open(self.ref_csf_threshold, 'r') as file:
            for line in file:
                key, value = line.strip().split()
                key = int(key)
                value = float(value)
                self.ref_csf_threshold_dict[str(key)] = value

        # Parse the CSF reference definition for parent CSF (e.g. ROHF) and LMCT CSF (user defined).
        with open(self.ref_define, 'r') as file:
            for line in file:
                key, value = line.strip().split()
                if key == 'GS':
                    root, csf, occ = value.split(',')
                    root = int(root)
                    csf = int(csf)
                    occ = int(occ)
                    self.ground_state_csf = (root, csf, occ)
                elif key == 'LMCT':
                    loss, gain = value.split(',')
                    loss = int(loss)
                    gain = int(gain)
                    self.lmct_gl = (loss, gain)

    def degree_analysis(self):
        def get_excitation_degree(ref_csf, curr_csf):
            ''' A function to compute the excitation degree between two CSFs '''
            electron_loss = 0
            electron_gain = 0
            loss_indices = []
            gain_indices = []

            for i, (ref, curr) in enumerate(zip(ref_csf, curr_csf)):
                if ref != curr:
                    if ref[0] == '2':
                        ref_value = 2
                    elif ref[0] == '1':
                        ref_value = 1
                    else:
                        ref_value = 0
                    if curr[0] == '2':
                        curr_value = 2
                    elif curr[0] == '1':
                        curr_value = 1
                    else:
                        curr_value = 0
                    if ref_value > curr_value:
                        electron_loss += ref_value - curr_value
                        loss_indices.append(i)
                    elif ref_value < curr_value:
                        electron_gain += curr_value - ref_value
                        gain_indices.append(i)
            if electron_loss == electron_gain:
                return electron_loss, (loss_indices, gain_indices)
            else:
                return None

        #Retrieve ground state CSF (user-defined).
        ground_state_csf = self.wavefunction[self.ground_state_csf[0]][self.ground_state_csf[1]][self.ground_state_csf[2]]

        state_configs = defaultdict(list)
        # Loop over all states (roots) and their associated CSFs
        for root, csfs in self.wavefunction.items():
            root_str = str(root)  # convert root number to string to use it as a key
            # Extract the reference CSFs (above a given weight, to take into account heavily multireference states) for the current state (root)
            for weight, curr_csf in csfs:
                if round(weight,2) >= self.ref_csf_threshold_dict[root_str]:  # If the weight is above the threshold, store it as a main CSF
                    excitation_degree, (loss_indices, gain_indices) = get_excitation_degree(ground_state_csf, curr_csf)
                    if excitation_degree == 1 and (loss_indices[0] == self.lmct_gl[0] and gain_indices[0] == self.lmct_gl[1]):
                        state_configs[root].append(("LMCT", loss_indices, gain_indices, weight, curr_csf))
                    else:
                        state_configs[root].append((excitation_degree, loss_indices, gain_indices, weight, curr_csf))

            #Loop over all CSFs for the current state (root)
            #Get all remaining excitations from the Ground state
            for weight, curr_csf in csfs:
                # Ignore the previously stored main CSFs and 'LMCT+Mono' excitations
                if curr_csf in [csf for _, _, _, _, csf in state_configs[root]]:
                    continue
                excitation_degree, (loss_indices, gain_indices) = get_excitation_degree(ground_state_csf, curr_csf)
                state_configs[root].append((excitation_degree, loss_indices, gain_indices, weight, curr_csf))

            summary = defaultdict(lambda: defaultdict(float))  # Initialize a dictionary with default values as 0.0
            inverse_original_orbitals = {v: k for k, v in self.orbitals.items()}  # To map back from indices to orbital names
    
            #Process the information and store in a dictionary (summary)
            for root, configs in state_configs.items():
                # Loop over all configurations for the current root
                summary[root]['Other'] = 1 - sum([config[3] for config in configs])
                for config in configs:
                    #Track the LMCT in all roots
                    if config[0] == 'LMCT':
                        summary[root]['LMCT'] += config[3]
                    #Classify the reference CSFs (single and double excitations from the ground state with a weight superior to the threshold)
                    elif config[0] in [1, 2] and config[3] >= self.ref_csf_threshold_dict[str(root)]:
                        if config[0] == 1:
                            for loss_indice, gain_indice in zip(config[1], config[2]):
                                key = r'$\\mathrm{'+ inverse_original_orbitals[loss_indice]+ r'\\rightarrow\\mathrm{' + inverse_original_orbitals[gain_indice] + r'$'
                                summary[root][key] += config[3]
                        else:  # case of double excitations (because some reference configurations are double excitations)
                            key = r'$' + ', '.join([r'\mathrm{%s' % inverse_original_orbitals[i] for i in config[1]]) + r' \rightarrow ' + ', '.join([r'\mathrm{%s' % inverse_original_orbitals[i] for i in config[2]]) + r'$'
                            summary[root][key] += config[3]
                    #Classify remaining excitations
                    elif 0 < config[0] <= 4:
                        summary[root][{1: 'Single', 2: 'Double', 3: 'Triple', 4: 'Quadruple'}[config[0]]] += config[3]
                        #Track the LMCT+Mono excitations in an additionnal category
                        if config[0] == 2:
                            if (self.lmct_gl[0] in config[1]) and (self.lmct_gl[1] in config[2]):
                                summary[root]['LMCT+Single'] += config[3]
                    #Track ROHF CSF
                    elif config[0] == 0:
                        summary[root]['ROHF'] += config[3]
                    #Sum the remaining CSFs in an 'others' category
                    else:
                        summary[root]['Other'] += config[3]

            # Write the data to the text file
            with open("excitation_classes.txt", 'w') as file:
                for root, value in summary.items():
                    for key, value in value.items():
                        file.write(f"({root}, '{key}'): {value}\n")

    ### VISUALIZATION ###
    def visualize(self, excitation_classes_filename, thresh_pie=0.04, thresh_bar=0.01, save_dir ='./plots'):
        excitation_classes = {}
        if not self.degree_analysis_vis:
            with open(excitation_classes_filename, 'r') as file:
                for line in file:
                    key_str, value_str = line.strip().split(': ')
                    value = float(value_str)
                    key_str = key_str.replace('\\\\', '\\').replace("'", "")  
                    key_tuple = tuple(key_str.strip('()').split(', '))
                    excitation_classes[int(key_tuple[0]),key_tuple[1]] = value
        elif self.degree_analysis_vis:
            with open(excitation_classes_filename, 'r') as file:
                for line in file:
                    key_str, value_str = line.strip().split(': ')
                    key_str = key_str.replace('\\\\', '\\').replace("'", "") 
                    value = float(value_str)
                    key_tuple = (int(key_str.strip('()').split(', ')[0]), key_str.strip('()').split(', ',1)[1])
                    excitation_classes[key_tuple] = value
       
        # Sort roots and get unique values
        roots = sorted(set(root for root, _ in excitation_classes.keys()))

        # Define colors for the reference CSF and static contributions
        colors = {'Local': 'paleturquoise', 'Others': 'gray', 'ROHF': 'plum', 'LMCT': 'gold', 'LMCT+Single': 'tomato', 'Other': 'gray'}

        # Create a new color for each unique excitation in the dataset
        excitation_types = set(excitation for _, excitation in excitation_classes.keys())

        # Create a color map with more unique colors by combining several colormaps
        color_list = plt.cm.tab20(np.linspace(0, 1, 20)).tolist() + plt.cm.Pastel1(np.linspace(0, 1, 9)).tolist() + plt.cm.tab20c(np.linspace(0, 1, 20)).tolist()  + plt.cm.Pastel2(np.linspace(0, 1, 8)).tolist()
        cm = mcolors.LinearSegmentedColormap.from_list('combined_colormap', color_list, N=len(color_list))

        for i, excitation_type in enumerate(sorted(excitation_types - set(colors.keys()))):
            # Check if the excitation type is a 'Root n' label
            if excitation_type.startswith('T'):
                colors[excitation_type] = 'paleturquoise'  # set color for all root types to blue
            else:
                colors[excitation_type] = cm(i)

        os.makedirs(save_dir, exist_ok=True)

        # Iterate over roots to create pie charts
        for root in roots:
            # Filter out excitations with weight below threshold
            labels, sizes = zip(*[(excitation, excitation_classes[(root, excitation)]) 
                                  for _, excitation in excitation_classes.keys() 
                                  if _ == root and excitation_classes[(root, excitation)] >= thresh_pie])
        
            # Calculate the weight of the neglected contributions
            neglected_weight = 1 - sum(sizes)
        
            # Add the neglected contributions as "Others"
            if neglected_weight > 0:
                labels += ('Others',)
                sizes += (neglected_weight,)
        
            fig, ax = plt.subplots(figsize=(10,5))
            wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=[colors[label] for label in labels], autopct='%1.1f%%')

            #ax.set_title(f'T{root + 1}', fontsize=24)  # font size of title
            for autotext in autotexts:
                autotext.set_bbox(dict(facecolor='none', alpha=0.5, edgecolor='none', boxstyle='round,pad=2'))
            
            # Add a black border to the first slice
            wedges[0].set_edgecolor('black')
            wedges[0].set_linewidth(1)  # Adjust as needed
            wedges[1].set_edgecolor('black')
            wedges[1].set_linewidth(1)  # Adjust as needed
        
            # Make the pie chart size more consistent by expanding the axes
            ax.set_aspect('equal')  # Equal aspect ratio ensures the pie chart is circular
            ax.axis('off')  # Hide the axes
        
            plt.setp(autotexts, size=16) # font size of pie chart values
            plt.setp(texts, size=16) # font size of pie chart labels
        
            # Save the figure, convert root to string before concatenating
            plt.margins(0,0)

            plt.savefig(os.path.join(save_dir, f'T_{root}_pie.png'), bbox_inches='tight',pad_inches = 0,transparent=True)
            plt.close()

        # Find the maximum weight value across all roots
        max_weight = max(value for value in excitation_classes.values() if value >= thresh_bar)

        # Iterate over roots to create bar charts
        for root in roots:
            # Filter out excitations with weight below threshold
            items = [(excitation, excitation_classes[(root, excitation)]) 
                      for _, excitation in excitation_classes.keys() 
                      if _ == root and excitation_classes[(root, excitation)] >= thresh_bar]

            # Sort items by size
            items.sort(key=lambda x: x[1])

            # Calculate the weight of the neglected contributions
            neglected_weight = 1 - sum(size for _, size in items)

            # Add the neglected contributions as "Others"
            if neglected_weight > 0:
                items.append(('Others', neglected_weight))

            labels, sizes = zip(*items)

            fig, ax = plt.subplots(figsize=(10,6))
            y_pos = np.arange(len(labels))
            ax.barh(y_pos, sizes, color=[colors[label] for label in labels])
            ax.set_yticks(y_pos)
            ax.set_yticklabels(labels, fontsize=12)  # font size of labels
            ax.set_xlabel('Weight', fontsize=12)  # font size of x label
            ax.set_title(f'T{root}', fontsize=14)  # font size of title

            # Set the x-axis limit to be the same for all plots
            ax.set_xlim(0, max_weight + 0.1)  # Added a small buffer for annotations

            # Annotate weights on the side of the bars
            for i, v in enumerate(sizes):
                ax.text(v + 0.01, i, f'{v:.3f}', color='black', fontsize=10)

            # Save the figure
            plt.savefig(os.path.join(save_dir, f'T_{root}_bar.png'), bbox_inches='tight')
            plt.close()
