import os
import ast 
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

class Wavefunction:
    def __init__(self, orca_output_file=None, orbitals_file=None, thresh_bar=0.01, thresh_pie=0.04):
        self.roots = []
        self.thresh_bar = thresh_bar
        self.thresh_pie = thresh_pie
        self.orbitals_file = orbitals_file
        self.orca_output_file = orca_output_file
        self.wavefunction = {}
        self.CSF = None
        self.orbitals = None

    #Read orbitals names and labels from a predefined file (Warning: ORCA is 0 based)
    def read_orbitals(self):
        # Read the content of the orbitals.dict file
        with open(self.orbitals_file, 'r') as file:
            file_content = file.read()

        # Use ast.literal_eval to parse the content of the file into a dictionary
        self.orbitals = ast.literal_eval(file_content)

    #Parse orca output file and extract CSF and weights
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
    
        pattern = re.compile(r"STATE\s+(\d+).+?\n((?:\s+\d+\.\d+.+?\n)+)")
        matches = pattern.finditer(data)

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

    def analyze(self):
        excitation_classes = defaultdict(float)

        # Mapping from numbers to strings
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

            #Write all contributions to a text file
            with open("excitation_classes.txt", "w") as file:
                for excitation_class, weight in excitation_classes.items():
                    file.write(f"{excitation_class}: {weight}\n")


    def visualize(self, excitation_classes_filename, thresh_pie=0.04, thresh_bar=0.01, save_dir ='./plots'):
        excitation_classes = {}
        with open(excitation_classes_filename, 'r') as file:
            for line in file:
                # Split the line into key and value
                key_str, value_str = line.strip().split(': ')
                
                # Convert value to float
                value = float(value_str)
                
                # Parse the key as a tuple and preserve the original string for later use as Latex
                key_str = key_str.replace('\\\\', '\\').replace("'", "")  
                key_tuple = tuple(key_str.strip('()').split(', '))

                # Add the key-value pair to the dictionary
                excitation_classes[int(key_tuple[0]),key_tuple[1]] = value
        
        # Sort roots and get unique values
        roots = sorted(set(root for root, _ in excitation_classes.keys()))

        # Define colors for the reference CSF and static contributions
        colors = {'Local': 'paleturquoise', 'Others': 'gray'}

        # Create a new color for each unique excitation in the dataset
        excitation_types = set(excitation for _, excitation in excitation_classes.keys())
        # Create a color map with more unique colors by combining several colormaps
        color_list = plt.cm.tab20(np.linspace(0, 1, 20)).tolist() + plt.cm.tab20b(np.linspace(0, 1, 20)).tolist() + plt.cm.tab20c(np.linspace(0, 1, 20)).tolist()
        cm = mcolors.LinearSegmentedColormap.from_list('combined_colormap', color_list, N=len(color_list))
    
        for i, excitation_type in enumerate(sorted(excitation_types - set(colors.keys()))):
            # Check if the excitation type is a 'Root n' label
            if excitation_type.startswith('T'):
                colors[excitation_type] = 'paleturquoise'  # set color for all root types to blue
            else:
                colors[excitation_type] = cm(i)

        # Create the save directory if it doesn't exist
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
            
            # Create an explode list
            #explode = [0.0 if i==0 else 0 for i in range(len(sizes))]
        
            fig, ax = plt.subplots(figsize=(10,5))
            wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=[colors[label] for label in labels], autopct='%1.1f%%')

            #ax.set_title(f'T{root + 1}', fontsize=24)  # font size of title
            for autotext in autotexts:
                autotext.set_bbox(dict(facecolor='none', alpha=0.5, edgecolor='none', boxstyle='round,pad=2'))
            
            #for i, autotext in enumerate(autotexts):
            #    x, y = autotext.get_position()
            #if i % 2 == 0:
            #    autotext.set_position((x*1.1, y*1.1))
            #else:
            #    autotext.set_position((x*1.1, y*1.1))
        
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

            plt.savefig(os.path.join(save_dir, f'T_{root}_loc_pie.png'), bbox_inches='tight',pad_inches = 0,transparent=True)
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
            plt.savefig(os.path.join(save_dir, f'T_{root}_loc_bar.png'), bbox_inches='tight')
            plt.close()
