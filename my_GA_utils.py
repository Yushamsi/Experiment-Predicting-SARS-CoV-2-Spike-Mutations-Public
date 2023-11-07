
# # Import Statements

# 
import random
import numpy as np
import pandas as pd
import pygad
import subprocess
import os
import blosum as bl
import random
import math
import re
import time

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align import AlignInfo
from Bio import AlignIO

from pymsaviz import MsaViz
from IPython.display import display

import unittest



# 
# # Supplementary Materials
# 

# 
# ### Data Extraction

# 
def check_and_create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f'Created folder [{directory}] for saving data.')


# 



def extract_sequences(input_path):
    sequences = []  # List to store the extracted sequences

    try:
        if os.path.isdir(input_path):  # If the input path is a directory
            # Iterate over files in the directory
            for file_name in os.listdir(input_path):
                file_path = os.path.join(input_path, file_name)

                # Check if the file is a FASTA file
                if file_name.endswith(".fasta") or file_name.endswith(".fa"):
                    # Iterate over sequences in the FASTA file
                    for record in SeqIO.parse(file_path, "fasta"):
                        sequence = str(record.seq)  # Convert sequence to string

                        # Process the sequence as needed
                        sequences.append(sequence)

        elif os.path.isfile(input_path):  # If the input path is a file
            # Read sequences from the fasta file and convert to strings
            for record in SeqIO.parse(input_path, "fasta"):
                sequences.append(str(record.seq))

        else:  # If the input path is neither a file nor a directory
            print("The provided input_path is neither a file nor a directory.")

    except Exception as e:
        print("An error occurred while reading the FASTA file(s). Please ensure that the provided path is either a folder containing FASTA files or a path to a single FASTA file.")
        print(f"Error details: {str(e)}")
        return None

    print(f"Sequences extracted from {input_path}")

    return sequences


# 


def write_list_to_fasta(sequences_list, output_fasta_file):
    # Convert string sequences to SeqRecords
    seq_records = [SeqRecord(Seq(seq), id=f"seq{i+1}") for i, seq in enumerate(sequences_list)]

    # Write SeqRecords to the output file
    with open(output_fasta_file, "w") as output_handle:
        SeqIO.write(seq_records, output_handle, "fasta")

    print(f"Sequences written to {output_fasta_file}")

# 
# ### Conversion of Data

# 
def map_numbers_to_amino_acids(numbers):
    
    amino_acids = {
        1: 'A', 2: 'R', 3: 'N', 4: 'D', 5: 'C', 6: 'Q', 7: 'E', 8: 'G',
        9: 'H', 10: 'I', 11: 'L', 12: 'K', 13: 'M', 14: 'F', 15: 'P',
        16: 'S', 17: 'T', 18: 'W', 19: 'Y', 20: 'V'
    }

    sequence = []
    for number in numbers:
        if isinstance(number, np.ndarray) and number.size > 1:
            sequence.append('X')  # 'X' represents an unknown amino acid for arrays of size > 1
        else:
            if isinstance(number, np.generic):
                number = number.item()  # Convert NumPy scalar to Python scalar
            number = int(number)  # Convert to int
            if number in amino_acids:
                sequence.append(amino_acids[number])
            else:
                print("ERROR map_numbers_to_amino_acids cannot convert number", number)
                sequence.append('X')  # 'X' represents an unknown amino acid

    sequence = ''.join(sequence)

    return sequence




# 
def map_amino_acids_to_numbers(sequence):

    amino_acids = {
        'A': 1, 'R': 2, 'N': 3, 'D': 4, 'C': 5, 'Q': 6, 'E': 7, 'G': 8,
        'H': 9, 'I': 10, 'L': 11, 'K': 12, 'M': 13, 'F': 14, 'P': 15,
        'S': 16, 'T': 17, 'W': 18, 'Y': 19, 'V': 20
    }

    numbers = []
    for amino_acid in sequence:
        if amino_acid in amino_acids:
            numbers.append(amino_acids[amino_acid])
        else:
            # You can handle the case of an unknown amino acid here, if needed.
            print("ERROR map_amino_acids_to_numbers cannot convert amino acid", amino_acid)
            numbers.append(0)

    return numbers




# 
# ### Preprocessing

# 
def convert_sequences_to_numbers(extracted_sequences):
    # Initialize an empty 2D array
    amino_acid_numbers = []

    # Iterate over each amino acid sequence
    for sequence in extracted_sequences:
        # Convert the amino acid sequence to numbers using the map_amino_acids_to_numbers function
        sequence_numbers = map_amino_acids_to_numbers(sequence)

        # Append the sequence numbers to the 2D array
        amino_acid_numbers.append(sequence_numbers)

    return amino_acid_numbers

# 
def preprocess_sequences(folder_path):
    print('Extracting Sequences...')
    extracted_sequences = extract_sequences(folder_path)

    amino_acid_sequences = extracted_sequences

    print('Processing Sequences...')
    amino_acid_numbers = convert_sequences_to_numbers(amino_acid_sequences)

    return amino_acid_numbers

# 
# # Genetic Algorithim Supplementary Materials
# 

# 
# ### Fitness Function

# 
def run_spikepro(sequence):
    # Set directory for SpikePro
    dir_path = "SpikePro_Requirements"

    # Compile the C++ code
    compile_command = "c++ SpikePro.cpp edlib/src/edlib.cpp CSVparser.cpp -o SpikePro -I edlib/include/ -std=c++11"
    subprocess.run(compile_command, shell=True, cwd=dir_path)

    # Save the sequence to a temporary file
    sequence_file = os.path.join(dir_path, "temp.fasta")
    with open(sequence_file, "w") as file:
        file.write(sequence)

    # Run the C++ code for the temporary sequence file
    run_command = f"./SpikePro temp.fasta go"
    result = subprocess.run(run_command, shell=True, capture_output=True, text=True, cwd=dir_path)

    # Extract the output from the subprocess result
    output = result.stdout

    # Parse the output and extract the fitness values
    fitness = 0.0
    lines = output.split('\n')
    for line in lines:
        if line.startswith("Spike protein predicted fitness:"):
            fitness = float(line.split(": Φ = ")[1].split(" ")[0])
            if not fitness:
                fitness = 0.0

    # Remove the temporary sequence file
    os.remove(sequence_file)

    return fitness


# 
# ### BLOSUM90

# 
# BLOSUM90 MATRIX FOR MUTATIONS
def sample_amino_acid_blosum(amino_acid):
    matrix = bl.BLOSUM(90)
    val = matrix[amino_acid]

    keys_to_remove = ['B', 'J', 'Z', 'X', '*']
    for key in keys_to_remove:
        if key in val:
            del val[key]

    # Compute the exponent of each value in the dictionary
    exponents = {key: math.exp(value) for key, value in val.items()}

    # Compute the sum of all exponentiated values
    sum_exponents = sum(exponents.values())

    # Convert values to probabilities that add up to 1
    probabilities = {key: value / sum_exponents for key, value in exponents.items()}

    # Randomly sample a value based on the probabilities
    sampled_value = random.choices(list(probabilities.keys()), weights=list(probabilities.values()))[0]

    return sampled_value



# 
# ### PyGAD Genetic Algorithim Specifics
# 



# Define fitness function using the run_spikepro function for PyGAD
def fitness_func(ga_instance, solution, solution_idx):
    sequence = map_numbers_to_amino_acids(solution)
    
    fitness_value = run_spikepro(sequence)

    return fitness_value


# Define mutation function using the BLOSUM function for PyGAD
def mutate_func(solutions, solution_idx):
    mutated_solutions = []
    for idx, solution in enumerate(solutions):
        max_mutations = 10
        mutated_solution = solution.copy() # No need for solution[0] if solution is not nested
        num_mutations = 0
        indexes = list(range(len(solution)))
        random.shuffle(indexes)  # Randomize the order of indexes

        # Iterate over each element (amino acid) in a randomized order
        for i in indexes:
            amino_acid = map_numbers_to_amino_acids([mutated_solution[i]])

            # Get the most likely amino acid mutations based on the BLOSUM90 matrix
            most_likely = sample_amino_acid_blosum(amino_acid)
            most_likely_number = map_amino_acids_to_numbers(most_likely)   

            # Mutate the amino acid in the solution if the maximum mutation count is not reached
            if num_mutations < max_mutations:
                mutated_solution[i] = most_likely_number[0]
                num_mutations += 1
                # print('Num Mutations:', num_mutations)
            else:
                break

        mutated_solutions.append(mutated_solution)

    output = np.array(mutated_solutions)

    return output






# # Genetic Algorithim


def run_ga_multiple_times(initial_population, num_generations, num_parents_mating, fitness_func, mutation_type, parent_selection_type, K_tournament, crossover_type, keep_parents, keep_elitism, save_solutions, stop_criteria, random_seed, suppress_warnings, num_runs, save_and_display_images=True):
    best_solutions = []
    best_amino_acid_sequences = []
    for i in range(num_runs):
        ga_instance = pygad.GA(initial_population=initial_population, 
                               num_generations=num_generations, 
                               num_parents_mating=num_parents_mating, 
                               fitness_func=fitness_func, 
                               mutation_type=mutation_type,
                               parent_selection_type=parent_selection_type, 
                               K_tournament=K_tournament,
                               crossover_type=crossover_type,
                               keep_parents=keep_parents,
                               keep_elitism=keep_elitism,
                               save_solutions=save_solutions,
                               stop_criteria=stop_criteria,
                               random_seed=random_seed,
                               suppress_warnings=suppress_warnings)
        ga_instance.run()

        if save_and_display_images:
            # Save plots
            os.makedirs('GA Performance Images', exist_ok=True)
            ga_instance.plot_fitness(plot_type="scatter", title=f'Run {i+1}: Generation vs Fitness', save_dir=f'GA Performance Images/Generation_vs_Fitness_Run_{i+1}')
            ga_instance.plot_new_solution_rate(title=f'Run {i+1}: Generation vs New Solution Rate', save_dir=f'GA Performance Images/Generation_vs_New_Solution_Rate_Run_{i+1}')

        # Get best solution
        solution, solution_fitness, solution_idx = ga_instance.best_solution()
        amino_acid_sequence = map_numbers_to_amino_acids(solution)
        best_solutions.append((i+1, amino_acid_sequence, solution_fitness, solution_idx))
        best_amino_acid_sequences.append(amino_acid_sequence)

    # Write the best sequences to a fasta file
    write_list_to_fasta(best_amino_acid_sequences, "GA Saved Results/best_sequences.fasta")

    print('Best sequences saved to GA Saved Results/best_sequences.fasta.')

    return best_amino_acid_sequences




# # Displaying Results 


# ### MSA File Generation, Visualisation and Data Save


def create_msa_files(sequences, output_file, consensus_file):
    # The specific sequence to be added at the beginning
    specific_seq = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCGLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
    specific_record = SeqRecord(Seq(specific_seq), id="Wuhan-Hu Spike Glycoprotein")

    # Convert list of sequences into SeqRecord objects
    seq_records = [SeqRecord(Seq(seq), id=f"seq{i+1}") for i, seq in enumerate(sequences)]

    # Add the specific sequence at the beginning
    seq_records.insert(0, specific_record)

    # Write the sequences to a temporary file
    SeqIO.write(seq_records, "temp.fa", "fasta")

    # Perform a multiple sequence alignment with Clustal Omega
    clustalomega_cline = ClustalOmegaCommandline("Clustal Omega/clustalo", infile="temp.fa", outfile=output_file, verbose=True, auto=True, force=True)
    clustalomega_cline()

    # Read the alignment
    msa = AlignIO.read(output_file, "fasta")

    # Compute the consensus sequence
    summary_align = AlignInfo.SummaryInfo(msa)
    consensus = summary_align.dumb_consensus()

    # Create a new SeqRecord for the consensus sequence
    consensus_record = SeqRecord(Seq(str(consensus)), id="consensus")

    # Write the consensus sequence to a separate file
    SeqIO.write(consensus_record, consensus_file, "fasta")











def visualize_alignment(file_path):
    mv = MsaViz(file_path, wrap_length=60, show_count=True, show_consensus=True)

    # Define the regions and their ranges
    regions = {
        "S1 subunit NTD": (13, 304, "red"),
        "S1 subunit RBD": (319, 541, "green"),
        "RBD Receptor Binding Motif": (438, 508, "blue"),
        "SD-1 & SD-2 sub-domains, S1/S2 cleavage region, S2 fusion subunit": (543, 1208, "cyan"),
        "S1/S2 cleavage region": (672, 709, "magenta"),
        "S2 Fusion Peptide (FP-1)": (798, 806, "yellow"),
        "S2 Internal Fusion Peptide  (FP-2)": (816, 833, "black"),
        "S2 Heptad-repeat-1": (918, 983, "purple"),
        "S2 Heptad-repeat-2": (1162, 1203, "orange"),
        "Transmembrane domain": (1214, 1234, "pink"),
        "S2 subunit, intra-virion / Cytoplasmic Tail (CT)": (1233, 1273, "brown"),
    }

    # Add text annotations for each region
    for region, (start, end, color) in regions.items():
        mv.add_text_annotation((start, end), region, text_color=color, range_color=color)

    output_file = file_path.rsplit('.', 1)[0] + "_with_colored_annotations.png"
    mv.savefig(output_file)
    print(f"Saved visualization to {output_file}")









# ### SpikePro Data Extraction



def run_spikepro_full(sequence):
    # Set directory for SpikePro
    dir_path = "SpikePro_Requirements"

    # Compile the C++ code
    compile_command = "c++ SpikePro.cpp edlib/src/edlib.cpp CSVparser.cpp -o SpikePro -I edlib/include/ -std=c++11"
    subprocess.run(compile_command, shell=True, cwd=dir_path)

    # Save the sequence to a temporary file
    sequence_file = os.path.join(dir_path, "temp.fasta")
    with open(sequence_file, "w") as file:
        file.write(sequence)

    # Run the C++ code for the temporary sequence file
    run_command = f"./SpikePro temp.fasta go"
    result = subprocess.run(run_command, shell=True, capture_output=True, text=True, cwd=dir_path)

    # Extract the output from the subprocess result
    output = result.stdout

    # Parse the output and extract the fitness values
    phi = float(re.findall(r"Φ = (\d+\.\d+)", output)[0])
    stab = float(re.findall(r"stab=(\d+\.\d+)", output)[0])
    ace2 = float(re.findall(r"ACE2=(\d+\.\d+)", output)[0])
    nab = float(re.findall(r"nAb=(\d+\.\d+)", output)[0])
    total_variants = int(re.findall(r"Total number of amino acid variants: (\d+)", output)[0])

    # Extracting all the variants
    variants = re.findall(r"(Variant [\w\d]+, occurrence in GISAID = \d+\.\d+%)", output)
    
    # Remove the temporary sequence file
    os.remove(sequence_file)

    return phi, stab, ace2, nab, total_variants, variants





def run_spikepro_best_sequences(sequences):
    # Initialize an empty list to store each sequence's data
    all_data = []

    # Initialize an empty dictionary to store each variant's data
    all_variants = {}

    # Loop over the list of sequences
    for i, sequence in enumerate(sequences, start=1):
        # Use the run_spikepro_full function to get the sequence's data
        phi, stab, ace2, nab, total_variants, variants = run_spikepro_full(sequence)
        
        # Store the sequence's data in a dictionary
        sequence_data = {
            'Seq Name': f'Seq {i}',
            'Sequence': sequence,
            'Phi': phi,
            'Stab': stab,
            'ACE2': ace2,
            'nAb': nab,
            'Total Variants': total_variants
        }

        # Append the dictionary to the list of all data
        all_data.append(sequence_data)

        # Loop over the list of variants
        for variant in variants:
            # If the variant is not in the dictionary, add it
            if variant not in all_variants:
                all_variants[variant] = {'Count': 1, 'Sequences': [f'Seq {i}']}
            # If the variant is already in the dictionary, increment its count and add the sequence to its list
            else:
                all_variants[variant]['Count'] += 1
                all_variants[variant]['Sequences'].append(f'Seq {i}')

    # Convert the list of dictionaries to a DataFrame
    df_sequences = pd.DataFrame(all_data)

    # Convert the dictionary to a DataFrame
    df_variants = pd.DataFrame.from_dict(all_variants, orient='index')

    return df_sequences, df_variants




def analyze_and_save_dataframes(initial_sequences, final_sequences):
    # Set the base path for the result files
    base_path = 'GA Saved Results'

    # Check if the directory exists, if not, create it
    if not os.path.exists(base_path):
        os.makedirs(base_path)
        print(f"Folder '{base_path}' has been created.")
        
    # Analyze the initial population
    print("Analyzing the initial population...")
    
    # Run spikepro on the sequences to get the two dataframes
    df_initial_population_sequences, df_initial_population_variants = run_spikepro_best_sequences(initial_sequences)

    # Save the data frames to csv files
    df_initial_population_sequences.to_csv(os.path.join(base_path, 'initial_population_sequences.csv'), index=True)
    print("Initial population sequences data saved to 'initial_population_sequences.csv'.")
    
    df_initial_population_variants.to_csv(os.path.join(base_path, 'initial_population_variants.csv'), index=True)
    print("Initial population variants data saved to 'initial_population_variants.csv'.")

    print("\nAnalyzing generated sequences...")

    # Call the run_spikepro_best_sequences function with the sequences
    df_sequences, df_variants = run_spikepro_best_sequences(final_sequences)

    # Save the data frames to csv files
    df_sequences.to_csv(os.path.join(base_path, 'sequences.csv'), index=True)
    print("Generated sequences data saved to 'sequences.csv'.")
    
    df_variants.to_csv(os.path.join(base_path, 'variants.csv'), index=True)
    print("Generated variant data saved to 'variants.csv'.")

    # Novelty analysis for variants
    print("Performing novelty analysis for all variants...")

    # Initialize a new column to indicate whether a variant exists in the initial population
    df_variants['Exists in Initial Population'] = df_variants.index.isin(df_initial_population_variants.index)

    # Add new column where variants exist in the initial population
    df_variants['Initial Population Sequences that Variant Exists In'] = df_variants.apply(lambda row: df_initial_population_variants.loc[row.name, 'Sequences'] if row['Exists in Initial Population'] else None, axis=1)

    df_variants_with_novelty_analysis = df_variants

    # Save the data frame to csv file
    df_variants_with_novelty_analysis.to_csv(os.path.join(base_path, 'variants_with_novelty_analysis.csv'), index=True)
    print("Variant data with novelty analysis saved to 'variants_with_novelty_analysis.csv'.")

    # Display the DataFrames
    print("\nInitial Population Sequences")
    display(df_initial_population_sequences)
    display(df_initial_population_variants)
    print("\nOutput Sequences")
    display(df_sequences)
    display(df_variants_with_novelty_analysis)





# # Running Everything
# 


def run_all_ga(folder_path, num_generations, num_parents_mating, fitness_func, mutation_type, parent_selection_type, K_tournament, crossover_type, keep_parents, keep_elitism, save_solutions, stop_criteria, random_seed, suppress_warnings, num_runs, save_and_display_images=True):
    start_time = time.time() # Record the start time

    print('Extracting and Processing Data...')
    # Process the sequences
    initial_population = preprocess_sequences(folder_path)
    print('Data Extraction and Processing Successful. Elapsed time:', time.time() - start_time, 'seconds.')

    print("\nChecking Directory for Saving Results...")
    check_and_create_directory('GA Saved Results')

    print('\nRunning Genetic Algorithm...')
    # Run the genetic algorithm multiple times
    best_amino_acid_sequences = run_ga_multiple_times(initial_population=initial_population, 
                                                      num_generations=num_generations, 
                                                      num_parents_mating=num_parents_mating, 
                                                      fitness_func=fitness_func,
                                                      mutation_type=mutation_type,
                                                      parent_selection_type=parent_selection_type, 
                                                      K_tournament=K_tournament,
                                                      crossover_type=crossover_type,
                                                      keep_parents=keep_parents,
                                                      keep_elitism=keep_elitism,
                                                      save_solutions=save_solutions,
                                                      stop_criteria=stop_criteria,
                                                      random_seed=random_seed,
                                                      suppress_warnings=suppress_warnings,
                                                      num_runs=num_runs, 
                                                      save_and_display_images=save_and_display_images)
    print('Genetic Algorithm Run Completed.')
    
    # Calculate total elapsed time in seconds
    total_seconds = time.time() - start_time
    # Convert to minutes and seconds
    minutes, seconds = divmod(total_seconds, 60)
    
    print('\nGenetic Algorithm Complete. Elapsed time: {:.2f} seconds ({:.0f} minutes and {:.2f} seconds).'.format(total_seconds, minutes, seconds))
    print('Genetic Algorithm Complete Please Run Results Below')

    return best_amino_acid_sequences




def run_all_results(folder_path, fasta_file_path):
    start_time = time.time() # Record the start time

    # Check and create directory if it doesn't exist
    check_and_create_directory('GA Saved Results')

    # Extract sequences from the folder
    initial_sequences = extract_sequences(folder_path)

    # Extract sequences from the fasta file
    final_sequences = extract_sequences(fasta_file_path)

    print('\nCalculating Results...')
    # Analyze and save data frames
    analyze_and_save_dataframes(initial_sequences, final_sequences)

    # Set the base path for the result files
    base_path = 'GA Saved Results'

    # Provide the file paths for the output alignment and consensus files
    output_alignment_file = os.path.join(base_path, "my_alignment.fa")
    output_consensus_file = os.path.join(base_path, "consensus.fa")

    # Call the create_msa_files function with the sequences and output file paths
    print("\nCreating Multiple Sequence Alignment Files...")
    create_msa_files(final_sequences, output_alignment_file, output_consensus_file)
    print("MSA files created and saved.")

    # Call the visualize_alignment function with the output alignment file path that was just created
    print("\nVisualizing Multiple Sequence Alignment Files...")
    visualize_alignment(output_alignment_file)
    print("MSA visualization complete.")

    # Calculate total elapsed time in seconds
    total_seconds = time.time() - start_time
    # Convert to minutes and seconds
    minutes, seconds = divmod(total_seconds, 60)

    print(f'\nResults Calculated. Total elapsed time: {total_seconds:.2f} seconds ({minutes:.0f} minutes and {seconds:.2f} seconds).')



# # Testing


class TestGeneticAlgorithmFunctions1(unittest.TestCase):
    def test_sample_amino_acid_blosum(self):
        amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        sampled_value = sample_amino_acid_blosum('A')
        self.assertIn(sampled_value, amino_acids, f"Unexpected value {sampled_value}")

    def test_run_spikepro(self):
        sequence = "LEKTTELLFLVMFLLTTKRTMFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCC"
        fitness = run_spikepro(sequence)
        self.assertIsInstance(fitness, float, f"Unexpected type {type(fitness)}")

    def test_map_amino_acids_to_numbers(self):
        sequence = 'ARNDC'
        expected_output = [1, 2, 3, 4, 5]
        self.assertEqual(map_amino_acids_to_numbers(sequence), expected_output, f"Unexpected output {map_amino_acids_to_numbers(sequence)}")

    def test_map_numbers_to_amino_acids(self):
        sequence = [1, 2, 3, 4, 5]
        expected_output = 'ARNDC'
        self.assertEqual(map_numbers_to_amino_acids(sequence), expected_output, f"Unexpected output {map_numbers_to_amino_acids(sequence)}")

# Run the tests
# unittest.main(argv=[''], exit=False)


class TestGeneticAlgorithmFunctions(unittest.TestCase):
    def test_extract_sequences(self):
        # Provide a valid folder_path
        folder_path = 'SpikePro_Requirements/Test Sequences'
        extracted_sequences = extract_sequences(folder_path)

        # Check if the function returns a list
        self.assertIsInstance(extracted_sequences, list, f"Unexpected type {type(extracted_sequences)}")

        # Check if the function returns the expected number of sequences
        self.assertEqual(len(extracted_sequences), 3, f"Unexpected number of sequences {len(extracted_sequences)}")

        # Check if each sequence is a string of capital letters
        for sequence in extracted_sequences:
            self.assertIsInstance(sequence, str, f"Unexpected type {type(sequence)} for sequence {sequence}")
            self.assertTrue(sequence.isupper(), f"Sequence {sequence} is not all uppercase")

        # Check if each sequence only contains valid amino acid letters
        valid_amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
        for sequence in extracted_sequences:
            for letter in sequence:
                self.assertIn(letter, valid_amino_acids, f"Invalid amino acid {letter} in sequence {sequence}")

    def test_convert_sequences_to_numbers(self):
        # Provide a list of valid sequences
        sequences = ["LEKTTELLFLVMFLLTTKRTMFVFLVLLPLVSSQCV", "LEKTTELLFLVMFLLTTKRTMFVFLVLLPLVSSQCV"]
        numbers = convert_sequences_to_numbers(sequences)

        # Check if the function returns a list
        self.assertIsInstance(numbers, list, f"Unexpected type {type(numbers)}")

        # Check if the function returns the expected number of sequences
        self.assertEqual(len(numbers), 2, f"Unexpected number of sequences {len(numbers)}")

        # Check if the function returns the correct conversion for each sequence
        expected_output = [[11, 7, 12, 17, 17, 7, 11, 11, 14, 11, 20, 13, 14, 11, 11, 17, 17, 12, 2, 17, 13, 14, 20, 14, 11, 20, 11, 11, 15, 11, 20, 16, 16, 6, 5, 20], [11, 7, 12, 17, 17, 7, 11, 11, 14, 11, 20, 13, 14, 11, 11, 17, 17, 12, 2, 17, 13, 14, 20, 14, 11, 20, 11, 11, 15, 11, 20, 16, 16, 6, 5, 20]]
        self.assertEqual(numbers, expected_output, f"Unexpected output {numbers}")

    def test_preprocess_sequences(self):
        # Provide a valid folder_path
        folder_path = 'SpikePro_Requirements/Test Sequences'
        processed_sequences = preprocess_sequences(folder_path)

        # Check if the function returns a list
        self.assertIsInstance(processed_sequences, list, f"Unexpected type {type(processed_sequences)}")

        # Check if the function returns the expected number of sequences
        self.assertEqual(len(processed_sequences), 3, f"Unexpected number of sequences {len(processed_sequences)}")

        # Check if each sequence is a list of numbers
        for sequence in processed_sequences:
            self.assertIsInstance(sequence, list, f"Unexpected type {type(sequence)} for sequence {sequence}")
            for number in sequence:
                self.assertIsInstance(number, int, f"Unexpected type {type(number)} for number {number} in sequence {sequence}")

# Run the tests
# unittest.main(argv=[''], exit=False)

##### Further Analysis ####
# Update map_residue_to_region to account for mutations falling under multiple regions
# ## Analysis of SARS-CoV-2 Spike Protein Sequences


import pandas as pd


import ast


# Update map_residue_to_region to account for mutations falling under multiple regions
def map_residue_to_multiple_regions(residue_number):
    regions_list = []
    
    if 13 <= residue_number <= 304:
        regions_list.append('S1 subunit NTD')
    if 319 <= residue_number <= 541:
        regions_list.append('S1 subunit RBD')
    if 438 <= residue_number <= 508:
        regions_list.append('RBD Receptor Binding Motif')
    if 543 <= residue_number <= 1208:
        regions_list.append('SD-1 & SD-2 sub-domains, S1/S2 cleavage region, S2 fusion subunit')
    if 672 <= residue_number <= 709:
        regions_list.append('S1/S2 cleavage region')
    if 798 <= residue_number <= 806:
        regions_list.append('S2 Fusion Peptide (FP-1)')
    if 816 <= residue_number <= 833:
        regions_list.append('S2 Internal Fusion Peptide (FP-2)')
    if 918 <= residue_number <= 983:
        regions_list.append('S2 Heptad-repeat-1')
    if 1162 <= residue_number <= 1203:
        regions_list.append('S2 Heptad-repeat-2')
    if 1214 <= residue_number <= 1234:
        regions_list.append('Transmembrane domain')
    if 1233 <= residue_number <= 1273:
        regions_list.append('S2 subunit, intra-virion / Cytoplasmic Tail (CT)')
    if not regions_list:
        regions_list.append('Unknown')
    
    return regions_list


def further_analyze_sequences(file_path_sequences, file_path_variants):
    # Load sequences.csv
    sequences_df = pd.read_csv(file_path_sequences)
    
    # Load variants_with_novelty_analysis.csv
    variants_df = pd.read_csv(file_path_variants)
    
    # Extract mutation name and occurrence in GISAID
    variants_df['Mutation'] = variants_df['Unnamed: 0'].apply(lambda x: x.split(',')[0].split(' ')[1])
    variants_df['Residue Number'] = variants_df['Mutation'].apply(lambda x: int(''.join(filter(str.isdigit, x))))
    
    # Convert the string representation of list to actual list
    variants_df['Sequences'] = variants_df['Sequences'].apply(ast.literal_eval)
    
    # Map residue number to spike protein regions
    variants_df['Mapped Regions'] = variants_df['Residue Number'].apply(map_residue_to_multiple_regions)
    
    # Create DataFrame for all regions
    all_regions = [
        'S1 subunit NTD', 'S1 subunit RBD', 'RBD Receptor Binding Motif',
        'SD-1 & SD-2 sub-domains, S1/S2 cleavage region, S2 fusion subunit',
        'S1/S2 cleavage region', 'S2 Fusion Peptide (FP-1)',
        'S2 Internal Fusion Peptide (FP-2)', 'S2 Heptad-repeat-1',
        'S2 Heptad-repeat-2', 'Transmembrane domain', 
        'S2 subunit, intra-virion / Cytoplasmic Tail (CT)', 'Unknown'
    ]
    mutation_counts_df = pd.DataFrame(0, index=sequences_df['Seq Name'], columns=all_regions)
    
    # Calculate the counts for each region
    for _, row in variants_df.iterrows():
        for region in row['Mapped Regions']:
            for seq in row['Sequences']:
                mutation_counts_df.loc[seq, region] += 1
    
    # Merge the mutation counts with the sequences DataFrame
    merged_df = sequences_df.merge(mutation_counts_df, left_on='Seq Name', right_index=True)
    
    # Calculate the correlation matrix
    correlation_df = merged_df[all_regions + ['Stab', 'ACE2', 'nAb']].corr()
    
    # Extract only the correlations between the regions and the fitness components
    correlation_table = correlation_df.loc[all_regions, ['Stab', 'ACE2', 'nAb']]
    
    # Summary statistics in a table-like format
    stats_table = pd.DataFrame({
        'Total Sequences': len(sequences_df),
        'Phi': [f"Mean: {sequences_df['Phi'].mean():.2f}", f"Min: {sequences_df['Phi'].min()}", f"Max: {sequences_df['Phi'].max()}"],
        'Total Variants': [f"Mean: {sequences_df['Total Variants'].mean():.2f}", f"Min: {sequences_df['Total Variants'].min()}", f"Max: {sequences_df['Total Variants'].max()}"],
        'Stab': [f"Mean: {sequences_df['Stab'].mean():.2f}", f"Min: {sequences_df['Stab'].min()}", f"Max: {sequences_df['Stab'].max()}"],
        'ACE2': [f"Mean: {sequences_df['ACE2'].mean():.2f}", f"Min: {sequences_df['ACE2'].min()}", f"Max: {sequences_df['ACE2'].max()}"],
        'nAb': [f"Mean: {sequences_df['nAb'].mean():.2f}", f"Min: {sequences_df['nAb'].min()}", f"Max: {sequences_df['nAb'].max()}"]
    }).T
    
    return stats_table.T, correlation_table





