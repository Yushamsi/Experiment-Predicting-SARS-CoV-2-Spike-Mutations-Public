
"PREDICTING SARS-COV-2 SPIKE PROTEIN MUTATIONS BASED ON FITNESS USING GENETIC ALGORITHMS"

---

### Description

This project strives to predict potential mutations in the SARS-CoV-2 spike protein genome using a genetic algorithm (GA). The GA was trained on six amino acid sequences of length 1270 from variants designated as "variants of concern." With the application of a custom fitness function (based on SpikePro) and a custom mutation function (employing a BLOSUM matrix), the GA was run until sequences met a threshold fitness of 2000. This process was repeated to generate ten distinct sequences. These sequences were then analyzed for multiple sequence alignment using PyMSAVis, and relevant information such as fitness, mutations, and mutation count was extracted.

### Background

This study was initiated to delve into the SARS-CoV-2 spike genome, targeting potential advantageous mutations and mutation hotspots. We harnessed an in silico algorithmic approach to predict potential viral genetic mutations, coupled with a tailored fitness function for SARS-CoV-2. Genetic Algorithms (GAs), renowned for mimicking natural selection processes, were used to generate potential viral genetic sequences. The emphasis was on identifying mutations in the spike protein that could enhance viral fitness. While GAs offer a wide array of sequence possibilities, it's paramount to remember that the results reflect the quality of the input data and foundational assumptions.

### Project Objectives

1. Conduct in-depth research on the SARS-CoV-2 spike genome to determine mutation hotspots.
2. Design a GA tailored to predicting viral genetic mutations.
3. Establish a SARS-CoV-2-specific fitness function.
4. Detect spike protein mutations that heighten viral fitness.
5. Engineer a reusable GA to produce multiple sequences for in-depth analysis.
6. Predict the most probable spike protein mutations using GA-generated sequences.
7. Identify specific spike protein mutations and their protein regions contributing to high-fitness variants.

### Repository Structure

- **my_GA_utils.py**: Contains all custom code for data extraction, processing, the GA, custom mutation and fitness functions. Key functions are `run_all_ga` for running the GA and `run_all_results` for extracting data from the generated sequences.
- **SpikePro Files**: Essential dependencies for the custom fitness function.
- **Clustal Omega Dependencies**: Required dependencies for multiple sequence alignment visualization. [Refer to Clustal Omega's official documentation](http://www.clustal.org/omega/).
- **Amino Acid FASTA Sequences**: The initial population sequences fed into the GA.
- **further_analysis.py**: Contains functions for deeper analysis and statistics on the generated sequences.
- **run.ipynb**: Jupyter notebook with the framework to customize and run the GA using different parameters and initial populations.

### Usage

**my_GA_utils.run_all_ga** Parameters:
- folder_path
- num_generations
- num_parents_mating
- fitness_func
- mutation_type
- parent_selection_type
- K_tournament
- crossover_type
- keep_parents
- keep_elitism
- save_solutions
- stop_criteria
- random_seed
- suppress_warnings
- num_runs
- save_and_display_images

**my_GA_utils.run_all_results** Parameters:
- folder_path
- fasta_file_path

**further_analysis.further_analyze_sequences** Parameters:
- file_path_sequences
- file_path_variants

### Dependencies

- The core of the fitness function relies on the [SpikePro algorithm](https://github.com/3BioCompBio/SpikeProSARS-CoV-2).
- For handling multiple sequence alignment, this repository depends on Clustal Omega. [Please refer to the official documentation for installation and usage details](http://www.clustal.org/omega/).

### Acknowledgments

The SpikePro algorithm files were sourced from the [original repository](https://github.com/3BioCompBio/SpikeProSARS-CoV-2).

---

