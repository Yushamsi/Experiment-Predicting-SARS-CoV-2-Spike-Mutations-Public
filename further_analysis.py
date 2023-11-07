
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



