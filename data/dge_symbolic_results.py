import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from mygene import MyGeneInfo



# Define the CORRECTED DGE analysis function
def run_pydeseq2_corrected(counts, metadata, control_group, treatment_group):
    """
    Runs pyDESeq2 analysis comparing a treatment group to a control group.
    The contrast is correctly set to Treatment vs. Control.
    """
    filter_metadata = metadata[(metadata['Condition'] == control_group ) | (metadata['Condition'] == treatment_group)]
    filtered_counts = counts[counts.index.isin(filter_metadata.index)]

    dds = DeseqDataSet(
        counts=filtered_counts,
        metadata=filter_metadata,
        design_factors="Condition",
        refit_cooks=True,
    )

    dds.deseq2()

    # Corrected contrast for conventional log2fc sign
    deseq_stats = DeseqStats(dds, contrast=["Condition", treatment_group, control_group])
    deseq_stats.summary()
    results_df = deseq_stats.results_df
    
    # Add a column for the comparison being made
    treatment_name = [k for k, v in conditions_dict.items() if v == treatment_group][0]
    results_df['treatment_protocol'] = treatment_name
    
    results_df = results_df.reset_index() # 'ensgene' becomes a column
    return results_df



# Preprocess the counts data
counts = pd.read_csv('airway_rawcounts_more_exps.csv')
counts = counts.set_index('ensgene')
row_sums = counts.sum(axis=1)
counts = counts[row_sums > 0]
counts = counts.T
counts = counts.loc[:, ~counts.columns.duplicated()]

# Prepare metadata
metadata = pd.read_csv('airway_metadata_more_exps.csv', index_col='srr_id')
metadata.index.name = None
metadata = metadata[['treatment_protocol']].rename(columns={'treatment_protocol': 'Condition'})
conditions_dict = {'Untreated': 'C', 'Albuterol': 'A', 'Albuterol_Dexamethasone': 'AD', 'Dexamethasone': 'D'}
metadata['Condition'] = metadata['Condition'].map(conditions_dict)
metadata = metadata.dropna(subset=['Condition'])
counts = counts[counts.index.isin(metadata.index)]

# Gene ID Conversion
print("Fetching gene symbols from Ensembl IDs using mygene.info...")
# Get a unique list of all Ensembl gene IDs from our counts data
ensembl_ids = counts.columns.tolist()

# Query the mygene API
mg = MyGeneInfo()
results = mg.querymany(
    ensembl_ids,
    scopes='ensembl.gene',
    fields='symbol',
    species='human',
    as_dataframe=True,
    df_index=True
)

# Create the mapping dictionary from Ensembl ID -> Gene Symbol
mapping_dict = results['symbol'].to_dict()
print(f"Successfully mapped {len(mapping_dict)} out of {len(ensembl_ids)} genes.")

# Run all comparisons and build final DataFrame
all_results_df = pd.DataFrame()
control_code = 'C'

for treatment_name, treatment_code in conditions_dict.items():
  if treatment_code == control_code:
    continue
  
  print(f"\nRunning comparison for Untreated vs. {treatment_name}...")
  # Run the DGE analysis with the corrected function
  results = run_pydeseq2_corrected(counts, metadata, control_code, treatment_code)
  
  # Append the results to our final dataframe
  all_results_df = pd.concat([all_results_df, results], ignore_index=True)

print("\nAll DGE analyses complete.")

# Add Symbols and Clean the Final DataFrame
print("Mapping Ensembl IDs to gene symbols in the final results...")
# Use the mapping_dict created in Step 3 to add the 'symbol' column
all_results_df['symbol'] = all_results_df['ensgene'].map(mapping_dict)
mg.http_client.close() # Close API connection to avoid error

# Data Cleaning (mirroring the steps from Week 3)
# Remove genes that could not be mapped to a symbol
print(f"Original number of rows: {len(all_results_df)}")
cleaned_df = all_results_df.dropna(subset=['symbol'])
print(f"Rows after dropping unmapped genes: {len(cleaned_df)}")

# Remove duplicate gene symbols within each comparison group, keeping the first
cleaned_df = cleaned_df.drop_duplicates(subset=['symbol', 'treatment_protocol'], keep='first')
print(f"Rows after dropping duplicate symbols per comparison: {len(cleaned_df)}")

# Remove Ensembl gene ID column
sym_counts_df = cleaned_df.drop(columns=['ensgene'])

# Save the final, corrected, and symbolic results to a CSV file
output_filename = 'complete_DGE_comparisons_symbolic.csv'
sym_counts_df.to_csv(output_filename, index=False)

print(f"\nSuccessfully saved the final symbolic DGE results to '{output_filename}'")