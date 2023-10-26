import pandas as pd
import glob

def process_and_merge_csvs(folder_path):
    all_files = glob.glob(folder_path + "/*.tsv")
    df_list = []

    for filename in all_files:
        # Step 1: Read each CSV file
        df = pd.read_csv(filename)

        # Step 2: Filter out rows with is_decoy set to 1
        filtered_df = df[df['is_decoy'] == 0]

        # Step 3: Sort the dataframe by protein_global_qvalue for the current file
        sorted_df = filtered_df.sort_values(by='protein_global_qvalue')

        # Step 4: Process for each unique q-value in the current file
        fake_entries = []
        unique_qvalues = sorted_df['protein_global_qvalue'].unique()

        fake_protein = 0

        for qvalue in unique_qvalues:
            count = len(sorted_df[(sorted_df['protein_global_qvalue'] <= qvalue) & (sorted_df['is_decoy'] == 0)])
            num_fake_entries = int(count * qvalue) - len(sorted_df[(sorted_df['protein_global_qvalue'] <= qvalue) & (sorted_df['is_decoy'] == 1)])

            for _ in range(num_fake_entries):
                fake_protein += 1
                fake_entry = {
                    'protein_accessions': 'fake_protein_' + str(fake_protein),
                    'protein_global_qvalue': qvalue,
                    'is_decoy': 1,
                    'condition': 'fake_condition',  # or whatever placeholder you'd like
                    'dataset_accession': sorted_df['dataset_accession'].iloc[1]  # or whatever placeholder you'd like
                }
                fake_entries.append(fake_entry)

        fake_df = pd.DataFrame(fake_entries)
        processed_df = pd.concat([sorted_df, fake_df], ignore_index=True).sort_values(by='protein_global_qvalue')
        df_list.append(processed_df)

    # Step 5: Merge all the processed dataframes
    merged_df = pd.concat(df_list, ignore_index=True)

    return merged_df


def update_qvalues(df):
    # Count total number of decoy proteins
    targets = (1-df["is_decoy"]).cumsum()
    decoys = df["is_decoy"].cumsum()

    # Calculate the ratio
    
    # Update the protein_global_qvalue for each protein
    df['protein_global_qvalue'] = decoys / targets

    df['protein_global_qvalue'] = df['protein_global_qvalue'][::-1].cummin()[::-1]

    return df


def remove_duplicates(df):
    # Sort dataframe by 'protein_accessions' and 'protein_global_qvalue'
    # to ensure we keep the lowest q-value when dropping duplicates
    df = df.sort_values(by=['protein_accessions', 'protein_global_qvalue'])
    
    # Drop duplicate entries based on 'protein_accessions' and 'condition'
    # Since the dataframe is sorted by q-value in ascending order, the first occurrence (lowest q-value) is retained
    df = df.drop_duplicates(subset=['protein_accessions'], keep='first')
    
    return df

# After merging and processing:
output_df = process_and_merge_csvs('in')

# Remove duplicate protein entries:
output_df = remove_duplicates(output_df)

# Update the qvalues:
output_df = update_qvalues(output_df)

# Save the final dataframe to a CSV
output_df.to_csv('out/final_processed_and_merged.csv', index=False)
