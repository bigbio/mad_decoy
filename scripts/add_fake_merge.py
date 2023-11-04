import logging

import pandas as pd
import glob


def process_and_merge_csvs(folder_path):
    all_files = glob.glob(folder_path + "/*.tsv")
    df_list = []

    for filename in all_files:
        # Step 1: Read each CSV file
        df = pd.read_csv(filename)

        ## remove all proteins with q-value Nan (this is a bug in the code)
        df = df[df['protein_global_qvalue'].notna()]
        ## continue if df is empty
        if df.empty:
            continue

        print(f"Processing file: {filename}")

        # Step 2: Filter out rows with is_decoy set to 1
        filtered_df = df[df['is_decoy'] == 0]

        # Step 3: Sort the dataframe by protein_global_qvalue for the current file
        sorted_df = filtered_df.sort_values(by='protein_global_qvalue')

        # Step 4: Process for each unique q-value in the current file
        fake_entries = []
        unique_qvalues = sorted_df['protein_global_qvalue'].unique()

        fake_protein = 0

        for qvalue in unique_qvalues:
            # Todo: Lukas - I think this two conditions are always true (sorted_df['is_decoy'] == 0) and
            # len(sorted_df[(sorted_df['protein_global_qvalue'] <= qvalue) & (sorted_df['is_decoy'] == 1)]) is always 0
            # because we filtered out rows with is_decoy set to 1 in line 24. Please double check.

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


def update_qvalues(df, filter_qvalue : float, filter_decoy :bool):

    df = df.sort_values(by='protein_global_qvalue')

    # Count the total number of decoy proteins
    targets = (1-df["is_decoy"]).cumsum()
    decoys = df["is_decoy"].cumsum()

    # Calculate the ratio
    
    # Update the protein_global_qvalue for each protein
    df['protein_adjusted_qvalue'] = decoys / targets

    df['protein_adjusted_qvalue'] = df['protein_adjusted_qvalue'][::-1].cummin()[::-1]
    if filter_qvalue is not None:
        df = df[df['protein_adjusted_qvalue'] <= float(filter_qvalue)]
    if filter_decoy:
        df = df[df['is_decoy'] == 1]

    return df


def remove_duplicates(df):
    # Sort dataframe by 'protein_accessions' and 'protein_global_qvalue'
    # to ensure we keep the lowest q-value when dropping duplicates
    df = df.sort_values(by=['protein_accessions', 'protein_global_qvalue'])
    
    # Drop duplicate entries based on 'protein_accessions' and 'condition'
    # Since the dataframe is sorted by q-value in ascending order, the first occurrence (lowest q-value) is retained
    df = df.drop_duplicates(subset=['protein_accessions'], keep='first')
    
    return df

# add the main script to start the tool from command line
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process and merge CSV files')
    parser.add_argument('--folder_path', dest='folder_path', required=True, help='Path to folder containing CSV files')
    parser.add_argument('--output_file', dest='output_file', required=True, help='Path to output CSV file')
    parser.add_argument('--filter_qvalue', dest='filter_qvalue', required=False, help='Filter q-value')
    parser.add_argument('--filter_decoy', dest='filter_decoy', action='store_true', required=False, help='Filter decoy')
    args = parser.parse_args()

    output_df = process_and_merge_csvs(args.folder_path)
    output_df = remove_duplicates(output_df)
    output_df = update_qvalues(output_df, args.filter_qvalue, args.filter_decoy)
    output_df.to_csv(args.output_file, index=False)
