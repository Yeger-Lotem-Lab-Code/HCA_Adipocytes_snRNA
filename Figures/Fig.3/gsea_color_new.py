import pandas as pd
import os


b = "adipo10"
folder_path = r'insert_your_input_folder_path_containting GSEA_enrichment.R containing csv file for each cluster with the enriched pathways with pv.adj and NES'
df = pd.read_csv(r'insert a path to a csv  with column of the gsea pathways related to each enrichment function, see example in Fig3 folder, "Processes selected for visualization.csv')

# Get a list of all files in the folder
file_list = os.listdir(folder_path)

non_empty_columns = {col: df[col].dropna().tolist() for col in df.columns}
path_num_dict = {}
results_df = pd.DataFrame(columns=['Pathways', 'Sum_of_Scores'])
# Iterate through the list of files
for file_name in file_list:
    # Construct the full path to the file
    full_path = os.path.join(folder_path, file_name)
    if file_name != '.DS_Store':
        cluster = file_name.split('_')
        if cluster[0] == 'cluster':
            a = cluster[1]
        else:
            a = cluster[0]
        df = pd.read_csv(full_path)
        df = df[df['p.adjust'] <= 0.05]
        count = 0
        save_df = pd.DataFrame()
        for key, values in non_empty_columns.items():
            print(key, values)
            filtered_df = df[df['Description'].isin(values)]
            save_df = save_df.append(filtered_df)
            sum_of_scores = filtered_df['NES'].sum()
            row_num = len(filtered_df)
            count= count+row_num

            if row_num != 0:
                Mean_score = sum_of_scores / row_num
            else:
                Mean_score = 0.0
            if sum_of_scores == 0:
                sum_of_scores = 0.0

            #print(num_path)

            results_df = results_df.append(
                {'Cluster': a, 'Pathways': key, 'Sum_of_Scores': sum_of_scores, 'Mean_score': Mean_score},
                ignore_index=True)


        path_num_dict[a]=len(save_df)
        save_df.to_csv(
            r'insert_your_output_path' + b + a + '.csv')
for i in range(len(results_df)):
    if path_num_dict[results_df.at[i, 'Cluster']] == 0:
        results_df.at[i, 'adj_score'] = 0
    else:
        results_df.at[i, 'adj_score'] = results_df.at[i, 'Sum_of_Scores'] / path_num_dict[results_df.at[i, 'Cluster']]


table = pd.pivot_table(results_df, values='adj_score', index=['Cluster'],
                       columns=['Pathways'])
table.fillna(0, inplace=True)
table.to_csv(
    r'insert_your_output_path' + b + 'adj_score.csv')

