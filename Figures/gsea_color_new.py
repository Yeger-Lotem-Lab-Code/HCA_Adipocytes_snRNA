import pandas as pd
import os

#folder_path = r'/Users/mayaziv/Dropbox (Personal)/HCA_paper/tables/gsea_sub/'
b = "redo_adipo10_com"
folder_path = r'/Users/mayaziv/OneDrive - post.bgu.ac.il/Adipose/redo_GSEA/adipo10_com/'
df = pd.read_csv(r'/Users/mayaziv/OneDrive - post.bgu.ac.il/Adipose/Summarized_objects_5_3/gsea_pathways_1.csv')

# Get a list of all files in the folder
file_list = os.listdir(folder_path)
#
# num_pathways = {'Angiogenesis': 18, 'ECM': 9, 'Adaptive immune': 38, 'Innate immune': 26, 'Lipid Metabolism': 14,
#                 'Ribosomal': 5}
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
            r'/Users/mayaziv/OneDrive - post.bgu.ac.il/Adipose/redo_GSEA/adipo10_com/' + b + a + '.csv')
for i in range(len(results_df)):
    if path_num_dict[results_df.at[i, 'Cluster']] == 0:
        results_df.at[i, 'adj_score'] = 0
    else:
        results_df.at[i, 'adj_score'] = results_df.at[i, 'Sum_of_Scores'] / path_num_dict[results_df.at[i, 'Cluster']]

# table = pd.pivot_table(results_df, values='Sum_of_Scores', index=['Cluster'],
#                        columns=['Pathways'])
# table.fillna(0, inplace=True)
# table.to_csv(
#     r'/Users/mayaziv/OneDrive - post.bgu.ac.il/Adipose/Summarized_objects_5_3/check_gsea/' + b + 'score_sum.csv')
table = pd.pivot_table(results_df, values='adj_score', index=['Cluster'],
                       columns=['Pathways'])
table.fillna(0, inplace=True)
table.to_csv(
    r'/Users/mayaziv/OneDrive - post.bgu.ac.il/Adipose/redo_GSEA/' + b + 'adj_score.csv')
table = pd.pivot_table(results_df, values='Mean_score', index=['Cluster'],
                       columns=['Pathways'])
table.fillna(0, inplace=True)

table.to_csv(
    r'/Users/mayaziv/OneDrive - post.bgu.ac.il/Adipose/Summarized_objects_5_3/check_gsea/' + b + 'mean_score.csv')

# Load the Excel file
df1 = pd.read_excel(
    os.path.join(
        r'/Users/mayaziv/Dropbox (Personal)/HCA_paper/tables/SX_tables - for revision/Table S6 SA1-7 and VA1-8 hAd1-7 gseGO results.xlsx'
    ),
    engine='openpyxl', sheet_name=None)
pathways_list= df

flattened_list = [item for sublist in pathways_list for item in sublist]

print(flattened_list)
keys = list(df1.keys())
for key in keys[1:]:
    print(key)
    df2 = df1[key]
    df2.columns = df2.iloc[0]
    df2 = df2[1:]
    # Check if values in 'Name' column are in the target_names list
    df2['Visualization'] = df2['Description'].isin(flattened_list).astype(int)
    df1[key] = df2
df2 = df1.copy()

with pd.ExcelWriter(
        r'/Users/mayaziv/Dropbox (Personal)/HCA_paper/tables/SX_tables/Table S6 SA1-7 and VA1-8 gseGO results_visualization.xlsx') as writer:
    for key in df1.keys():
        print(key)
        df1[key].to_excel(writer, sheet_name=key)
    df2.to_excel(writer, sheet_name='Sheet_name_2')

###new_ visualization
# Flatten the lists in the DataFrame
flattened_list = [item for sublist in df.values.tolist() for item in sublist]

# Convert DataFrame to dictionary
df_dict = {col: df[col].dropna().tolist() for col in df.columns}



keys = list(df1.keys())
for key in keys[1:]:  # Skip the first sheet
    df2 = df1[key].copy()
    df2.columns = df2.iloc[0]
    df2 = df2[1:].reset_index(drop=True)
    # Check if values in 'Description' column are in the flattened_list
    df2.loc[:, 'Visualization'] = df2['Description'].isin(flattened_list).astype(int)
    df1[key] = df2

df2 = df1.copy()

# Save the processed DataFrame to Excel
output_path = '/path/to/output.xlsx'
with pd.ExcelWriter(output_path) as writer:
    for key in df1.keys():
        df1[key].to_excel(writer, sheet_name=key)
    df2.to_excel(writer, sheet_name='Sheet_name_2')