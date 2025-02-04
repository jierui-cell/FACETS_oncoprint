import pandas as pd 

def select_oncogenic_only(df):
    """
    Filter the dataframe to only include oncogenic mutations.

    Parameters:
    - df: DataFrame containing mutation data.

    Returns:
    - DataFrame containing only oncogenic mutations.
    """
    return df[df['ONCOGENIC'].isin(['Likely Oncogenic', 'Oncogenic'])]

def merge_mutations_fusion_cna(df_mutations, df_fusion, df_cna, samples, oncogenic_only): 
    # Filter mutation data for the specified samples
    df_mutations_subset = df_mutations[df_mutations['Tumor_Sample_Barcode'].isin(samples)]
    df_mutations_subset = df_mutations_subset.rename(columns={'Tumor_Sample_Barcode': 'SAMPLE_ID'})

    # Filter fusion data for the specified samples
    df_fusion_subset = df_fusion[df_fusion['Sample_ID'].isin(samples)]
    df_fusion_subset = df_fusion_subset.rename(columns={'Sample_ID': 'SAMPLE_ID'})

    # Filter CNA data for the specified samples
    df_cna_subset = df_cna[df_cna['SAMPLE_ID'].isin(samples)]

    # filter for oncogenic alterations if required 
    if oncogenic_only: 
        df_mutations_subset = select_oncogenic_only(df_mutations_subset)
        df_fusion_subset = select_oncogenic_only(df_fusion_subset)
        df_cna_subset = select_oncogenic_only(df_cna_subset)


    # Merge mutation, fusion, and CNA data
    df_fusion_cna = pd.merge(df_fusion_subset, df_cna_subset, on='SAMPLE_ID', how='outer') 
    df_fusion_cna = df_fusion_cna.fillna(0).infer_objects(copy=False)
    df_merged = pd.merge(df_mutations_subset, df_fusion_cna, on='SAMPLE_ID', how='outer')
    df_merged = df_merged.fillna(0).infer_objects(copy=False)

    return df_merged
def generate_complete_alteration_table(df_mutations, df_fusion, df_cna, samples, oncogenic_only):
    # alteration includes mutation, fusion, and CNA 
    df_fusion_cna_mutations = merge_mutations_fusion_cna(df_mutations, df_fusion, df_cna, samples, oncogenic_only=oncogenic_only) 
    genes = list(set(list(df_fusion['Site1_Hugo_Symbol'].unique()) + list(df_fusion['Site2_Hugo_Symbol'].unique()) + list(df_mutations['Hugo_Symbol'].unique()) + list(df_cna['HUGO_SYMBOL'].unique())))

    # the columns for the gene table are now a bit tricky 
    # we need to have the gene name and the type of alteration from the following list:
    # ['Mutation', 'Fusion', 'CNA']
    # eg: RIT1 will have 3 columns: RIT1_Mutation, RIT1_Fusion, RIT1_CNA   
    columns = [f"{gene}_{alteration}" for gene in genes for alteration in ['Mutation', 'Fusion', 'CNA']]
    df_gene_table = pd.DataFrame(0, index=samples, columns=columns)

    # fill in the table
    for i in range(df_fusion_cna_mutations.shape[0]):
        row = df_fusion_cna_mutations.iloc[i]
        if row['Site1_Hugo_Symbol'] != 0:
            df_gene_table.loc[row['SAMPLE_ID'], f"{row['Site1_Hugo_Symbol']}_Fusion"] = 1
        if row['Site2_Hugo_Symbol'] != 0:
            df_gene_table.loc[row['SAMPLE_ID'], f"{row['Site2_Hugo_Symbol']}_Fusion"] = 1
        if row['Hugo_Symbol'] != 0:
            df_gene_table.loc[row['SAMPLE_ID'], f"{row['Hugo_Symbol']}_Mutation"] = 1
        if row['HUGO_SYMBOL'] != 0:
            df_gene_table.loc[row['SAMPLE_ID'], f"{row['HUGO_SYMBOL']}_CNA"] = 1

    return df_gene_table

# how to avoid the table to be too spare? 
def generate_mutation_table(df_mutations, samples, oncogenic_only, genes): 
    df_mutations_subset = df_mutations[df_mutations['Tumor_Sample_Barcode'].isin(samples)]
    df_mutations_subset = df_mutations_subset.rename(columns={'Tumor_Sample_Barcode': 'SAMPLE_ID'})
    
    if oncogenic_only: 
        df_mutations_subset = select_oncogenic_only(df_mutations_subset)
    
    columns = [f"{gene}_Mutation" for gene in genes]
    
    # Initialize the DataFrame with empty strings to avoid dtype issues
    df_gene_table = pd.DataFrame('', index=samples, columns=columns, dtype='object')
    
    # Fill in the mutation data
    for i in range(df_mutations_subset.shape[0]):
        row = df_mutations_subset.iloc[i]
        if row['Hugo_Symbol'] not in genes:
            continue
        if df_gene_table.loc[row['SAMPLE_ID'], f"{row['Hugo_Symbol']}_Mutation"] == '':
            df_gene_table.loc[row['SAMPLE_ID'], f"{row['Hugo_Symbol']}_Mutation"] = row['Variant_Classification']
        else:
            df_gene_table.loc[row['SAMPLE_ID'], f"{row['Hugo_Symbol']}_Mutation"] += ";" + row['Variant_Classification']
    
    return df_gene_table


def generate_fusion_table(df_fusion, samples, oncogenic_only, genes): 
    df_fusion_subset = df_fusion[df_fusion['Sample_ID'].isin(samples)]
    df_fusion_subset = df_fusion_subset.rename(columns={'Sample_ID': 'SAMPLE_ID'})
    
    if oncogenic_only: 
        df_fusion_subset = select_oncogenic_only(df_fusion_subset)
    
    columns = [f"{gene}_Fusion" for gene in genes]
    
    # Initialize the DataFrame with empty strings to avoid dtype issues
    df_gene_table = pd.DataFrame('', index=samples, columns=columns, dtype='object')
    
    # Fill in the fusion data
    for i in range(df_fusion_subset.shape[0]):
        row = df_fusion_subset.iloc[i]
        df_gene_table.loc[row['SAMPLE_ID'], f"{row['Site1_Hugo_Symbol']}_Fusion"] = 'Structural Variant'
        if row['Site2_Hugo_Symbol'] in genes:
            df_gene_table.loc[row['SAMPLE_ID'], f"{row['Site2_Hugo_Symbol']}_Fusion"] = 'Structural Variant'

    # if want to go in more details about the fusion type
    # for i in range(df_fusion_subset.shape[0]):
    #     row = df_fusion_subset.iloc[i]
    #     if pd.isna(row['Class']): 
    #         continue
    #     if df_gene_table.loc[row['SAMPLE_ID'], f"{row['Site1_Hugo_Symbol']}_Fusion"] == '':
    #         df_gene_table.loc[row['SAMPLE_ID'], f"{row['Site1_Hugo_Symbol']}_Fusion"] = row['Class']
    #     else:
    #         df_gene_table.loc[row['SAMPLE_ID'], f"{row['Site1_Hugo_Symbol']}_Fusion"] += ";" + row['Class']

    return df_gene_table


def generate_cna_table(df_cna, samples, oncogenic_only, genes): 
    df_cna_subset = df_cna[df_cna['SAMPLE_ID'].isin(samples)]
    
    if oncogenic_only: 
        df_cna_subset = select_oncogenic_only(df_cna_subset)
    
    columns = [f"{gene}_CNA" for gene in genes]
    
    # Initialize the DataFrame with empty strings to avoid dtype issues
    df_gene_table = pd.DataFrame('', index=samples, columns=columns, dtype='object')
    
    # Fill in the CNA data
    for i in range(df_cna_subset.shape[0]):
        row = df_cna_subset.iloc[i]
        if row['HUGO_SYMBOL'] not in genes:
            continue
        if df_gene_table.loc[row['SAMPLE_ID'], f"{row['HUGO_SYMBOL']}_CNA"] == '':
            df_gene_table.loc[row['SAMPLE_ID'], f"{row['HUGO_SYMBOL']}_CNA"] = row['ALTERATION']
        else:
            df_gene_table.loc[row['SAMPLE_ID'], f"{row['HUGO_SYMBOL']}_CNA"] += ";" + row['ALTERATION']
    
    return df_gene_table


def generate_alteration_table_with_details_oncoprint(df_mutation_table, df_fusion_table, df_cna_table, target_genes): 
    # merge the 3 tables together by a way that for each gene we only have one column
    # and within this column we have the type of alteration combined
    # for example: TP53 will have one column with the following content: 'Missense_Mutation;Structural Variant;Amplification'
    # columns are target_genes 
    columns = target_genes
    df_gene_table = pd.DataFrame('', index=df_mutation_table.index, columns=columns, dtype='object')

    for gene in target_genes:
        if f"{gene}_Mutation" in df_mutation_table.columns: 
            df_gene_table[gene] = df_mutation_table[f"{gene}_Mutation"]
        if f"{gene}_Fusion" in df_fusion_table.columns: 
            df_gene_table[gene] += ";" + df_fusion_table[f"{gene}_Fusion"]
        if f"{gene}_CNA" in df_cna_table.columns: 
            df_gene_table[gene] += ";" + df_cna_table[f"{gene}_CNA"]

    # remove the extra ";" 
    # if the string equals to ';;', replace it with ''
    df_gene_table = df_gene_table.map(lambda x: '' if x == ';;' else x)
    # if there are two consecutive ';' in the string, replace them with one ';' 
    df_gene_table = df_gene_table.map(lambda x: x.replace(';;', ';')) 
    # if there is a ';' at the end of the string, remove it
    df_gene_table = df_gene_table.map(lambda x: x[:-1] if (len(x) > 0) and (x[-1] == ';') else x)
    # if there is a ';' at the beginning of the string, remove it
    df_gene_table = df_gene_table.map(lambda x: x[1:] if (len(x) > 0) and (x[0] == ';') else x) 
                
    return df_gene_table

# only includes True False for each value instead of detailed alteration types
def generate_alteration_table_simplified(df_mutation_table, df_fusion_table, df_cna_table, target_genes): 
    df_gene_table = generate_alteration_table_with_details_oncoprint(df_mutation_table, df_fusion_table, df_cna_table, target_genes)
    df_gene_table = df_gene_table.map(lambda x: True if x != '' else False)
    return df_gene_table

def transform_alterations_table(df_alteration_table):
    """
    Transforms a wide table of alterations into a long-format table
    with columns: [sample, track, event].
    
    - df_alteration_table index or first column should be the sample ID.
    - Each column beyond that is a gene (track).
    - Cell values can be empty or semicolon-separated strings of alterations.
    
    Returns a DataFrame with columns:
        ["sample", "track", "event"] 
    in long format, where each event is a separate row.
    """

    # Ensure the sample IDs are a column named 'sample'
    # (If your sample IDs are currently the DataFrame's index, we reset and rename)
    if df_alteration_table.index.name is None:
        # If the index has no name but truly stores sample IDs,
        # you might do: df_alteration_table.index.name = 'sample'
        df_alteration_table.index.name = 'sample'
    df_long = df_alteration_table.reset_index().melt(
        id_vars='sample',        # keep sample as an identifier
        var_name='track',        # gene name (e.g., ERBB2, TP53, etc.)
        value_name='event'       # alteration(s)
    )

    # Drop rows where event is NaN or an empty string
    df_long = df_long.dropna(subset=['event'])
    df_long = df_long[df_long['event'].str.strip() != '']

    # Split semicolon-separated events into multiple rows
    # (e.g. "Nonsense_Mutation;Frame_Shift_Del" => 2 rows)
    rows = []
    for _, row in df_long.iterrows():
        # Split on semicolons
        alterations = [e.strip() for e in row['event'].split(';') if e.strip()]
        for alt in alterations:
            rows.append({
                'sample': row['sample'],
                'track':  row['track'],
                'event':  alt
            })

    # Create the final DataFrame
    df_result = pd.DataFrame(rows, columns=['sample', 'track', 'event'])
    
    return df_result
