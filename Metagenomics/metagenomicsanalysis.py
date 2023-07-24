# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 16:59:06 2023

@author: Pratyay Sengupta

This code is prepared for the metagenomics analysis

ISS_bracken_species: Metagenomics (sheet: PMA bracken)
Metagenomics_Eb_presence.tsv: Metagenomics (sheet: Filtered table)
    (after manual intervention of family level data, obtained from modified_db)
"""

# Import required packages
import pandas as pd
from ete3 import NCBITaxa
from accessoryfunctions import *

if __name__ =="__main__":
    # Section 1 
    # Set file path
    file_path = r"S:\Microbiome Comparison\Data\1. Raw data\ISS_bracken_species.tsv"
    
    # Reading the dataframe
    df = read_tsv_file(file_path)
    
    # Segregating reads from the main bracken file
    df_num, df_frac = segregate_dataframe(df)
    
    # Fetching taxonomy details and minor modifications of the table
    taxadf = fetch_taxonomy_details(df_num)
    df_num = remove_columns(taxadf)
    df_num = convert_name_to_index(df_num)
    
    # Converting absolute abundance to relative abundance
    df_rel = convert_to_relative_abundance(df_num)
    
    # Filtering data based on abundance and prevalence
    df_filtered = filter_rows_by_threshold(df_rel, 
                                           threshold=0.01, percentage=0.05)
    
    # Readding the taxa details to the filtered dataframe
    df_filtered = adding_taxa_back(df_filtered, taxadf)
    df_with_phylum = add_lineage_to_dataframe(df_filtered, 'taxa_id')
    
    # Modify table based on the organism of interest
    modified_df, nonzero_indices = process_dataframe(df_with_phylum)
    
    
    # Section 2: 
    # Load the Ebugandensis presence file
    Eb_file_path = r"S:\IITM-JPL_Enterobacter\metagenomics\Metagenomics_Eb_presence.tsv"
    
    # Reading the dataframe
    Eb_df = read_tsv_file(Eb_file_path)
    
    Eb_family_table = group_by_family_and_sum(Eb_df)
