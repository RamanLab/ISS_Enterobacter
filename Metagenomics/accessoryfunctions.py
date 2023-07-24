# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 17:06:42 2023

@author: Pratyay
"""

# Import required packages
import pandas as pd
from ete3 import NCBITaxa, Tree
import numpy as np

# Segment 1: Data preprocessing for every location (e.g. ISS/ JPL/ KSC)
# Required functions
def read_tsv_file(file_path):
    """
    This function will load the dataframe using given path
    
    Parameters
    ----------
    file_path : str
        path to the bracken file

    Returns
    -------
    dataframe : pandas dataframe
        loaded bracken dataframe from the given path

    """
    try:
        # Read the TSV file using Pandas
        dataframe = pd.read_csv(file_path, sep='\t')
        return dataframe
    except FileNotFoundError:
        print("Error: File not found.")
    except Exception as e:
        print("An error occurred:", str(e))
        
        
def segregate_dataframe(dataframe):
    """
    This function will segregate the dataframe into two other dataframes
    One is with the absolute abundance and one with relative abundance

    Parameters
    ----------
    dataframe : pandas dataframe
        loaded bracken dataframe from the given path

    Returns
    -------
    dataframe_num : pandas dataframe
        bracken dataframe with absolute abundance
    dataframe_frac : pandas dataframe
        bracken dataframe with relative abundance

    """
    
    columns_ending_with_num = [col for col in dataframe.columns 
                               if col.endswith('num')]
    columns_ending_with_frac = [col for col in dataframe.columns 
                                if col.endswith('frac')]
    common_columns = list(set(dataframe.columns) - set(columns_ending_with_num)
                          - set(columns_ending_with_frac))

    dataframe_num = dataframe[common_columns + columns_ending_with_num]
    dataframe_frac = dataframe[common_columns + columns_ending_with_frac]

    return dataframe_num, dataframe_frac


def fetch_taxonomy_details(dataframe):
    """
    This function will fetch the taxonomy for the given dataframe with column
    'taxonomy_id' as specified from bracken using ete3 and select only rows
    with 'Bacteria' and 'Fungi'

    Parameters
    ----------
    dataframe : pandas dataframe
        loaded bracken dataframe from the given path

    Returns
    -------
    filtered_dataframe : pandas dataframe
        Updated with 'taxonomy_details' and filtered for only 'Bacteria' 
        and 'Fungi'

    """
    # Initialize NCBI taxonomy database
    ncbi = NCBITaxa()

    # Fetch taxonomy details for each taxonomy ID
    taxonomy_details = []
    for tax_id in dataframe['taxonomy_id']:
        try:
            lineage = ncbi.get_lineage(tax_id)
            names = ncbi.get_taxid_translator(lineage)
            taxonomy = [names[taxid] for taxid in lineage]
            taxonomy_details.append(taxonomy)
        except Exception as e:
            print(f"Error fetching taxonomy for tax_id: {tax_id}. Error: {str(e)}")

    # Add taxonomy details as a new column in the DataFrame
    dataframe['taxonomy_details'] = taxonomy_details

    # Keep rows where taxonomy list contains 'bacteria' or 'fungi'
    filtered_dataframe = dataframe[dataframe['taxonomy_details'].apply(
        lambda x: 'Bacteria' in x)]
        #or 'Fungi' in x)]

    return filtered_dataframe


def remove_columns(dataframe):
    """
    This function removes columns which are not necessary hereafter like
    taxonomy_id, taxonomy_lvl, and our newly generated 'taxononmy_details'

    Parameters
    ----------
    dataframe : pandas dataframe
       dataframe updated with 'taxonomy_details'

    Returns
    -------
    filtered_dataframe : pandas dataframe
        dataframe with mentioned columns deleted

    """
    columns_to_remove = [col for col in dataframe.columns 
                         if col == 'taxonomy_id' 
                         or col == 'taxonomy_lvl'
                         or col == 'taxonomy_details']
    updated_dataframe = dataframe.drop(columns=columns_to_remove)
    
    return updated_dataframe


def convert_name_to_index(dataframe):
    """
    This function will assign the name column as index column

    Parameters
    ----------
    dataframe : pandas dataframe
        NA

    Returns
    -------
    dataframe : pandas dataframe
        name column as index

    """
    dataframe.set_index('name', inplace=True)
    
    return dataframe


def filter_rows_by_threshold(dataframe, threshold=0.01, percentage=0.1):
    """
    This function will filtered out rows aka the taxa based on the given 
    cutoffs

    Parameters
    ----------
    dataframe : pandas dataframe
        name column as index
    threshold : float, optional
        relative abundance cutoff. The default is 0.01.
    percentage : float, optional
        prevalence cutoff. The default is 0.1.

    Returns
    -------
    filtered_dataframe : pandas dataframe
        filtered dataframe

    """
    # Calculate the minimum number of columns that should meet the threshold
    min_columns = int(len(dataframe.columns) * percentage)

    # Filter rows based on the threshold
    filtered_dataframe = dataframe[dataframe.ge(threshold).sum(axis=1) 
                                   >= min_columns]

    return filtered_dataframe

def convert_to_relative_abundance(dataframe):
    """
    This function will convert the absoulte abundance to relative abundance
    
    Parameters
    ----------
    dataframe : pandas dataframe
        absoulte abundance data

    Returns
    -------
    new_dataframe : pandas dataframe
        relative abundance data (dividing every value with column sum)

    """
    column_sums = dataframe.sum()
    new_dataframe = dataframe.div(column_sums)
    return new_dataframe

def adding_taxa_back(dataframe, dataframe_taxainfo):
    """
    This function readd the taxonomy details from the main taxa dataframe

    Parameters
    ----------
    dataframe : pandas dataframe
        relative abundance data
    dataframe_taxainfo : pandas dataframe
        dataframe with taxonomic details

    Returns
    -------
    dataframe : pandas dataframe
        relative abundance data with taxa details

    """
    typedata, taxaid = [], []
    for items in dataframe.index:
        taxa = dataframe_taxainfo[dataframe_taxainfo['name'] 
                    == items]['taxonomy_details'].to_string(index=False)
        
        tax_id = dataframe_taxainfo[dataframe_taxainfo['name'] 
                    == items]['taxonomy_id'].to_string(index=False)
        
        taxaid.append(tax_id)
        
        if 'Bacteria' in taxa:
            typedata.append('Bacteria')
        elif 'Fungi' in taxa:
            typedata.append('Fungi')
        else:
            typedata.append('')
    
    dataframe['taxa'] = typedata
    dataframe['taxa_id'] = taxaid
    
    return dataframe

def add_lineage_to_dataframe(dataframe, tax_id_column):
    """
    Adding all lineage information to the dataframe

    Parameters
    ----------
    dataframe : pandas dataframe
       taxa table
    tax_id_column : string
       column having tax ids

    Returns
    -------
    dataframe : pandas dataframe
       taxa table with lineage information

    """
    ncbi = NCBITaxa()
    lineage_info = []
    for tax_id in dataframe[tax_id_column]:
        lineage = ncbi.get_lineage(tax_id)
        names = ncbi.get_taxid_translator(lineage)
        lineage_names = [names[taxid] for taxid in lineage]
        lineage_str = ', '.join(lineage_names)
        lineage_info.append(lineage_str)
    dataframe['lineage'] = lineage_info
    return dataframe

def process_dataframe(dataframe):
    """
    This function identifies co-existing microbes for given organism, 
    here E. bugandensis

    Parameters
    ----------
    dataframe : pandas dataframe
        taxa table with lineage information

    Returns
    -------
    dataframe : pandas dataframe
        identifies samples which have non-zero values for the selected organism
    nonzero_indices : dictionary
        co-existing microbes in the samples

    """
    # Drop columns 'taxa', 'taxa_id', 'lineage'
    dataframe = dataframe.drop(columns=['taxa', 'taxa_id', 'lineage'])
       
    # Drop columns with zero values where 'Enterobacter bugandensis' is zero
    dataframe = dataframe.loc[:, (dataframe.loc['Enterobacter bugandensis'] != 0)]
    
    # Create a dictionary with non-zero indices for each column
    nonzero_indices = {}
    for column in dataframe.columns:        
        indices = []
        for i, value in dataframe[column].iteritems():
            if value > 0.01:
                indices.append(i)
        nonzero_indices[column] = indices
    
    return dataframe, nonzero_indices

def group_by_family_and_sum(dataframe):
    """
    This function is used for the prepare the table for the family level
    relative abundance

    Parameters
    ----------
    dataframe : pandas dataframe
        dataframe with family information

    Returns
    -------
    grouped_df : pandas dataframe
        dataframe concatenated (sum) across family level

    """
    # Set the 'name' column as the index
    dataframe.set_index('name', inplace=True)

    # Convert non-numeric columns to numeric (if necessary)
    numeric_columns = dataframe.select_dtypes(include=[int, float]).columns
    dataframe[numeric_columns] = dataframe[numeric_columns].apply(pd.to_numeric,
                                                                  errors='coerce')

    # Group by 'Family' column and calculate the sum of numerical columns
    grouped_df = dataframe.groupby('Family').sum()

    return grouped_df
