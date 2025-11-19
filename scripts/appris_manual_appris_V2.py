#!/usr/bin/env python
"""
File name: appris_manual_appris.py
Author: Daniel Cerdán Vélez
Email: dcerdan@cnio.es
Created: 2025-05-25
Version: 1.0

Description:
    This script uses the APPRIS manual annotations table 
    (APPRIS_manual_{date}.tsv) to correct the genes that 
    have been predicted using the APPRIS pipeline.

Args:
    file_loc (str): The file location of the spreadsheet
    print_cols (bool): A flag used to print the columns to the console
        (default is False)

Returns:
    list: a list of strings representing the header columns
"""


# Import statements
import argparse
import re
import os
import sys

import pandas as pd
import numpy as np
from loguru import logger

appris_header = ['Ensembl Gene ID', 'Gene name (HGNC)', 'Transcript ID', 'Translation ID', \
            'Is translated?', 'Transcript type', 'Not found tag', 'CCDS ID', 'Transcript support level', \
            'Protein length', 'Functional residues (firestar)', 'Structure score (Matador)', \
            'Conservation (CORSAIR)','Domain Score (SPADE)', 'Trans-membrane helices (THUMP)', \
            'Signal sequence (CRASH)', 'Trifid Score', 'Peptides', 'APPRIS score', 'APPRIS Annotation']
            
def main():
    parser = argparse.ArgumentParser(
        description="Command-line arguments parser", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-s", "--species", default="GRCh38", help="Species.")
    parser.add_argument("-i", "--input_folder", help="Folder with the annotation that must be modified.")
    parser.add_argument("-f", "--manual_file", help="TSV file with the manual annotations for PRINCIPAL, ALTERNATIVE OR MINOR.")
    parser.add_argument("--loglevel", default="INFO", help="Log level")
    parser.add_argument("--logfile", default="log", help="Log file")
    args = parser.parse_args()

    if os.path.dirname(args.logfile) != "":
        create_dir(os.path.dirname(args.logfile))

    logger.add(args.logfile, level=args.loglevel, colorize=False, backtrace=True, diagnose=True)
    logger.info(f"APPRIS_MANUAL has started and its output will be ready here")

    manual_df = pd.read_csv(args.manual_file, header = 0, sep ='\t')
    # Remove possible whitespaces that could be problematic
    manual_df_trimmed = manual_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    manual_df_trimmed.groupby('Gene_ID').apply(lambda gene_df: fix_gene(gene_df, args.input_folder))


def create_dir(dirpath: str) -> str:
    """mkdir -p python equivalent

    Arguments:
        dirpath {str} -- Path to create the new folder

    Returns:
        absolute_path {str} -- Absolute path of the new folder
    """
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)



def fix_gene(gene_df, input_folder):

    if len(gene_df[['Species', 'Chromosome']].drop_duplicates()) != 1:
        logger.error(f'Diferent Species/Chromosome for the same Gene ID {gene}')
    
    gene = gene_df.iloc[0]['Gene_ID']
    gene_HGNC = gene_df.iloc[0]['Gene_Name']
    chrom = gene_df.iloc[0]['Chromosome']

    for file in ['appris', 'appris.label']:

        # We just going to show logs for 'appris' file, because they will be the same for 'appris.label'
        if file == "appris":
            logger.info(f"INFO. Checking gene {gene} ({gene_HGNC}), in chr {chrom}, for APPRIS Manual.)")

        appris_file = f"{input_folder}/{chrom}/{gene}/{file}"
        appris_df = pd.read_csv(appris_file, sep = '\t', names = appris_header, dtype=str).replace(np.nan, "", regex=True)
        appris_df['Protein length'] = appris_df['Protein length'].replace("", np.nan, regex=True).astype('Int64')
        
        ## CHECK LEGTHS
        gene_length_df = gene_df.sort_values(by = 'Transcript_ID')[['Transcript_ID', 'Length']].reset_index(drop = True)
        appris_length_df = appris_df[appris_df['Transcript ID'].isin(gene_length_df['Transcript_ID'])] \
                            .sort_values(by = 'Transcript ID')[['Transcript ID', 'Protein length']].reset_index(drop = True)

        df_diff = pd.concat([gene_length_df,appris_length_df.rename(columns={'Transcript ID':'Transcript_ID', 'Protein length':'Length'})]).drop_duplicates(keep=False)
        # If the lengths do not match -> Don't change anything. Check it manually.
        if len(df_diff) != 0:
            for i in df_diff['Transcript_ID'].drop_duplicates(keep='first'):
                if file == "appris":
                    logger.error(f"ERROR. In gene {gene} ({gene_HGNC}), in chr {chrom}, transcript {i} has change its length or has been removed. Check please.")

        else:
        ## PRINCIPAL

            # If exists any PRINCIPAL:M to update:
            gene_df_PRINC = gene_df[gene_df['Tag'].str.contains('PRINCIPAL', na=False)]
            if len(gene_df_PRINC) != 0:
                # Set old PRINCIPALS to 'ALTERNATIVE:M'
                appris_df.loc[appris_df['APPRIS Annotation'].str.contains('PRINCIPAL', na=False), 'APPRIS Annotation'] = 'ALTERNATIVE:M'

                # Set manual PRINCIPALS to 'PRINCIPAL:M'
                appris_df.loc[appris_df['Transcript ID'].isin(gene_df_PRINC['Transcript_ID']), 'APPRIS Annotation'] = 'PRINCIPAL:M'

            # Check other transcript with same length that manual PRINCIPALS.
            same_length_df = appris_df.loc[appris_df['Protein length'].isin(gene_df_PRINC['Length'])]
            same_length_df = same_length_df[same_length_df['APPRIS Annotation'].str.contains('PRINCIPAL:M', na=False) == False]
            if len(same_length_df) != 0:
                for i in same_length_df['Transcript ID'].drop_duplicates(keep='first'):
                    if file == "appris":
                        logger.warning(f"WARNING. In gene {gene} ({gene_HGNC}), in chr {chrom}, transcript {i} has the same length that PRINCIPAL:M. Check please.")

        ## ALTERNATIVE
            
            # If new ALTERNATIVE:M is already a PRINCIPAL we do not take it into account
            gene_df_ALT = gene_df[gene_df['Tag'].str.contains('ALTERNATIVE',na=False)]
            APPRIS_PRINC = appris_df[appris_df['APPRIS Annotation'].str.contains('PRINCIPAL', na=False)]['Transcript ID']
            gene_df_ALT_filtered = gene_df_ALT[~gene_df_ALT['Transcript_ID'].isin(APPRIS_PRINC)]

            if len(gene_df_ALT_filtered) != 0:
                appris_df.loc[appris_df['Transcript ID'].isin(gene_df_ALT_filtered['Transcript_ID']), 'APPRIS Annotation'] = 'ALTERNATIVE:M'

        ## MINOR

            # We check that a PRINCIPAL:M exists before change a MINOR:M
            gene_df_MINOR = gene_df[gene_df['Tag'].str.contains('MINOR',na=False)]
            if len(gene_df_MINOR) != 0:
                APPRIS_PRINC_M = appris_df[appris_df['APPRIS Annotation'].str.contains('PRINCIPAL:M', na=False)]
                if len(APPRIS_PRINC_M) != 0:
                    appris_df.loc[appris_df['Transcript ID'].isin(gene_df_MINOR['Transcript_ID']), 'APPRIS Annotation'] = 'MINOR:M'
                else:
                    if file == "appris":
                        logger.error(f"ERROR. Gene {gene} ({gene_HGNC}), in chr {chrom}, has a MINOR:M but there isn't a PRINCIPAL:M. Changes not saved. Check please.")
                

            appris_df.to_csv(f"{appris_file}", sep='\t', header=False, index=False)

            # Finally we remove the extra tabs that have been added at the end of the file
            with open(f"{appris_file}", "r") as content_file:
                content = content_file.read()

            content_aux=re.sub('\t+\n','\n',content)

            with open(f"{appris_file}", "w") as output_file:
                output_file.write(content_aux)
        

if __name__ == "__main__":
    main()