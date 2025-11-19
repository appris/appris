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
    parser.add_argument("-d", "--trifid", help="Trifid predicctions file.")
    parser.add_argument("--loglevel", default="info", help="Log level")
    parser.add_argument("--logfile", default="info", help="Log file")
    args = parser.parse_args()

    if os.path.dirname(args.logfile) != "":
        create_dir(os.path.dirname(args.logfile))

    logger.add(args.logfile, level=args.loglevel, colorize=False, backtrace=True, diagnose=True)
    logger.info(f"APPRIS_MANUAL has started and its output will be ready here")

    manual_df = pd.read_csv(args.manual_file, header = 0, sep ='\t')
    trifid_df = pd.read_csv(args.trifid, compression='gzip', header=0, sep='\t', quotechar='"')

    manual_df.groupby('Gene_ID').apply(lambda gene_df: fix_gene(gene_df, args.input_folder, trifid_df))


def create_dir(dirpath: str) -> str:
    """mkdir -p python equivalent

    Arguments:
        dirpath {str} -- Path to create the new folder

    Returns:
        absolute_path {str} -- Absolute path of the new folder
    """
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)



def fix_gene(gene_df, input_folder, trifid_df):

    if len(gene_df[['Species', 'Chromosome']].drop_duplicates()) != 1:
        logger.error(f'Diferent Species/Chromosome for the same Gene ID {gene}')
    
    gene = gene_df.iloc[0]['Gene_ID']
    gene_HGNC = gene_df.iloc[0]['Gene_Name']
    chrom = gene_df.iloc[0]['Chromosome']

    for file in ['appris', 'appris.label']:

        logger.info(f"INFO. Checking gene {gene} (chr {chrom}) for file {file}...")
        appris_file = f"{input_folder}/{chrom}/{gene}/{file}"
        appris_df = pd.read_csv(appris_file, sep = '\t', names = appris_header, dtype=str).replace(np.nan, "", regex=True)
        appris_df['Protein length'] = appris_df['Protein length'].replace("", np.nan, regex=True).astype('Int64')
        
        ## CHECK LEGTHS
        gene_length_df = gene_df.sort_values(by = 'Transcript_ID')[['Transcript_ID', 'Length']].reset_index(drop = True)
        appris_length_df = appris_df[appris_df['Transcript ID'].isin(gene_df['Transcript_ID'])] \
                            .sort_values(by = 'Transcript ID')[['Transcript ID', 'Protein length']].reset_index(drop = True)
        df_diff = pd.concat([gene_length_df,appris_length_df.rename(columns={'Transcript ID':'Transcript_ID', 'Protein length':'Length'})]).drop_duplicates(keep=False)
        
        # If the lengths do not match -> Don't change anything. Check it manually.
        if len(df_diff) != 0:
            for i in df_diff['Transcript_ID'].drop_duplicates(keep='first'):
                logger.error(f"ERROR. In gene {gene} ({gene_HGNC}), in chr {chrom}, transcript {i} has change its length or has been removed. Check please.")

        else:
        ## PRINCIPALS
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
                    logger.warning(f"WARNING. In gene {gene} ({gene_HGNC}), in chr {chrom}, transcript {i} has the same length that the manual PRINCIPAL:M. Check please.")

        ## ALTERNATIVES AND MINORS

            # Change previous annotations to 'ALTERNATIVE:M' and 'MINOR:M'
            gene_df_ALT = gene_df[gene_df['Tag'].str.contains('ALTERNATIVE',na=False)]
            gene_df_MINOR = gene_df[gene_df['Tag'].str.contains('MINOR',na=False)]

            # If there are any ALTERNATIVE/MINOR to update:
            if ((len(gene_df_ALT) != 0) or (len(gene_df_MINOR) != 0)):
                appris_df.loc[appris_df['Transcript ID'].isin(gene_df_ALT['Transcript_ID']), 'APPRIS Annotation'] = 'ALTERNATIVE:M'
                appris_df.loc[appris_df['Transcript ID'].isin(gene_df_MINOR['Transcript_ID']), 'APPRIS Annotation'] = 'MINOR:M'

                # If we override the only PRINCIPAL(s) we took as PRINCIPAL the ALTERNATIVE with more 
                # 'norm_trifid_score' between the ALTERNATIVES that has not been annotated manually
                if not appris_df['APPRIS Annotation'].str.contains('PRINCIPAL', na=False).any():
                    ALT_df = appris_df[appris_df['APPRIS Annotation'].str.contains('ALTERNATIVE:[1,2]', regex=True, na=False)]
                    merge_df = ALT_df.merge(trifid_df, how='inner', left_on='Transcript ID', right_on='transcript_id')
                    max_df = merge_df[merge_df['norm_trifid_score'].values == merge_df['norm_trifid_score'].values.max()]

                    # If we found one or more ALTERNATIVE:[1,2] candidates; we took the one(s) 
                    # with the more common length and then save the new file.
                    if len(max_df) != 0:
                        max_df['Count'] = max_df.groupby('Protein length')['Transcript ID'].transform('count')
                        max_df = max_df.sort_values('Count', ascending=False).reset_index()
                        transcripts_list = max_df[max_df['Protein length'] == max_df.iloc[0]['Protein length']]['Transcript ID']

                        for i in transcripts_list:
                            logger.warning(f"WARNING. Gene {gene} ({gene_HGNC}), in chr {chrom}, ended with no PRINCIPAL. Transcript {i} has been taken as PRINCIPAL:M instead. Check please.")
                            appris_df.loc[appris_df['Transcript ID'] == i, 'APPRIS Annotation'] = 'PRINCIPAL:M'

                        appris_df.to_csv(f"{appris_file}", sep='\t', header=False, index=False)
                        logger.info(f"INFO. Gene {gene} (chr {chrom}) has been updated with APPRIS Manual.")

                    # If we are not able to found any ALTERNATIVE:[1,2] candidates to PRINCIPAL:M we
                    # show a warning ERROR in the log but we dont modify the gene
                    else:
                        logger.error(f"ERROR. Gene {gene} ({gene_HGNC}), in chr {chrom}, ended with no PRINCIPAL and could not find any candidate. Changes not saved. Check please.")
                
                # If after update ALTERNATIVE:M and MINOR:M we still have any PRINCIPAL, we just 
                # save the new manually annotated file.
                else:
                    appris_df.to_csv(f"{appris_file}", sep='\t', header=False, index=False)
                    logger.info(f"INFO. Gene {gene} ({gene_HGNC}), in chr {chrom}, has been updated with APPRIS Manual.")
            
            # If there are not any ALTERNATIVE/MINOR to update we just 
            # save the new manually annotated file.
            else:
                appris_df.to_csv(f"{appris_file}", sep='\t', header=False, index=False)
                logger.info(f"INFO. Gene  ({gene_HGNC}), in chr {chrom}, has been updated with APPRIS Manual.")

            # Finally we remove the extra tabs that have been added at the end of the file
            with open(f"{appris_file}", "r") as content_file:
                content = content_file.read()

            content_aux=re.sub('\t+\n','\n',content)

            with open(f"{appris_file}", "w") as output_file:
                output_file.write(content_aux)


if __name__ == "__main__":
    main()