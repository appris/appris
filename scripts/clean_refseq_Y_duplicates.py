import argparse
import os, re
from datetime import datetime


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--input', '-i',
                    help='input annotation file in GTF format.', 
                    type=str)

args = parser.parse_args()

log_path = os.path.dirname(os.path.abspath(args.input))
with open(os.path.join(log_path,"clean_refseq_Y_duplicates.log"), "w") as log_file:

    # Log file
    log_file.write("# ---------------")
    log_file.write("Running clean_refseq_Y_duplicates.py\n\n")
    script_path = os.path.join(os.getcwd(), __file__)
    log_file.write(" ".join(["python3", script_path, "\\", "\n"]))
    log_file.write(" ".join(["--input", args.input, "\n"]))
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    log_file.write(" ".join(["\nDate:", dt_string]))
    log_file.write(" ".join(["\nUser:", os.getlogin()]))


    duplicates_set = set()
    x_gene_set = set()
    regex = r'gene=([A-Za-z0-9-\.\@_]*);'
    num_dup = 0

    filename, file_extension = os.path.splitext(args.input)
    
    with open(args.input, "r") as in_file, open(args.input.replace(file_extension, '.clean'+file_extension), "w") as out_file:
        for line in in_file:
            if line.startswith("NC_000023.11"):
                gene_list = re.findall(regex, line)

                if len(gene_list) == 0:
                    out_file.write(line)
                elif len(gene_list) == 1:
                    gene = gene_list[0]
                    if gene not in x_gene_set:
                        x_gene_set.add(gene)
                    out_file.write(line)
                else:
                    print(line)
                    break

            elif line.startswith("NC_000024.10"):
                gene_list = re.findall(regex, line)

                if len(gene_list) == 0:
                    out_file.write(line)
                elif len(gene_list) == 1:
                    gene = gene_list[0]
                    if gene not in x_gene_set:
                        out_file.write(line)
                    else:
                        duplicates_set.add(gene)
            else:
                out_file.write(line)

    # Print duplicated genes in log file
    sorted_dup = sorted(duplicates_set)
    log_file.write("\n\n#### Duplicated genes found in chr Y (NC_000024.10):\n")
    for gene in sorted_dup:
        # write each item on a new line
        log_file.write("\t%s\n" % gene)
    log_file.write(" ".join(['#### Number of duplicated genes found:', str(len(sorted_dup))]))
