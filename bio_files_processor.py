import os


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):
    """
    Reads the input fasta-file, in which the sequence can be split into several lines,
    then saves it into a new fasta-file in which each sequence fits into one line

    :param input_fasta: Function convert_multiline_fasta_to_oneline
    :param output_fasta: new fasta-file in which each sequence fits into one line
    """
    if not os.path.isfile(input_fasta):
        raise ValueError('No such fasta file in directory')
    if output_fasta is None:
        output_fasta = 'converted_' + input_fasta
    else:
        output_fasta = output_fasta + '.fasta'
    with open(input_fasta, mode='r') as file, open(output_fasta, mode='w') as out_file:
        out_file.write(file.readline())
        for line in file:
            if '>' in line:
                out_file.write('\n' + line)
            else:
                out_file.write(line.strip())


def bring_indexes(all_genes: list[str], matching_genes: list[str], out_neighboor_genes, n_before: int = 1,
                  n_after: int = 1) -> list[str]:
    """
    The function receives as input a list of all genes found in the file,
    a list of target genes - genes of interest, and then selects genes neighboring the genes of interest,
    located at a given distance before or after the gene of interest

    :param all_genes: list of all genes that was found in .gbk file
    :param matching_genes: list of genes - matches genes of interest
    :param out_neighboor_genes: list to write output - genes thar are neighbours to genes of interest
    :param n_before: number of genes before (>0). Default value is 1
    :param n_after: number of genes after (>0). Default value is 1
    :return: out_neighboor_genes - list to write output - genes thar are neighbours to genes of interest
    """
    for gene in matching_genes:
        gene_index = all_genes.index(gene)
        for index in range(gene_index - n_before if n_before <= gene_index else 0,
                           gene_index + n_after + 1 if n_after + gene_index <= (len(all_genes) - 1) else len(
                               all_genes)):
            if (all_genes[index] != gene and
                        all_genes[index] not in out_neighboor_genes and
                        all_genes[index] not in matching_genes:
                out_neighboor_genes.append(all_genes[index])
    return out_neighboor_genes


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: list[str],
                                   n_before: int = 1, n_after: int = 1,
                                   output_fasta: str = 'output_neighboor_genes_fasta'):
    """
    Select some number of genes before and after each gene of interest from .gbk file and
    save their protein sequence (translation) to a fasta file

    :param input_gbk:  path to the input GBK file
    :param genes: list of genes of interest, next to which neighbors are searched
    :param n_before: number of genes before (>0). Default value is 1
    :param n_after: number of genes after (>0). Default value is 1
    :param output_fasta: name of the output file. The .fasta extension is added to the name
    :return:
    """
    if not os.path.isfile(input_gbk):
        raise ValueError('No such file in directory')
    with open(input_gbk, mode='r') as file, open(output_fasta + '.fasta', mode='w') as out_file:
        dict_gene_seq = {}
        for line in file:
            if '/gene' in line:
                gene_name = line.strip().split('=')[1]
                test = 0
                seq = []
                while test == 0:
                    new_line = file.readline().strip()
                    if '/translation' in new_line:
                        seq.append(new_line.split('="')[1])
                    elif not new_line.startswith('/') and new_line.isupper():
                        if '"' in new_line:
                            seq.append(new_line[:-1])
                            test = 1
                        else:
                            seq.append(new_line)
                dict_gene_seq[gene_name.strip('"')] = ''.join(seq)
        all_genes = list(dict_gene_seq.keys())
        matching_genes = []
        for gene_to_find in genes:
            matching_genes += [s for s in all_genes if gene_to_find in s]

        out_neighboor_genes = bring_indexes(all_genes, matching_genes, [], n_before, n_after)

        for neighboor_gene in out_neighboor_genes:
            out_file.write('>' + neighboor_gene + "\n")
            out_file.write(dict_gene_seq[neighboor_gene] + "\n")


def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta: str):
    """
    Shift the initial position in the one-line fasta with one record

    :param input_fasta: one-line fasta with only one record
    :param shift: an integer (can be negative) - how much to shift the initial position in the file
    :param output_fasta: one-line fasta with one record with changed start position
    """
    with open(input_fasta, mode='r') as file, open(output_fasta, mode='w') as out_file:
        name_line = file.readline()
        seq_line = file.readline().strip()
        if len(seq_line) < abs(shift):
            check_point = len(seq_line) * (abs(shift) // len(seq_line))
            shift = shift - check_point if shift > 0 else shift + check_point
        new_seq_line = seq_line[shift:] + seq_line[:shift]
        out_file.write(name_line)
        out_file.write(new_seq_line)


def parse_blast_output(input_file: str, output_file: str = None):
    """
    Reads the given txt file, for each QUERY select the first line from the Description column.
    Save the set of obtained proteins to a new file in a single column sorted alphabetically

    :param input_file: path to the input txt file
    :param output_file: name of the output file
    """
    if output_file is None:
        output_file = 'parse_blast_blast_results.txt'
    with open(input_file, mode='r') as file, open(output_file, mode='w') as out_file:
        protein_seqs = []
        for line in file:
            if 'Sequences producing significant alignments' in line:
                test = 0
                while test == 0:
                    new_line = file.readline()
                    if 'Description' in new_line:
                        protein_line = file.readline().strip()
                        test = 1
                        if "  " in protein_line:
                            protein_seqs.append(protein_line.split("  ")[0])
        protein_seqs.sort()
        for protein_seq in protein_seqs:
            out_file.write(protein_seq + "\n")
