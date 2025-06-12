from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.PDB import PDBList, DSSP, PDBParser
from os.path import exists
from pathlib import Path
import requests

base_path = 'Desktop/ProteoSync'

res_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T',
            'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A',
            'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'MSE': 'M'}

dssp_dict = {'H': 'A', 'B': 'B', 'E': 'B', 'I': '-', 'S': '-', 'G': '-', 'T': '-', 'C': '-', '-': '-'}


def structure_search(seq_file: str, rec_count: int = 0, search_range: (int, int) = (0, 0)) -> list[(str, str)]:
    """Runs a BLAST search of the local BLAST-formatted PDB sequence database for sequences that are closest to the
    given sequence. Downloads the corresponding .pdb files from the PDB, and extracts the secondary structure for the
    hit sequences. Aligns the secondary structure to the query sequence and returns as a string.

    Recursively searches large regions left unmodelled by the found structures, to account for structures that model
    single domains of a larger protein.

    Parameters:
        -   seq_file: str, the path to the file that contains the query sequence
        -   rec_count: int, current depth of recursion
        -   search_range: (int, int), range of the sequence to search. (0, 0) defaults to searching the entire sequence.

    Returns:
        List of tuples each containing the PBD code of the hit sequence and a string representing the secondary
        structure of that chain aligned to the query sequence.
    """


    with open(seq_file, 'r') as query_file:
        query = query_file.read()
        query_file.close()

    # Clause preventing infinite recursion
    if rec_count >= 2:
        return []

    with open('Desktop/ProteoSync/temp_files/sub_query.txt', 'w') as subfile:
        if not search_range == (0, 0):
            sub_query = query[search_range[0] - 1:search_range[1] - 1]
            subfile.write(sub_query)
        else:
            subfile.write(query)
        subfile.close()

    # Searches the PDB database
    blast = NcbiblastpCommandline(cmd='blastp', query=base_path+'/temp_files/sub_query.txt',
                                  db=base_path+'/databases/pdb/pdb.fasta',
                                  out=base_path+'/temp_files/raw_data_pdb.txt', outfmt=0, num_threads=8,
                                  num_alignments=2, num_descriptions=2)
    blast()

    # Pulls the PDB codes of the hits
    with open(base_path+'/temp_files/raw_data_pdb.txt') as file:
        result_lines = file.readlines()
        file.close()

    code = ''
    chain = ''
    pdb_codes = []
    for line in result_lines:
        if line != '' and line[0] == '>':
            # If the line contains the name of a BLAST hit, pull the pdb code from it
            code = line[1:6].strip()
            chain = line[6:8].strip().replace('_', '')
        elif line.find('Identities') != -1:
            i = line.find('%')
            score_str = line[i - 3:i].replace('(', '')
            score = int(score_str)
            if score >= 45:
                pdb_codes.append([code, chain])

    if len(pdb_codes) == 0:
        if rec_count == 0:
            print("No sufficiently similar sequences were found in the PDB.")
        return []

    # Generates strings representing the secondary structure of each hit
    structure_strs = []
    for pdb_code in pdb_codes:
        code, chain = pdb_code
        # Downloads the pdb file, if it doesn't exist already
        if not exists(base_path+'/downloads/pdb_structures/pdb' + code + '.pdb'):
            pdb_list = PDBList(server='http://files.pdbj.org')
            filename = pdb_list.retrieve_pdb_file(pdb_code=code, pdir=base_path+'/downloads/pdb_structures',
                                                  file_format='pdb')
            # Some large structures are not available in .pdb format. These structures cause an error when trying to
            # move the file, upon which the current structure will be skipped.
            try:
                path = Path(filename)
                path.rename(path.with_suffix('.pdb'))
            except (Exception, ):
                print('An error occured while downloading this structure. It is possible that it is not available in '
                      '.pdb format. Skipping this structure...')
                continue
        filename = base_path+'/downloads/pdb_structures/pdb' + code + '.pdb'

        # Locates unmodelled areas in the PDB file and extracts the protein sequence. If there are any sections added to
        # the front of the sequence with position labelled as <= 0, count the number of residues in said section, so it
        # can be accounted for when removing unmodelled sections from the structure string.
        with open(filename, 'r') as pdb_file:
            lines = pdb_file.readlines()
            pdb_file.close()

        unmodelled_sections = []  # List of unmodelled sections
        row_tally = -1
        remark_num = 0
        seq_str = ''  # Protein sequence
        chain_start = 0  # The residue number the chain starts at
        lowest_index = 0  # lowest residue index, including those before chain_start
        lowest_diff = 0

        for line in lines:
            if 'MISSING RESIDUES' in line:
                # Finds the start of the missing residues section, starts collection of unmodelled section data
                row_tally = 5
                remark_num = int(line[6:10])
            elif row_tally > 0:
                # Counts down 5 rows from the beginning of unmodelled residue section
                row_tally -= 1
            elif row_tally == 0:
                if int(line[6:10]) == remark_num and line[19] == chain:
                    res_num = int(line[21:26])
                    if len(unmodelled_sections) != 0 and unmodelled_sections[-1][1] == res_num - 1:
                        unmodelled_sections[-1][1] = res_num
                    else:
                        unmodelled_sections.append([res_num, res_num])
                elif int(line[6:10]) != remark_num:
                    row_tally = -1  # Reached end of unmodelled residue section
            elif (line[0:6] == 'DBREF ' or line[0:6] == 'DBREF1') and line[12] == chain:
                chain_start = int(line[14:18])
                lowest_index = chain_start
            elif line[0:6] == 'SEQADV':
                # Finds residues marked as positions < chain_start
                index_str = line[17:22].strip().replace('-', '')
                if index_str.isnumeric():
                    index = int(line[17:22])
                    if line[16] == chain and index - chain_start < lowest_diff:
                        lowest_diff = index - chain_start
                        lowest_index = index
            elif line[0:6] == 'SEQRES' in line:
                # Pulls protein sequence
                if line[11] == chain:
                    seq_line = line[19:].rstrip()
                    i = 0
                    while i < len(seq_line):
                        res = seq_line[i:i + 3]
                        if res in res_dict.keys():
                            seq_str += res_dict[res]
                        else:
                            seq_str += 'X'  # Accounts for uncommon residue modifications
                        i += 4

        if lowest_index < chain_start:
            start_offset = lowest_index
        else:
            start_offset = chain_start

        # Remove unmodelled sections from the sequence
        if len(unmodelled_sections) > 0:
            unmodelled_seq = ''
            unmodelled_seq += seq_str[0:unmodelled_sections[0][0]-start_offset]
            for i in range(len(unmodelled_sections) - 1):
                start_i, end_i = unmodelled_sections[i][1], unmodelled_sections[i+1][0]
                unmodelled_seq += seq_str[1+start_i-start_offset:end_i-start_offset]
            unmodelled_seq += seq_str[1+unmodelled_sections[-1][1]-start_offset:]
        else:
            unmodelled_seq = seq_str

        # Parses the PDB file to get a structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_code, filename)
        model = structure[0]

        # Runs DSSP on the structure model to get the secondary structure.
        # Some structures cannot be properly processed by DSSP, such as structures with their residues stripped. This
        # causes an error, upon which the current structure is skipped.
        try:
            dssp = DSSP(model, filename)
        except (Exception, ):
            print('An error occurred while processing this structure\'s pdb file. Skipping this structure...\n')
            continue

        # Converts the secondary structure to a string
        dssp_structure_str = ''
        for key in list(dssp.keys()):
            if key[0] == chain:
                dssp_structure_str += dssp[key][2]

        # Converts the string from DSSP format
        structure_str = ''
        for ch in dssp_structure_str:
            structure_str += dssp_dict[ch]

        # Aligns query sequence with the sequence of the PDB hit, aligns structure string to query sequence
        with open(base_path+'/temp_files/s_aln.txt', 'w') as s_aln:
            s_aln.write('> Query_sequence\n')
            s_aln.write(query + '\n')
            s_aln.write('> Query_seq_2\n')
            s_aln.write(query + '\n')
            s_aln.write('> Hit_sequence\n')
            s_aln.write(unmodelled_seq + '\n')
            s_aln.close()

        cmd = ClustalwCommandline('clustalw', infile=base_path+'/temp_files/s_aln.txt')
        stdout, stderr = cmd()
        # Extracts the name of the output file from the readout
        s = stdout.rfind('[')
        e = stdout.rfind(']')
        output_file = stdout[s + 1:e]

        with open(output_file, 'r') as output:
            aln_lines = output.readlines()
            output.close()

        query_aln = ''
        hit_aln = ''
        for line in aln_lines:
            if 'Query_sequence' in line:
                query_aln += line[20:].strip()
            elif 'Hit_sequence' in line:
                hit_aln += line[20:].strip()

        # Align the secondary structure string to the query sequence

        i = 0  # Index in aligned sequences
        st_i = 0  # Index in structure_str
        final_struc_str = ''  # Final aligned structure string

        while i < len(query_aln):
            # It is assumed that it is impossible for both aligned sequences to have '-' at the same index
            if query_aln[i] == '-':
                # For sections that are missing from the query sequence, skip that part of the structure string
                st_i += 1
            elif hit_aln[i] == '-':
                # For sections missing from the hit sequence, list section as unmodelled.
                final_struc_str += 'X'
            else:
                # Otherwise, plot the secondary structure as normal
                if st_i < len(structure_str):
                    final_struc_str += structure_str[st_i]
                    st_i += 1
                else:
                    final_struc_str += 'X'
            i += 1

        # Creates a new unmodelled sections list that is aligned to the structure string
        unmodelled_sections_new = []
        for i in range(0, len(final_struc_str)):
            if final_struc_str[i] == 'X':
                if len(unmodelled_sections_new) != 0 and unmodelled_sections_new[-1][1] == i:
                    unmodelled_sections_new[-1][1] = i + 1
                else:
                    unmodelled_sections_new.append([i + 1, i + 1])

        structure_strs.append((code, final_struc_str, unmodelled_sections_new))

    code, top_str, top_unmodelled = structure_strs[0]

    return_list = []
    for struc in structure_strs:
        new_code = struc[0]
        is_in = False
        for old_struc in return_list:
            if old_struc[0] == new_code:
                is_in = True
                break
        if not is_in:
            return_list.append(struc)

    # Recursively reruns any large unmodelled sections on their own

    if rec_count == 0:
        rerun_list = []
        for sec in top_unmodelled:
            if sec[1] - sec[0] >= 20:
                rerun_list.append(sec)

        for rerun_sec in rerun_list:
            rerun_list = structure_search(seq_file, rec_count + 1, rerun_sec)
            for new_struc in rerun_list:
                new_code = new_struc[0]
                is_in = False
                for old_struc in return_list:
                    if old_struc[0] == new_code:
                        is_in = True
                        break
                if not is_in:
                    return_list.append(new_struc)

    return return_list


def alpha_struc_search(seq_file: str, uniprot_id: str) -> str:
    """Downloads an alphafold prediction from UniProt listed by the given UniProt id. Returns a string representing its
    secondary structure aligned to the query sequence in seq_file.

    Parameters:
        -   seq_file: str, path to the file containing the query sequence
        -   uniprot_id: uniprot

    Returns:
        A string representing the secondary structure of the alphafold prediction, aligned to the query sequence
    """

    with open(seq_file, 'r') as query_file:
        query = query_file.read()
        query_file.close()

    filename = base_path+'/downloads/AF_structures/AF-' + uniprot_id + ".pdb"

    if not exists(base_path+'/downloads/AF_structures/AF-' + uniprot_id + '.pdb'):
        # Downloads the alphafold file
        link_pattern = 'http://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v3.pdb'
        url = link_pattern.format(uniprot_id)
        response = requests.get(url)

        if not exists(filename):
            print('An error occured while downloading the AlphaFold prediction. '
                  'It is possible that a prediction is not available from UniProt.')
            return ''

        with open(filename, 'w') as pdb_file:
            for line in response.iter_lines():
                pdb_file.write(str(line)[2:len(line) - 2] + '\n')
            pdb_file.close()

    with open(filename, 'r') as pdb_file:
        lines = pdb_file.readlines()
        pdb_file.close()

    seq_str = ''  # Protein sequence
    seq_start = False  # Whether the sequence section of the file has been reached
    for line in lines:
        if 'SEQRES' in line:
            # Pulls protein sequence
            seq_start = True
            seq_line = line[19:].rstrip()
            i = 0
            while i < len(seq_line):
                res = seq_line[i:i + 3]
                seq_str += res_dict[res]
                i += 4
        elif seq_start:
            # If past the sequence section of the file, break the loop.
            break

    # Parses the PDB file to get a structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('alpha', filename)
    model = structure[0]

    # Runs DSSP on the structure model to get the secondary structure
    dssp = DSSP(model, filename)

    # Converts the secondary structure to a string
    dssp_structure_str = ''
    for x in range(0, len(dssp)):
        a_key = list(dssp.keys())[x]
        dssp_structure_str += dssp[a_key][2]

    # Converts the string from DSSP format and adds unmodelled sections
    structure_str = ''
    i = 0
    while i < len(dssp_structure_str):
        structure_str += dssp_dict[dssp_structure_str[i]]
        i += 1

    # Aligns query sequence with the sequence of the PDB hit, aligns structure string to query sequence
    with open(base_path+'/temp_files/s_aln.txt', 'w') as s_aln:
        s_aln.write('> Query_sequence\n')
        s_aln.write(query + '\n')
        s_aln.write('> Hit_sequence\n')
        s_aln.write(seq_str + '\n')
        s_aln.close()

    cmd = ClustalwCommandline('clustalw', infile=base_path+'/temp_files/s_aln.txt')
    stdout, stderr = cmd()
    # Extracts the name of the output file from the readout
    s = stdout.rfind('[')
    e = stdout.rfind(']')
    output_file = stdout[s + 1:e]

    with open(output_file, 'r') as output:
        aln_lines = output.readlines()
        output.close()

    query_aln = ''
    hit_aln = ''
    for line in aln_lines:
        if 'Query_sequence' in line:
            query_aln += line[20:].strip()
        elif 'Hit_sequence' in line:
            hit_aln += line[20:].strip()

    i = 0  # Index in aligned sequences
    st_i = 0  # Index in structure_str
    final_struc_str = ''  # Final aligned structure string

    while i < len(query_aln):
        # It is assumed that it is impossible for both aligned sequences to have '-' at the same index
        if query_aln[i] == '-':
            # For sections that are missing from the query sequence, skip that part of the structure string
            st_i += 1
        elif hit_aln[i] == '-':
            # For sections missing from the hit sequence, assume that the corresponding section of the query
            # sequence is not structurally relevant.
            if st_i < len(structure_str) and structure_str[st_i] == 'X':
                final_struc_str += structure_str[st_i]
                st_i += 1
            else:
                final_struc_str += 'X'
        else:
            # Otherwise, plot the secondary structure as normal
            if st_i < len(structure_str):
                final_struc_str += structure_str[st_i]
                st_i += 1
            else:
                final_struc_str += '-'
        i += 1

    return final_struc_str