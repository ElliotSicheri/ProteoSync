"""
Helps identify unconserved segments in a protein sequence by:
    -   BLAST searching the sequence against several small species-specific databases
    -   Creating a clustalw sequence alignment with the top hits from each search
    -   Finding the closest ortholog to the sequence from the PBD database and aligning its secondary structure
    -   Finding the AlphaFold model for the sequence and aligning its secondary structure
"""

from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast.Applications import NcbiblastpCommandline, NcbimakeblastdbCommandline
from Bio.PDB import PDBList, DSSP, PDBParser
from datetime import datetime, date
from os.path import exists
from pathlib import Path
import os
import shutil
import requests

res_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T',
            'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A',
            'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'MSE': 'M'}

dssp_dict = {'H': 'A', 'B': 'B', 'E': 'B', 'I': '-', 'S': '-', 'G': '-', 'T': '-', 'C': '-', '-': '-'}


def blast_search(seq_file: str, threshold: int, path_list: list[str]) -> (dict[str: str], dict[str: int]):
    """BLAST searches the sequence from the given file against a series of species-specific databases and writes the
    top hits from each in an output file if they have a % identity above the given threshold.

    Does not search for any species found in the given exclude_list. Species in include_list have their top hit written
    to the output file, whether it has % identity above the given threshold or not.

    Returns a dict of species that returned a hit above the threshold mapped to the accession code of the hit sequence,
    and a list of species that did not return a hit above the threshold.

    Parameters:
        - seq_file: str, the name of the file that contains the query sequence
        - threshold: int, % identity threshold for which BLAST hits are included
        - path_list: list[str], list of paths which contain a database to be checked
    """

    formatted_results = ''  # Properly formatted result string
    successful_searches = {}  # Species that returned a hit above the threshold
    failed_searches = {}  # Species that did not return a hit above the threshold
    for path in path_list:
        blast = NcbiblastpCommandline(cmd='blastp', query=seq_file, db=path + '/' + os.path.basename(path) + '.fasta',
                                      out='Desktop/ProteoSync/temp_files/raw_data.txt', outfmt=0,
                                      num_threads=8, num_alignments=5, num_descriptions=5)
        blast()
        sequence, code, score = _return_first_hit()
        species = os.path.basename(path)
        if score > threshold:
            formatted_results += '> ' + species + '_[' + str(score) + '] \n'
            formatted_results += sequence + '\n'
            successful_searches[species] = code
        else:
            failed_searches[species] = score

    with open('Desktop/ProteoSync/temp_files/formatted_results.txt', 'w') as file:
        # Clears the contents of formatted_results.txt and writes the formatted data to the file
        file.truncate(0)
        file.write(formatted_results)

    return successful_searches, failed_searches


def _return_first_hit() -> (str, str, int):
    """Returns the sequence, accession code and % identity of the first hit in the BLAST search output file."""

    hit = ''
    code = ''
    score = 0
    started = False
    segment_dict = {}

    with open('Desktop/ProteoSync/temp_files/raw_data.txt', 'r') as file:
        result_lines = file.readlines()
        file.close()

    for line in result_lines:
        if line != '' and line[0] == '>':
            # If the line contains the name of a BLAST hit after the sequence has been added, stop adding sequence lines
            # to the output string
            if started:
                break
            else:
                c_start = line.find(' ')
                c_start = line.find(' ', c_start + 1)
                code = line[2:c_start]
                started = True
        elif 'Identities' in line:
            i = line.find('%')
            score_str = line[i - 3:i].replace('(', '')
            score = int(score_str)
        elif line != '' and line[0:5] == 'Sbjct':
            # If the line contains the sequence of the hit, add it to the output string
            start_value = int(line[7:11])
            end = line.find(' ', 12)
            if end == 12:
                end = line.find(' ', 13)
            seq_line = line[12:end].lstrip()
            segment_dict[start_value] = seq_line
        elif 'No hits found' in line:
            return '', '', 0

    keys = list(segment_dict.keys())
    while len(keys) > 0:
        smallest = min(keys)
        hit += segment_dict[smallest]
        keys.remove(smallest)

    return hit.replace('-', ''), code, score


def get_alignment(input_file: str) -> str:
    """Creates a sequence alignment from the sequences in the input file, then writes it to an output file. Returns the
    name of the output file."""

    cmd = ClustalwCommandline('clustalw', infile=input_file, seqnos='ON')
    stdout, stderr = cmd()
    # print(stdout)
    # Extracts the name of the output file from the readout
    s = stdout.rfind('[')
    e = stdout.rfind(']')
    output_file = stdout[s + 1:e]
    # Returns the aligned output file
    return output_file


def structure_search(seq_file: str, rec_count: int = 0, search_range: (int, int) = (0, 0)) -> list[(str, str)]:
    """Runs a BLAST search of the PDB database for sequences that are closest to the given sequence. Returns the PBD
    codes and a string representing the secondary structure of these proteins."""

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
    blast = NcbiblastpCommandline(cmd='blastp', query='Desktop/ProteoSync/temp_files/sub_query.txt',
                                  db='Desktop/ProteoSync/databases/pdb/pdb.fasta',
                                  out='Desktop/ProteoSync/temp_files/raw_data_pdb.txt', outfmt=0, num_threads=8,
                                  num_alignments=2, num_descriptions=2)
    blast()

    # Pulls the PDB codes of the hits
    with open('Desktop/ProteoSync/temp_files/raw_data_pdb.txt') as file:
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
            if score >= 25:
                pdb_codes.append([code, chain])
        elif 'No hits found' in line:
            return []

    # Generates strings representing the secondary structure of each hit
    structure_strs = []
    for pdb_code in pdb_codes:
        code, chain = pdb_code
        # Downloads the pdb file, if it doesn't exist already
        if not exists('Desktop/ProteoSync/downloads/pdb_structures/pdb' + code + '.pdb'):
            pdb_list = PDBList(server='http://files.pdbj.org')
            filename = pdb_list.retrieve_pdb_file(pdb_code=code, pdir='Desktop/ProteoSync/downloads/pdb_structures',
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
        filename = 'Desktop/ProteoSync/downloads/pdb_structures/pdb' + code + '.pdb'

        # Locates unmodelled areas in the PDB file and extracts the protein sequence. If there are any sections added to
        # the front of the sequence with position labelled as <= 0, count the number of residues in said section, so it
        # can be accounted for when adding unmodelled sections to the structure string.
        with open(filename, 'r') as pdb_file:
            lines = pdb_file.readlines()
            pdb_file.close()

        unmodelled_sections = []  # List of unmodelled sections
        row_tally = -1
        remark_num = 0
        seq_str = ''  # Protein sequence
        chain_start = 0  # The residue number the chain starts at
        lowest_index = 1  # lowest residue index, including those before chain_start
        lowest_diff = 0

        for line in lines:
            if 'MISSING RESIDUES' in line:
                # Finds the start of the missing residues section, starts collection of unmodelled section data
                row_tally = 5
                remark_num = int(line[6:10])
            elif row_tally > 0:
                # Counts down 5 rows from the beginning of section
                row_tally -= 1
            elif row_tally == 0:
                if int(line[6:10]) == remark_num and line[19] == chain:
                    res_num = int(line[21:26])
                    if len(unmodelled_sections) != 0 and unmodelled_sections[-1][1] == res_num - 1:
                        unmodelled_sections[-1][1] = res_num
                    else:
                        unmodelled_sections.append([res_num, res_num])
                elif int(line[6:10]) != remark_num:
                    row_tally = -1
            elif (line[0:6] == 'DBREF ' or line[0:6] == 'DBREF1') and line[12] == chain:
                chain_start = int(line[14:18])
            elif line[0:6] == 'SEQADV':
                # Finds residues marked as position <= 0
                index_str = line[17:22].strip().replace('-', '')
                if index_str.isnumeric():
                    index = int(line[17:22])
                    if index - chain_start < 0:
                        if index - chain_start < lowest_diff:
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
                            seq_str += 'X'
                        i += 4

        unmodelled_residues = set()
        for section in unmodelled_sections:
            unmodelled_residues = unmodelled_residues.union(set(range(section[0], section[1] + 1)))

        unmodelled_seq = seq_str
        #removes unmodelled sections from the sequence
        for i in unmodelled_residues:
            unmodelled_seq = unmodelled_seq[0:i-lowest_index] + 'X' + unmodelled_seq[1+i-lowest_index:]
        unmodelled_seq = unmodelled_seq.replace('X', '')

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
        for x in range(0, len(dssp)):
            a_key = list(dssp.keys())[x]
            dssp_structure_str += dssp[a_key][2]

        # Converts the string from DSSP format
        structure_str = ''
        s_i = 0
        while s_i < len(dssp_structure_str):
            structure_str += dssp_dict[dssp_structure_str[s_i]]
            s_i += 1

        # Aligns query sequence with the sequence of the PDB hit, aligns structure string to query sequence
        with open('Desktop/ProteoSync/temp_files/s_aln.txt', 'w') as s_aln:
            s_aln.write('> Query_sequence\n')
            s_aln.write(query + '\n')
            s_aln.write('> Query_seq_2\n')
            s_aln.write(query + '\n')
            s_aln.write('> Hit_sequence\n')
            s_aln.write(unmodelled_seq + '\n')
            s_aln.close()

        cmd = ClustalwCommandline('clustalw', infile='Desktop/ProteoSync/temp_files/s_aln.txt')
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

    if rec_count == 0:
        # Reruns any large unmodelled sections on their own
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
    """Downloads an alphafold prediction listed by the given UniProt id. Returns a string representing its secondary
    structure aligned to the query string."""

    with open(seq_file, 'r') as query_file:
        query = query_file.read()
        query_file.close()

    filename = "Desktop/ProteoSync/downloads/AF_structures/AF-" + uniprot_id + ".pdb"

    if not exists('Desktop/ProteoSync/downloads/AF_structures/AF-' + uniprot_id + '.pdb'):
        # Downloads the alphafold file
        link_pattern = 'http://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v3.pdb'
        url = link_pattern.format(uniprot_id)
        response = requests.get(url) #verify = false

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
    with open('Desktop/ProteoSync/temp_files/s_aln.txt', 'w') as s_aln:
        s_aln.write('> Query_sequence\n')
        s_aln.write(query + '\n')
        s_aln.write('> Hit_sequence\n')
        s_aln.write(seq_str + '\n')
        s_aln.close()

    cmd = ClustalwCommandline('clustalw', infile='Desktop/ProteoSync/temp_files/s_aln.txt')
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


def update_database() -> None:
    """Updates the local PDB database."""

    # Downloads the sequence file
    print('Downloading sequence file...')
    pdb_list = PDBList(server='http://files.pdbj.org')
    pdb_list.get_seqres_file('pdb.fasta')

    # Writes today's date to the last_update file
    datefile = open('Desktop/ProteoSync/last_update.txt', 'w')
    datefile.truncate(0)
    today_date = date.today().strftime('%y %m %d')
    datefile.write('20' + today_date)
    datefile.close()

    # Removes old database files
    print('Deleting old database files...')
    os.remove('Desktop/ProteoSync/databases/pdb/pdb.fasta.phr')
    os.remove('Desktop/ProteoSync/databases/pdb/pdb.fasta.pin')
    os.remove('Desktop/ProteoSync/databases/pdb/pdb.fasta.psq')

    # Moves the file to the correct directory
    original = 'pdb.fasta'
    final = 'Desktop/ProteoSync/databases/pdb/pdb.fasta'
    shutil.move(original, final)

    # Formats the sequence file into a BLAST database
    print('Formatting database...')
    db_line = NcbimakeblastdbCommandline(dbtype='prot', input_file='Desktop/ProteoSync/databases/pdb/pdb.fasta',
                                         title='PDB_Database')
    db_line()

    os.remove('Desktop/ProteoSync/databases/pdb/pdb.fasta')


def assemble_output_file(alignment_file: str, threshold: int, successes: dict[str:str], fails: dict[str:int],
                         pdb_list: list[(str, str)] = [], alpha_struc_str: str = '',
                         filename: str = '') -> str:
    """Given the results of other functions from this package, parses a full alignment file. Returns the file name."""

    pdb_indexes = []
    for i in range(0,len(pdb_list)):
        pdb_indexes.append(0)

    if filename == '':
        # Gets current date and time to name output file
        current_date = date.today().strftime("%b%d_%Y")
        now = datetime.now()
        current_time = now.strftime("%H:%M")
        output_name = 'Alignment_' + current_date + '_' + current_time + '.txt'
        pymol_output_name = 'Alignment_' + current_date + '_' + current_time + '_pymol.txt'
        output = open('Desktop/ProteoSync/output/alignment_output/' + output_name, 'w')
        pymol_output = open('Desktop/ProteoSync/output/pymol_scripts/' + pymol_output_name, 'w')
    else:
        if exists('Desktop/ProteoSync/output/alignment_output/' + filename + '.txt'):
            ext = 1
            while exists('Desktop/ProteoSync/output/alignment_output/' + filename + '(' + str(ext) + ').txt'):
                ext += 1
            filename = filename + '(' + str(ext) + ')'

        output_name = filename + '.txt'
        pymol_output_name = filename + '_pymol.txt'
        output = open('Desktop/ProteoSync/output/alignment_output/' + output_name, 'w')
        pymol_output = open('Desktop/ProteoSync/output/pymol_scripts/' + pymol_output_name, 'w')

    with open(alignment_file, 'r') as file:
        alignment = file.readlines()
        file.close()

    # Identifies the beginning and end of the alignment section
    line = alignment[3].strip()
    if line[len(line) - 1] != '-':
        seq_end = alignment[3].rfind(' ')
        seq_start = alignment[3].rfind(' ', 0, seq_end) + 1
    else:
        seq_start = alignment[3].rfind(' ') + 1
        seq_end = len(alignment[3]) - 1

    identity_list = []
    two_dot_list = []
    one_dot_list = []
    no_dot_list = []

    # Writes results in output file
    output.write('% identity threshold: ' + str(threshold) + '\n \n')
    cons_index = 0
    struc_index = 0
    a_struc_index = 0
    highest_score = 0
    comp_seq = ''
    for line in alignment[3:]:
        if line != '\n' and line[0] != ' ':
            s = line.find('[') + 1
            score_str = line[s:s + 3].replace(']', '')
            score = int(score_str)
            if score >= highest_score:
                highest_score = score
                comp_seq = line[seq_start:seq_end]
        elif line[0] == ' ':
            conservation = line[seq_start:seq_end]
            for i in range(0, len(comp_seq)):
                if comp_seq[i] == ' ':
                    break
                elif comp_seq[i] != '-':
                    match conservation[i]:
                        case '*':
                            identity_list.append(cons_index + 1)
                        case ':':
                            two_dot_list.append(cons_index + 1)
                        case '.':
                            one_dot_list.append(cons_index + 1)
                        case ' ':
                            no_dot_list.append(cons_index + 1)
                    cons_index += 1

            if pdb_list != []:
                for i in range (0, len(pdb_list)):
                    struct = pdb_list[i]
                    output.write('STRUCTURE [' + struct[0].upper() + ']:' + ' ' * (seq_start - 17))
                    for char in comp_seq:
                        if char == '-':
                            output.write(' ')
                        elif pdb_indexes[i] < len(struct[1]):
                            output.write(struct[1][pdb_indexes[i]])
                            pdb_indexes[i] += 1
                    output.write('\n')

            if alpha_struc_str != '':
                output.write('ALPHAFOLD STRUCTURE:' + ' ' * (seq_start - 20))
                for char in comp_seq:
                    if char == '-':
                        output.write(' ')
                    elif a_struc_index < len(alpha_struc_str):
                        output.write(alpha_struc_str[a_struc_index])
                        a_struc_index += 1
                output.write('\n')
        output.write(line)
    output.write('\n')

    if len(successes) > 0:
        output.write('Accession codes:\n')
        for success in successes:
            output.write('\t- ' + success.replace('_', ' ') + ': ' + successes[success] + '\n')
    output.write('\n')

    if len(fails) > 0:
        output.write('\nThe following species lack a sufficiently similar protein sequence: \n')
        for species in fails.keys():
            output.write('\t- ' + species.replace('_', ' ') + ' (' + str(fails[species]) + '%)\n')
    output.write('\n')

    if len(identity_list) > 0:
        pymol_output.write('color density, resi ' + str(identity_list[0]))
        if len(identity_list) > 1:
            for i in identity_list[1:]:
                pymol_output.write('+' + str(i))
        pymol_output.write('\n')

    if len(two_dot_list) > 0:
        pymol_output.write('color marine, resi ' + str(two_dot_list[0]))
        if len(two_dot_list) > 1:
            for i in two_dot_list[1:]:
                pymol_output.write('+' + str(i))
        pymol_output.write('\n')

    if len(one_dot_list) > 0:
        pymol_output.write('color lightblue, resi ' + str(one_dot_list[0]))
        if len(one_dot_list) > 1:
            for i in one_dot_list[1:]:
                pymol_output.write('+' + str(i))
        pymol_output.write('\n')

    if len(no_dot_list) > 0:
        pymol_output.write('color white, resi ' + str(one_dot_list[0]))
        if len(no_dot_list) > 1:
            for i in no_dot_list[1:]:
                pymol_output.write('+' + str(i))
        pymol_output.write('\n')

    return output_name

