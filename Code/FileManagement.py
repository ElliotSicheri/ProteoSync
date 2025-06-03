from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.PDB import PDBList
from datetime import datetime, date
from os.path import exists
import os
import shutil
import requests

base_path = 'Desktop/ProteoSync'
# base_path = '.'

def update_database() -> None:
    """Updates the local PDB database."""

    # Downloads the sequence file
    print('Downloading sequence file...')
    pdb_list = PDBList(server='http://files.pdbj.org')
    pdb_list.get_seqres_file('pdb.fasta')

    # Writes today's date to the last_update file
    datefile = open(base_path+'/databases/pdb_last_update.txt', 'w')
    datefile.truncate(0)
    today_date = date.today().strftime('%y %m %d')
    datefile.write('20' + today_date)
    datefile.close()

    # Removes old database files
    print('Deleting old database files...')
    os.remove(base_path+'/databases/pdb/pdb.fasta.phr')
    os.remove(base_path+'/databases/pdb/pdb.fasta.pin')
    os.remove(base_path+'/databases/pdb/pdb.fasta.psq')

    # Moves the file to the correct directory
    original = 'pdb.fasta'
    final = base_path+'/databases/pdb/pdb.fasta'
    shutil.move(original, final)

    # Formats the sequence file into a BLAST database
    print('Formatting database...')
    db_line = NcbimakeblastdbCommandline(dbtype='prot', input_file=base_path+'/databases/pdb/pdb.fasta',
                                         title='PDB_Database')
    db_line()

    os.remove(base_path+'/databases/pdb/pdb.fasta')


def get_fastas_from_uniprots(uniprots: list[str]) -> list[str]:
    """Takes a list of UniProt codes and returns a list of corresponding protein sequences"""
    seq_list = []
    for code in uniprots:
        filepath = base_path + '/downloads/fastas/' + code + '.fasta'
        if not exists(filepath):
            # If we have not already downloaded the sequence, do so
            url = 'http://www.uniprot.org/uniprot/' + code + '.fasta'
            response = requests.get(url)

            with open(filepath, 'w') as file:
                # for line in response.iter_lines():
                #     file.write(str(line)[2:len(line) + 2] + '\n')
                # file.close()
                lines = list(response.iter_lines())
                if lines:
                    file.write(lines[0].decode('utf-8') + '\n')
                    sequence = ''.join(line.decode('utf-8') for line in lines[1:])
                    file.write(sequence)

        seq = ''
        with open(filepath, 'r') as file:
            next(file)
            for line in file:
                seq += line
            file.close()
        seq_list.append(seq)

    return seq_list


def assemble_output_file(alignment_file: str, threshold: int, len_threshold: int, successes: dict[str:str],
                         fails: dict[str:int], pdb_list: list[(str, str)] = [], exclude_list: list[str] = [],
                         alpha_struc_str: str = '', filename: str = '') -> str:
    """Given the results of other functions from this package, parses a full alignment file. Returns the file name."""

    pdb_indexes = []
    for i in range(0, len(pdb_list)):
        pdb_indexes.append(0)

    if filename == '':
        # Gets current date and time to name output file
        current_date = date.today().strftime("%b%d_%Y")
        now = datetime.now()
        current_time = now.strftime("%H:%M")
        filename = 'Run_' + current_date + '_' + current_time

        if not exists(base_path+'/output/' + filename):
            os.makedirs(base_path+'/output/' + filename)

        output_name = current_date + '_' + current_time + '_Alignment.txt'
        pymol_output_name = current_date + '_' + current_time + '_Pymol.txt'
        output = open(base_path+'/output/' + filename + '/' + output_name, 'w')
        pymol_output = open(base_path+'/output/' + filename + '/' + pymol_output_name, 'w')
    else:
        if exists(base_path+'/output/' + filename):
            ext = 1
            while exists(base_path+'/output/' + filename + '(' + str(ext) + ')'):
                ext += 1
            filename = filename + '(' + str(ext) + ')'

        os.makedirs(base_path+'/output/' + filename)

        output_name = filename + '_Alignment.txt'
        pymol_output_name = filename + '_Pymol.txt'
        output = open(base_path+'/output/' + filename + '/' + output_name, 'w')
        pymol_output = open(base_path+'/output/' + filename + '/' + pymol_output_name, 'w')

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
    output.write('% identity threshold: ' + str(threshold) + '\n')
    output.write('% length variability threshold: ' + str(len_threshold) + '\n\n')

    if len(exclude_list) > 0:
        output.write('The following species/groups were manually excluded from the search: \n')
        for s in exclude_list:
            output.write('   - ' + s + '\n')

    output.write('\n')

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
        pymol_output.write('color white, resi ' + str(no_dot_list[0]))
        if len(no_dot_list) > 1:
            for i in no_dot_list[1:]:
                pymol_output.write('+' + str(i))
        pymol_output.write('\n')

    return filename

