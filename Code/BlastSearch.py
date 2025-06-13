"""
Contains functions for performing BLAST searches against species database sets
"""

from Bio.Blast.Applications import NcbiblastpCommandline
import os

base_path = 'Desktop/ProteoSync'


def blast_search(seq_file: str, i_threshold: int, len_threshold: int, path_list: list[str]) \
        -> (dict[str: str], dict[str: int], dict[str: (int, int)]):
    """
    BLAST searches the query sequence from seq_file against a set of BLAST-formatted databases given by path_list.
    Writes the top hit sequence from each database in an output file if they have a % identity above the given identity
    threshold and have % length relative to the query sequence within the given length threshold.

    Parameters:
        seq_file: str, the path to the file that contains the query sequence
        i_threshold: int, % identity threshold for which BLAST hits are included
        len_threshold: int, % length threshold with respect to the length query sequence
        path_list: list[str], list of paths which contain databases to be searched

    Returns:
        -   dict of datasets that returned a hit above the thresholds mapped to the accession code of the hit sequence
        -   dict of datasets that did not return a hit above the identity threshold, mapped to their % identity
        -   dict of datasets that returned a hit above the identity threshold but not within the length threshold,
            mapped to % identity and % length with respect to the query sequence
    """

    formatted_results = ''  # Properly formatted result string
    successful_searches = {}  # Species that returned a hit above the threshold
    failed_iden = {}  # Species that did not return a hit above the threshold
    failed_len = {}  # Species that passed % identity threshold, but not length threshold

    low_lim = 100 - len_threshold
    high_lim = 100 + len_threshold

    with open(seq_file) as f:
        query_len = len(f.read().replace('\n', ''))
        f.close()

    for path in path_list:
        blast = NcbiblastpCommandline(cmd='blastp', query=seq_file, db=path + '/' + os.path.basename(path) + '.fasta',
                                      out=base_path+'/temp_files/raw_data.txt', outfmt=0,
                                      num_threads=8, num_alignments=5, num_descriptions=5)
        blast()
        sequence, code, score = _return_first_hit()
        species = os.path.basename(path)
        len_score = int((len(sequence) / query_len) * 100)

        if score < i_threshold:
            failed_iden[species] = score
        elif low_lim > len_score or len_score > high_lim:
            failed_len[species] = (score, len_score)
        else:
            formatted_results += '> ' + species + '_[' + str(score) + '] \n'
            formatted_results += sequence + '\n'
            successful_searches[species] = code

    with open(base_path+'/temp_files/formatted_results.txt', 'w') as file:
        # Clears the contents of formatted_results.txt and writes the formatted data to the file
        file.truncate(0)
        file.write(formatted_results)

    return successful_searches, failed_iden, failed_len


def _return_first_hit() -> (str, str, int):
    """
    Returns the sequence, accession code and % identity of the first hit in the BLAST search output file.
    """

    hit = ''
    code = ''
    score = 0
    started = False
    segment_dict = {}

    with open(base_path+'/temp_files/raw_data.txt', 'r') as file:
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
