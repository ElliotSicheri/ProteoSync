"""
Contains functions for performing BLAST searches against species database sets
"""

from Bio.Blast.Applications import NcbiblastpCommandline
import os

base_path = '.'


def blast_search(seq_file: str, i_threshold_low: int, i_threshold_high: int, len_threshold: int, path_list: list[str],
                 num_hits: int = 1) \
        -> (dict[str: str], dict[str: int], dict[str: (int, int)]):
    """
    BLAST searches the query sequence from seq_file against a set of BLAST-formatted databases given by path_list.
    Writes the top hit sequence from each database in an output file if they have a % identity above the given identity
    threshold and have % length relative to the query sequence within the given length threshold.

    Parameters:
        seq_file: str, the path to the file that contains the query sequence
        i_threshold_low: int, lower % identity threshold for which BLAST hits are included
        i_threshold_high: int, upper % identity threshold for which BLAST hits are included
        len_threshold: int, % length threshold with respect to the length query sequence
        path_list: list[str], list of paths which contain databases to be searched
        num_hits: int, max number of sequences to include from each species.

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
        query_seq = f.read().replace('\n', '')
        query_len = len(query_seq)
        f.close()

    formatted_results += '> QUERY_SEQUENCE \n'
    formatted_results += query_seq + '\n'

    for path in path_list:
        blast = NcbiblastpCommandline(cmd='blastp', query=seq_file, db=path + '/' + os.path.basename(path) + '.fasta',
                                      out=base_path+'/temp_files/raw_data.txt', outfmt=0,
                                      num_threads=8, num_alignments=30, num_descriptions=30)
        blast()
        hits_list = _return_first_hits(num_hits)
        species = os.path.basename(path)
        i = 0
        hits_list = sorted(hits_list, key=lambda x: x[-1], reverse=True)
        for sequence, code, score in hits_list:
            len_score = int((len(sequence) / query_len) * 100)

            s_tag = species
            if i > 0:
                s_tag += f'_{i + 1}'

            if score < i_threshold_low or score > i_threshold_high:
                failed_iden[s_tag] = score
            elif low_lim > len_score or len_score > high_lim:
                failed_len[s_tag] = (score, len_score)
            else:
                formatted_results += '> ' + s_tag + '_[' + str(score) + '] \n'
                formatted_results += sequence + '\n'
                successful_searches[s_tag] = code
            i += 1

    with open(base_path+'/temp_files/formatted_results.txt', 'w') as file:
        # Clears the contents of formatted_results.txt and writes the formatted data to the file
        file.truncate(0)
        file.write(formatted_results)

    return successful_searches, failed_iden, failed_len


def _return_first_hits(num_hits: int = 1) -> list[(str, str, int)]:
    """
    Returns the sequence, accession code and % identity of the first num_hits hits in the BLAST search output file.
    """

    hit_list = []

    hits = 0

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
            if hits == num_hits:
                keys = list(segment_dict.keys())
                while len(keys) > 0:
                    smallest = min(keys)
                    hit += segment_dict[smallest]
                    keys.remove(smallest)
                hit_list.append((hit.replace('-', ''), code, score))
                break

            if started:
                keys = list(segment_dict.keys())
                while len(keys) > 0:
                    smallest = min(keys)
                    hit += segment_dict[smallest]
                    keys.remove(smallest)
                hit_list.append((hit.replace('-', ''), code, score))
                hit = ''
                code = ''
                score = 0
                segment_dict = {}

            c_start = line.find(' ')
            c_start = line.find(' ', c_start + 1)
            code = line[2:c_start]
            started = True

            hits += 1
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
            return [('', '', 0)]

    if hits != num_hits:
        # We've reached the end of the file. append the last hit.
        keys = list(segment_dict.keys())
        while len(keys) > 0:
            smallest = min(keys)
            hit += segment_dict[smallest]
            keys.remove(smallest)
        hit_list.append((hit.replace('-', ''), code, score))

    return hit_list
