from Bio.Blast.Applications import NcbiblastpCommandline
import os

base_path = 'Desktop/ProteoSync'


def blast_search(seq_file: str, threshold: int, len_threshold: int, path_list: list[str]) -> (dict[str: str], dict[str: int]):
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

    low_lim = 1 - len_threshold / 100
    high_lim = 1 + len_threshold / 100

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
        if score > threshold and low_lim < len(sequence) / query_len < high_lim:
            formatted_results += '> ' + species + '_[' + str(score) + '] \n'
            formatted_results += sequence + '\n'
            successful_searches[species] = code
        else:
            failed_searches[species] = score

    with open(base_path+'/temp_files/formatted_results.txt', 'w') as file:
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
