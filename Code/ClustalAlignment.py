from Bio.Align.Applications import ClustalwCommandline


def get_alignment(input_file: str) -> str:
    """Creates a sequence alignment from the sequences in the input file, then writes it to an output file.

    Parameters:
        input_file: str, path to the file containing the sequences to be aligned.

    Returns:
        Path to the output file.
    """

    cmd = ClustalwCommandline('clustalw', infile=input_file, seqnos='ON')
    stdout, stderr = cmd()
    # Extracts the name of the output file from the readout
    s = stdout.rfind('[')
    e = stdout.rfind(']')
    output_file = stdout[s + 1:e]
    # Returns the aligned output file
    return output_file
