"""
Controller class for alignment manager
"""

import AlignGUI
import TaxFileManager
import traceback
import alignment_manager_v3

import os
import certifi
os.environ["SSL_CERT_FILE"] = certifi.where()

class AlignController:

    def __init__(self):
        self.successes = {}
        self.fails = {}
        self.threshold = 0
        self.alignment_file = ''
        self.struct = ''
        self.pdb_codes = []
        self.alpha_struct = ''
        self.struc_list = []

    def update_database(self) -> int:
        """Updates the local PDB database. Returns an error code.
        Error codes:
        0: Function executed successfully
        1: an error occurred"""
        try:
            alignment_manager_v3.update_database()
        except Exception as error:
            print("Error updating database:\n")
            # print(error)
            traceback.print_exc()
            print('\n')
            return 1

        return 0

    def run_blast(self, sequence: str, tax_tree: TaxFileManager.TaxTreeNode, threshold: int = 40) -> int:
        """
        Runs a blast search with the given sequence against all species databases selected in the given tax tree.
        Stores the results in this object.

        Returns an error code:
            0: Function was successful
            1: No hits were found
            2: Other error
        """

        print("Running BLAST searches...\n")

        seq_file = open('Desktop/ProteoSync/temp_files/seq.txt', 'w')
        seq_file.truncate(0)
        seq_file.write(sequence)
        seq_file.close()

        # Collect list of paths to search from tax settings
        path_list = TaxFileManager.get_path_list(tax_tree)

        # Runs local BLAST searches on species databases
        try:
            successes, fails = alignment_manager_v3.blast_search('Desktop/ProteoSync/temp_files/seq.txt',
                                                                 threshold, path_list)
        except Exception as error:
            print("Error running BLAST searches:\n")
            traceback.print_exc()
            # print(type(error), error)
            print('\n')
            return 2

        if len(successes) <= 1:
            return 1

        # If no hits were found, updates the user and terminates
        blast_file = open('Desktop/ProteoSync/temp_files/formatted_results.txt')
        blast_results = blast_file.read()
        if blast_results == '':
            return 1

        self.successes = successes
        self.fails = fails
        self.threshold = threshold

        return 0

    def run_alignment(self) -> None:
        """Runs an alignment on the output file of run_blast. Stores the name of the alignment file in this object."""
        print("Aligning species...\n")
        self.alignment_file = alignment_manager_v3.get_alignment('Desktop/ProteoSync/temp_files/formatted_results.txt')

    def run_pdb_search(self) -> int:
        """
        Runs a pdb structure search on the sequence file and concatenates a string representing the protein
        structure. Stores the results in this object.

        Returns an error code:
        0: Function executed successfully.
        1: An error occurred.
        """
        print("Searching pdb database...\n")

        try:
            struc_list = alignment_manager_v3.structure_search('Desktop/ProteoSync/temp_files/seq.txt')
        except Exception as error:
            print("Error running PDB search:\n")
            traceback.print_exc()
            # print(type(error), error)
            print('\n')
            return 1

        self.struc_list = struc_list

        return 0

    def run_alpha_search(self, uniprot_id: str) -> int:
        """Runs an alphafold search on the given uniprot id. Stores the result in this object.

        Returns an error code:
        0: Function returned successfully
        1: An error occurred
        """
        print("Analyzing AlphaFold model...\n")
        try:
            self.alpha_struct = alignment_manager_v3.alpha_struc_search('Desktop/ProteoSync/temp_files/seq.txt',
                                                                        uniprot_id)
        except Exception as error:
            print('Error running AlphaFold analysis:')
            traceback.print_exc()
            # print(type(error), error)
            print('\n')
            return 1

        return 0

    def assemble_output(self, filename: str) -> str:
        """Assembles output file from other function outputs. Returns output file name."""
        try:
            output_name = alignment_manager_v3.assemble_output_file(self.alignment_file, self.threshold, self.successes,
                                                                    self.fails, self.struc_list, self.alpha_struct, filename)
            print('Results recorded in ' + output_name + '\n')
            return output_name
        except Exception as error:
            print('Error parsing output:\n')
            traceback.print_exc()
            # print(type(error), error)
            print('\n')
            return ''

    def clear(self):
        """Clears results stored in this object."""
        self.successes = {}
        self.fails = {}
        self.threshold = 0
        self.alignment_file = ''
        self.struct = ''
        self.pdb_codes = []
        self.alpha_struct = ''


if __name__ == "__main__":
    cont = AlignController()
    gui = AlignGUI.AlignGUI(cont)
