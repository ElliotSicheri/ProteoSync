"""
ProteoSync controller class

ProteoSync helps identify unconserved segments in a protein sequence by:
    -   BLAST searching the sequence against several small species-specific databases
    -   Creating a clustalw sequence alignment with the top hits from each search
    -   Finding the closest ortholog to the sequence from the PBD database and aligning its secondary structure
    -   Finding the AlphaFold model for the sequence and aligning its secondary structure
"""

import AlignGUI
import TaxFileManager
import traceback
import FileManagement
import BlastSearch
import StructSearch
import ClustalAlignment
import shutil

import os
from os.path import exists
import certifi
os.environ["SSL_CERT_FILE"] = certifi.where()
base_path = 'Desktop/ProteoSync'


class AlignController:

    def __init__(self):
        self.successes = {}
        self.fails = {}
        self.threshold = 0
        self.len_threshold = 0
        self.alignment_file = ''
        self.struct = ''
        self.pdb_codes = []
        self.alpha_struct = ''
        self.uniprot = ''
        self.struc_list = []
        self.exclude_list = []

    def update_database(self) -> int:
        """Updates the local PDB database. Returns an error code.
        Error codes:
        0: Function executed successfully
        1: an error occurred"""
        try:
            FileManagement.update_database()
        except Exception as error:
            print("Error updating database:\n")
            # print(error)
            traceback.print_exc()
            print('\n')
            return 1

        return 0

    def run_blast(self, sequence: str, tax_tree: TaxFileManager.TaxTreeNode, threshold: int = 40, len_threshold: int = 40) -> int:
        """
        Runs a blast search with the given sequence against all species databases selected in the given tax tree.
        Stores the results in this object.

        Returns an error code:
            0: Function was successful
            1: No hits were found
            2: Other error
        """

        print("Running BLAST searches...\n")

        seq_file = open(base_path+'/temp_files/seq.txt', 'w')
        seq_file.truncate(0)
        seq_file.write(sequence)
        seq_file.close()

        # Collect list of paths to search from tax settings
        path_list, exclude_list = TaxFileManager.get_path_list(tax_tree)
        self.exclude_list = exclude_list

        # Runs local BLAST searches on species databases
        try:
            successes, fails = BlastSearch.blast_search(base_path+'/temp_files/seq.txt', threshold, len_threshold, path_list)
        except Exception as error:
            print("Error running BLAST searches:\n")
            traceback.print_exc()
            # print(type(error), error)
            print('\n')
            return 2

        if len(successes) <= 1:
            return 1

        # If no hits were found, updates the user and terminates
        blast_file = open(base_path+'/temp_files/formatted_results.txt')
        blast_results = blast_file.read()
        if blast_results == '':
            return 1

        self.successes = successes
        self.fails = fails
        self.threshold = threshold
        self.len_threshold = len_threshold

        return 0

    def run_alignment(self) -> None:
        """Runs an alignment on the output file of run_blast. Stores the name of the alignment file in this object."""
        print("Aligning species...\n")
        self.alignment_file = ClustalAlignment.get_alignment(base_path+'/temp_files/formatted_results.txt')

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
            struc_list = StructSearch.structure_search(base_path+'/temp_files/seq.txt')
        except Exception as error:
            print("Error running PDB search:\n")
            traceback.print_exc()
            # print(type(error), error)
            print('\n')
            return 1

        self.struc_list = struc_list
        self.pdb_codes = [s[0] for s in struc_list]

        return 0

    def run_alpha_search(self, uniprot_id: str) -> int:
        """Runs an alphafold search on the given uniprot id. Stores the result in this object.

        Returns an error code:
        0: Function returned successfully
        1: An error occurred
        """
        print("Analyzing AlphaFold model...\n")
        try:
            self.alpha_struct = StructSearch.alpha_struc_search(base_path+'/temp_files/seq.txt', uniprot_id)
            self.uniprot = uniprot_id
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
            output_name = FileManagement.assemble_output_file(self.alignment_file, self.threshold, self.len_threshold,
                                                              self.successes, self.fails, self.struc_list,
                                                              self.exclude_list, self.alpha_struct, filename)

            # Copy over structure files, if they exist
            if exists(base_path+'/downloads/AF_structures/AF-' + self.uniprot + ".pdb"):
                shutil.copy(base_path+'/downloads/AF_structures/AF-' + self.uniprot + ".pdb",
                            base_path+'/output/' + output_name)

            for code in self.pdb_codes:
                if exists(base_path+'/downloads/pdb_structures/pdb' + code + ".pdb"):
                    shutil.copy(base_path+'/downloads/pdb_structures/pdb' + code + ".pdb",
                                base_path+'/output/' + output_name)

            print('Results recorded in ' + output_name + '\n')

            return output_name
        except Exception as error:
            print('Error parsing output:\n')
            traceback.print_exc()
            # print(type(error), error)
            print('\n')
            return ''

    def get_fastas_from_uniprots(self, uniprots: list[str]) -> list[str]:
        return FileManagement.get_fastas_from_uniprots(uniprots)

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
