"""
ProteoSync controller class.
"""

import AlignGUI
import TaxFileManager
import FileManagement
import BlastSearch
import StructSearch
import ClustalAlignment
import traceback
import shutil

import os
from os.path import exists
import certifi
os.environ["SSL_CERT_FILE"] = certifi.where()
base_path = 'Desktop/ProteoSync'


class AlignController:

    def __init__(self):
        self.successes = {}
        self.fails_i = {}
        self.fails_l = {}
        self.i_threshold = 0
        self.len_threshold = 0
        self.alignment_file = ''
        self.struct = ''
        self.pdb_codes = []
        self.alpha_struct = ''
        self.uniprot = ''
        self.struc_list = []
        self.exclude_list = []

    def update_database(self) -> int:
        """
        Updates the local BLAST-formatted PDB sequence database.

        Returns:
            -   0 if function executed successfully
            -   1 if an error occurred
        """
        try:
            FileManagement.update_database()
        except Exception as error:
            print("Error updating database:\n")
            # print(error)
            traceback.print_exc()
            print('\n')
            return 1

        return 0

    def run_blast(self, sequence: str, tax_tree: TaxFileManager.TaxTreeNode, i_threshold: int, len_threshold: int) -> int:
        """
        Runs a blast search with the given sequence against all species databases selected in the given tax tree.
        Stores the results in this object.

        Parameters:
            sequence: str, Query sequence to search against
            tax_tree: TaxTreeNode, Root node of the tax tree of databases to search
            i_threshold: int, % identity threshold for filtering BLAST results
            len_threshold: int, % length threshold for filtering BLAST results

        Returns:
            -   0 if function was successful
            -   1 if no hits were found from the BLAST search
            -   2 if an error occurred
        """

        print("Running BLAST searches...")

        seq_file = open(base_path+'/temp_files/seq.txt', 'w')
        seq_file.truncate(0)
        seq_file.write(sequence)
        seq_file.close()

        # Collect list of paths to search from tax settings
        path_list, exclude_list = TaxFileManager.get_path_list(tax_tree)
        self.exclude_list = exclude_list

        # Runs local BLAST searches on species databases
        try:
            successes, fails_i, fails_l = BlastSearch.blast_search(base_path+'/temp_files/seq.txt', i_threshold, len_threshold, path_list)
        except Exception as error:
            print("Error running BLAST searches:\n")
            traceback.print_exc()
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
        self.fails_i = fails_i
        self.fails_l = fails_l
        self.i_threshold = i_threshold
        self.len_threshold = len_threshold

        return 0

    def run_alignment(self) -> None:
        """Runs an alignment on the output file of run_blast. Stores the name of the alignment file in this object."""
        print("Aligning species...")
        self.alignment_file = ClustalAlignment.get_alignment(base_path+'/temp_files/formatted_results.txt')

    def run_pdb_search(self) -> int:
        """
        Runs a pdb structure search on the sequence file and concatenates a string representing the protein
        structure. Stores the results in this object.

        Returns:
            -   0 if function executed successfully
            -   1 if an error occurred
        """
        print("Searching pdb database...")

        try:
            struc_list = StructSearch.structure_search(base_path+'/temp_files/seq.txt')
        except Exception as error:
            print("Error running PDB search:\n")
            traceback.print_exc()
            print('\n')
            return 1

        self.struc_list = struc_list
        self.pdb_codes = [s[0] for s in struc_list]

        return 0

    def run_alpha_search(self, uniprot_id: str) -> int:
        """Runs an alphafold search on the given uniprot id. Stores the result in this object.

        Returns:
            -   0 if function executed successfully
            -   1 if an error occurred
        """
        print("Analyzing AlphaFold model...")
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
        """
        Assembles output file from other function outputs.
        Returns:
            Path to output file.
        """
        try:
            output_name = FileManagement.assemble_output_file(self.alignment_file, self.i_threshold, self.len_threshold,
                                                              self.successes, self.fails_i, self.fails_l,
                                                              self.struc_list, self.exclude_list, self.alpha_struct,
                                                              filename)

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
        self.i_threshold = 0
        self.alignment_file = ''
        self.struct = ''
        self.pdb_codes = []
        self.alpha_struct = ''


if __name__ == "__main__":
    cont = AlignController()
    gui = AlignGUI.AlignGUI(cont)
