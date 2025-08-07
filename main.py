from Code.AlignController import AlignController
from Code.AlignGUI import AlignGUI
from os.path import exists
import os

if __name__ == "__main__":
    if not exists('./downloads'):
        os.mkdir('./downloads')
        os.mkdir('./downloads/fastas')
        os.mkdir('./downloads/pdb_structures')
        os.mkdir('./downloads/AF_structures')

    if not exists('./output'):
        os.mkdir('./output')

    if not exists('./temp_files'):
        os.mkdir('./temp_files')

    cont = AlignController()
    gui = AlignGUI(cont)
