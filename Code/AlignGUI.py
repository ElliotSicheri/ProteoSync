"""
GUI for running Alignment Manager.
"""

import AlignController
import TaxFileManager
from datetime import date
import tkinter

base_path = 'Desktop/ProteoSync'
# base_path = '.'


class CheckboxTree:
    """Tree-based structure containing and connecting the checkbuttons on the taxonomy GUI."""

    def __init__(self, checkbox: tkinter.Checkbutton, node: TaxFileManager.TaxTreeNode):
        checkbox.configure(command=self.on_click)

        self.checkbox = checkbox
        self.children = []
        self.node = node

    def on_click(self):
        if self.node.include.get():
            self.child_select()
        else:
            self.child_deselect()

    def child_select(self):
        self.checkbox.select()
        for child in self.children:
            child.checkbox.configure(state="normal")
            child.child_select()

    def child_deselect(self):
        self.checkbox.deselect()
        for child in self.children:
            child.child_deselect()
            child.checkbox.configure(state="disabled")


class AlignGUI:
    def __init__(self, aln_cont: AlignController):
        self.aln_cont = aln_cont
        self.line_count = 1.0

        window = tkinter.Tk()
        window.title('ProteoSync')
        window.geometry('730x520')
        window['bg'] = '#363535'
        window.resizable(width=False, height=False)
        self.window = window

        self.window2 = None
        self.file_tree = None

        seq = tkinter.Label(text='Protein sequences:', bg='#363535', fg='white')
        seq.place(x=10, y=5)

        seq_entry = tkinter.Text(fg='black', bg='white', width=49, height=11, highlightbackground='#363535')
        seq_entry.place(x=10, y=30)
        self.seq_entry = seq_entry

        uniprot_label = tkinter.Label(text='OR UniProt IDs:', bg='#363535', fg='white')
        uniprot_label.place(x=370, y=5)

        uniprot_entry = tkinter.Text(fg='black', bg='white', width=49, height=11, highlightbackground='#363535')
        uniprot_entry.place(x=370, y=30)
        self.uniprot_entry = uniprot_entry

        pdb = tkinter.IntVar()
        self.pdb = pdb
        pdb_check = tkinter.Checkbutton(text='Search PDB database?', bg='#363535', fg='white', variable=self.pdb)
        pdb_check.place(x=10, y=180)
        pdb_check.select()

        thresh_label = tkinter.Label(text='% Identity Threshold:', bg='#363535', fg='white')
        thresh_label.place(x=10, y=217)

        threshold_slider = tkinter.Scale(from_=0, to=100, orient='horizontal', length=250, bg='#363535', fg='white')
        threshold_slider.place(x=210, y=200)
        threshold_slider.set(50)
        self.threshold_slider = threshold_slider

        len_thresh_label = tkinter.Label(text='% Length Difference Threshold:', bg='#363535', fg='white')
        len_thresh_label.place(x=10, y=257)

        len_threshold_slider = tkinter.Scale(from_=0, to=100, orient='horizontal', length=250, bg='#363535', fg='white')
        len_threshold_slider.place(x=210, y=240)
        len_threshold_slider.set(50)
        self.len_threshold_slider = len_threshold_slider

        output_name_label = tkinter.Label(text='Output file name (optional):', bg='#363535', fg='white')
        output_name_label.place(x=10, y=290)

        output_name_entry = tkinter.Entry(fg='black', bg='white', width=33, highlightbackground='#363535')
        output_name_entry.place(x=185, y=290)
        self.output_name_entry = output_name_entry

        start_button = tkinter.Button(text='Start!', width=8, height=2, bg='#616161', fg='black',
                                      highlightbackground='#363535', command=self.run_program)
        start_button.place(x=10, y=315)

        date_file = open(base_path+'/databases/pdb_last_update.txt', 'r')
        last_date_str = date_file.read()
        last_date = date(int(last_date_str[0:4]), int(last_date_str[5:7]), int(last_date_str[8:10]))
        today = date.today()
        delta = today - last_date
        days = str(delta.days)

        update_label = tkinter.Label(text='It has been ' + days + ' day(s) since the', bg='#363535', fg='white')
        update_label.place(x=505, y=190)

        update_label_2 = tkinter.Label(text='local PDB database has been updated.', bg='#363535', fg='white')
        update_label_2.place(x=470, y=210)

        pdb_button = tkinter.Button(text='Update PDB database', width=15, height=2, bg='#616161', fg='black',
                                highlightbackground='#363535', command=self.update_database)
        pdb_button.place(x=510, y=240)

        output_log_label = tkinter.Label(text='Output Log:', bg='#363535', fg='white')
        output_log_label.place(x=10, y=355)

        output_field = tkinter.Text(fg='black', bg='white', width=100, height=10, highlightbackground='#363535')
        output_field.place(x=10, y=375)
        self.output_field = output_field

        self.file_tree = TaxFileManager.make_tax_tree(base_path+'/databases/species')
        height = self.file_tree.get_size()

        window2 = tkinter.Toplevel(self.window)
        window2['bg'] = '#363535'
        window2.geometry('328x700')
        window2.title('Taxonomy settings')
        self.window2 = window2

        button_frame = tkinter.Frame(window2)

        select_button = tkinter.Button(button_frame, text='Select all', width=7, height=1, bg='#616161', fg='black',
                                       highlightbackground='#363535', command=self.select_all)
        select_button.pack(side=tkinter.LEFT)

        deselect_button = tkinter.Button(button_frame, text='Deselect all', width=7, height=1, bg='#616161', fg='black',
                                         highlightbackground='#363535', command=self.deselect_all)
        deselect_button.pack(side=tkinter.LEFT)

        button_frame.pack()

        scrollbar = tkinter.Scrollbar(window2)
        scrollbar.pack(side=tkinter.RIGHT, fill=tkinter.Y)

        canvas = tkinter.Canvas(window2, width=500, height=700, bg='#363535', yscrollcommand=scrollbar.set,
                                yscrollincrement=1, scrollregion=(0, 0, 500, height * 20 + 5), borderwidth=0,
                                highlightthickness=0)
        canvas.pack(side=tkinter.LEFT, fill=tkinter.BOTH)
        scrollbar.config(command=canvas.yview)
        self.canvas = canvas

        frame = tkinter.Frame(canvas, bg='#363535', width=500, height=height * 20 + 10, borderwidth=0,
                              highlightthickness=0)
        self.frame = frame

        frame.pack(side=tkinter.TOP)

        height, width, root = self.place_checkboxes(self.file_tree, 5, 5)
        self.root_box = root
        self.frame.configure(width=width)
        canvas.create_window(0, 0, anchor='nw', window=frame)

        window2.geometry(str(width+scrollbar.winfo_width()+10) + "x" + str(window2.winfo_height()))
        window2.update()

        window.mainloop()

    def run_program(self) -> None:
        """Collects user input from window fields and passes it to the controller."""
        self.line_count = 1.0
        self.output_field.delete(1.0, tkinter.END)

        seq_str = self.get_seq()
        uni_str = self.get_uniprot()
        threshold = self.get_threshold()
        len_threshold = self.get_len_threshold()

        uni_lst = []
        if uni_str != '\n':
            uni_lst = uni_str.split()
            seq_lst = self.aln_cont.get_fastas_from_uniprots(uni_lst)
        else:
            seq_lst = seq_str.split()

        filename = self.get_filename()

        for i in range(len(seq_lst)):
            # logic to decide what to name this run
            if filename != '':
                run_name = filename
                if i > 0:
                    run_name += '_' + str(i + 1)
            elif uni_str != '\n':
                run_name = uni_lst[i]
            else:
                run_name = ''

            if i > 0:
                self.printout('\n')

            if run_name == '':
                self.printout("Starting run " + str(i+1) + '\n')
                print("Starting run " + str(i+1))
            else:
                self.printout("Starting run " + run_name + '\n')
                print("Starting run " + run_name)

            if uni_str != '\n':
                self._run_program(seq_lst[i], uni_lst[i], threshold, len_threshold, run_name)
            else:
                self._run_program(seq_lst[i], '', threshold, len_threshold, run_name)

    def _run_program(self, seq_str: str, uni_str: str, threshold: int, len_threshold: int, run_name: str) -> None:
        self.aln_cont.clear()

        self.printout("Running BLAST searches...\n")
        error = self.aln_cont.run_blast(seq_str, self.file_tree, threshold, len_threshold)
        if error == 1:
            self.printout("1 or fewer sufficiently similar sequences were found, no alignment could be made.\n")
            return
        elif error == 2:
            self.printout("There was an error while attempting BLAST searches. See console window for more details.\n")
            return

        # Aligns the top hits
        self.printout('Creating alignment...\n')
        self.aln_cont.run_alignment()

        # Reads output alignment from file
        if self.pdb.get() == 1:
            self.printout('Searching PDB database...\n')
            error = self.aln_cont.run_pdb_search()
            if error == 1:
                self.printout("An error occurred while searching the PDB database. See console window for more "
                              "details.\n")
                self.printout("Continuing...\n")

        if uni_str != '':
            self.printout('Analyzing AlphaFold model...\n')
            error = self.aln_cont.run_alpha_search(uni_str)
            if error == 1:
                self.printout("An error occurred while analyzing the AlphaFold model. See console window for more "
                              "details.\n")
                self.printout("Continuing...\n")

        output_name = self.aln_cont.assemble_output(run_name)

        if output_name != '':
            self.printout('Results recorded in ' + output_name + '\n')
        else:
            self.printout('There has been an error, output could not be parsed. See console window for more details.\n')

    def place_checkboxes(self, file_node: TaxFileManager.TaxTreeNode, height: int, indent: int) \
            -> (int, int, CheckboxTree):
        """Places checkboxes on the taxonomy select screen"""
        curr_height = height
        checkbox = tkinter.Checkbutton(self.frame, bg='#363535', fg='white',
                                       variable=file_node.include)
        checkbox.place(x=indent, y=curr_height)
        checkbox.select()

        checktree = CheckboxTree(checkbox, file_node)

        label = tkinter.Label(self.frame, text=file_node.name, bg='#363535', fg='white')
        label.place(x=indent + checkbox.winfo_width() + 20, y=curr_height)

        self.frame.update_idletasks()
        max_width = label.winfo_width() + label.winfo_x()

        if not file_node.is_leaf:
            for child in file_node.children:
                curr_height += 20
                curr_height, width, child = self.place_checkboxes(child, curr_height, indent + 25)
                if width > max_width:
                    max_width = width
                checktree.children.append(child)

        return curr_height, max_width, checktree

    def update_database(self):
        self.line_count = 1.0
        self.output_field.delete(1.0, tkinter.END)
        self.printout('Updating local PDB database...\n')
        code = self.aln_cont.update_database()
        if code == 0:
            print("Done!\n")
            self.printout('Done!\n')
        elif code == 1:
            print('An error occurred while updating the PDB database.\n')
            self.printout('An error occurred while updating the PDB database.\n')

    def get_seq(self) -> str:
        """Return contents of the sequence entry field."""
        seq_str = self.seq_entry.get(1.0, tkinter.END)
        return seq_str.replace('\n', '')

    def get_threshold(self) -> int:
        """Return the current setting of the threshold slider."""
        return int(self.threshold_slider.get())

    def get_len_threshold(self) -> int:
        """Return the current setting of the length threshold slider."""
        return int(self.len_threshold_slider.get())

    def get_uniprot(self) -> str:
        """Return the contents of the Uniprot code entry field."""
        return self.uniprot_entry.get(1.0, tkinter.END)

    def get_filename(self) -> str:
        return self.output_name_entry.get()

    def printout(self, message: str):
        """Outputs a message to the user."""
        self.output_field.insert(self.line_count, message)
        self.line_count += 1
        self.window.update()

    def select_all(self):
        """Selects all checkboxes on the taxonomy menu."""
        self.root_box.child_select()

    def deselect_all(self):
        """Deselects all checkboxes on the taxonomy menu."""
        for child in self.root_box.children:
            child.child_deselect()


