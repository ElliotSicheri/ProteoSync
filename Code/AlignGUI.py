"""
ProteoSync GUI classes.
"""

import Code.AlignController as AlignController
import Code.TaxFileManager as TaxFileManager
import tkinter
from tkinter import filedialog as fd
from datetime import date

base_path = '.'

color_schemes = {
    'Blues': ['#bfbfff', '#0080ff', '#191999'],
    'Reds': ['#f7bcbc', '#ff5252', '#8a1313'],
    'Greens': ['#a6e6a6', '#51a657', '#0a4a0e'],
    'Grays': ['#b3b3b3', '#666666', '#000000'],
    'High contrast': ['#00ffff', '#ffff00', '#ff0000']
}


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


# def _validate_numeric_input(text):
#     return text.isdecimal() or text == ""  # Allow empty string for backspace/delete


class AlignGUI:
    def __init__(self, aln_cont: AlignController):
        self.aln_cont = aln_cont
        self.line_count = 1.0

        window = tkinter.Tk()
        window.title('ProteoSync')
        window.geometry('730x570')
        window['bg'] = '#363535'
        window.resizable(width=False, height=False)
        self.window = window

        self.tax_window = None
        self.color_window = None
        # self.proj_window = None
        self.file_tree = None
        # self.proj_seq_entry = None
        # self.num_lines_entry = None
        #
        # self.proj_file = tkinter.StringVar(value='Select file...')
        self.color = tkinter.StringVar(value='Blues')

        seq_label = tkinter.Label(text='Protein sequences:', bg='#363535', fg='white')
        seq_label.place(x=10, y=5)
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

        high_thresh_label = tkinter.Label(text='Upper % Identity Threshold:', bg='#363535', fg='white')
        high_thresh_label.place(x=10, y=217)
        high_threshold_slider = tkinter.Scale(from_=0, to=100, orient='horizontal', length=250, bg='#363535', fg='white')
        high_threshold_slider.place(x=210, y=200)
        high_threshold_slider.set(90)
        self.high_threshold_slider = high_threshold_slider

        low_thresh_label = tkinter.Label(text='Lower % Identity Threshold:', bg='#363535', fg='white')
        low_thresh_label.place(x=10, y=257)
        low_threshold_slider = tkinter.Scale(from_=0, to=100, orient='horizontal', length=250, bg='#363535', fg='white')
        low_threshold_slider.place(x=210, y=240)
        low_threshold_slider.set(50)
        self.low_threshold_slider = low_threshold_slider

        len_thresh_label = tkinter.Label(text='% Length Variability Threshold:', bg='#363535', fg='white')
        len_thresh_label.place(x=10, y=297)
        len_threshold_slider = tkinter.Scale(from_=0, to=100, orient='horizontal', length=250, bg='#363535', fg='white')
        len_threshold_slider.place(x=210, y=280)
        len_threshold_slider.set(50)
        self.len_threshold_slider = len_threshold_slider

        output_name_label = tkinter.Label(text='Output file name (optional):', bg='#363535', fg='white')
        output_name_label.place(x=10, y=330)
        output_name_entry = tkinter.Entry(fg='black', bg='white', width=33, highlightbackground='#363535')
        output_name_entry.place(x=185, y=330)
        self.output_name_entry = output_name_entry

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

        start_button = tkinter.Button(text='Start!', width=8, height=2, bg='#616161', fg='black',
                                      highlightbackground='#363535', command=self.run_program)
        start_button.place(x=10, y=365)

        tax_settings_button = tkinter.Button(text='Taxonomy settings', width=12, height=2, bg='#616161', fg='black',
                                             highlightbackground='#363535', command=self.open_tax_window)
        tax_settings_button.place(x=130, y=365)

        color_settings_button = tkinter.Button(text='Color settings', width=10, height=2, bg='#616161', fg='black',
                                               highlightbackground='#363535', command=self.open_color_window)
        color_settings_button.place(x=285, y=365)

        # projection_button = tkinter.Button(text='Custom projection', width=11, height=2, bg='#616161', fg='black',
        #                                    highlightbackground='#363535', command=self.open_custom_projection_window)
        # projection_button.place(x=420, y=365)

        output_log_label = tkinter.Label(text='Output Log:', bg='#363535', fg='white')
        output_log_label.place(x=10, y=406)
        output_field = tkinter.Text(fg='black', bg='white', width=100, height=10, highlightbackground='#363535')
        output_field.place(x=10, y=425)
        self.output_field = output_field

        self.open_tax_window()

        window.mainloop()

    def open_tax_window(self):
        if self.tax_window is not None and self.tax_window.winfo_exists():
            self.tax_window.lift()
            self.tax_window.focus_force()
            return

        self.file_tree = TaxFileManager.make_tax_tree(base_path + '/databases/species')
        height = self.file_tree.get_size()

        tax_window = tkinter.Toplevel(self.window)
        tax_window['bg'] = '#363535'
        tax_window.geometry('328x700')
        tax_window.title('Taxonomy settings')
        self.tax_window = tax_window

        button_frame = tkinter.Frame(tax_window)

        select_button = tkinter.Button(button_frame, text='Select all', width=7, height=1, bg='#616161', fg='black',
                                       highlightbackground='#363535', command=self.select_all)
        select_button.pack(side=tkinter.LEFT)

        deselect_button = tkinter.Button(button_frame, text='Deselect all', width=7, height=1, bg='#616161', fg='black',
                                         highlightbackground='#363535', command=self.deselect_all)
        deselect_button.pack(side=tkinter.LEFT)

        button_frame.pack()

        scrollbar = tkinter.Scrollbar(tax_window)
        scrollbar.pack(side=tkinter.RIGHT, fill=tkinter.Y)

        canvas = tkinter.Canvas(tax_window, width=500, height=700, bg='#363535', yscrollcommand=scrollbar.set,
                                yscrollincrement=1, scrollregion=(0, 0, 500, height * 20 + 5), borderwidth=0,
                                highlightthickness=0)
        canvas.pack(side=tkinter.LEFT, fill=tkinter.BOTH)
        scrollbar.config(command=canvas.yview)
        # self.canvas = canvas

        frame = tkinter.Frame(canvas, bg='#363535', width=500, height=height * 20 + 10, borderwidth=0,
                              highlightthickness=0)
        self.frame = frame
        frame.pack(side=tkinter.TOP)

        height, width, root = self.place_checkboxes(self.file_tree, 5, 5)
        self.root_box = root
        self.frame.configure(width=width)
        canvas.create_window(0, 0, anchor='nw', window=frame)

        tax_window.geometry(str(width + scrollbar.winfo_width() + 10) + "x" + str(tax_window.winfo_height()))
        tax_window.update()

    def open_color_window(self):
        if self.color_window is not None and self.color_window.winfo_exists():
            self.color_window.lift()
            self.color_window.focus_force()
            return

        color_window = tkinter.Toplevel(self.window)
        color_window['bg'] = '#363535'
        color_window.geometry('250x' + str(30 * len(color_schemes)))
        color_window.title('Color settings')
        self.color_window = color_window

        for scheme_name, colors in color_schemes.items():
            row = tkinter.Frame(color_window, bg='#363535')
            row.pack(anchor='w', pady=2, padx=5)

            rb = tkinter.Radiobutton(row, text=scheme_name, variable=self.color, value=scheme_name, bg='#363535',
                                     fg='white', selectcolor='#363535')
            rb.pack(side='left')

            for color in colors:
                canvas = tkinter.Canvas(row, width=15, height=15, highlightthickness=1, highlightbackground='black')
                canvas.create_rectangle(0, 0, 15, 15, fill=color, outline=color)
                canvas.pack(side='left', padx=2)

    # def open_custom_projection_window(self):
    #     if self.proj_window is not None and self.proj_window.winfo_exists():
    #         self.proj_window.lift()
    #         self.proj_window.focus_force()
    #         return
    #
    #     proj_window = tkinter.Toplevel(self.window)
    #     proj_window['bg'] = '#363535'
    #     proj_window.geometry('500x200')
    #     proj_window.title('Custom alignment projection')
    #     self.proj_window = proj_window
    #
    #     frame = tkinter.Frame(proj_window, bg='#363535')
    #     frame.pack(anchor='w', pady=2, padx=5)
    #     select_button = tkinter.Button(frame, text='Select file', width=7, height=1, bg='#616161', fg='black',
    #                                    highlightbackground='#363535', command=self._select_proj_file)
    #     select_button.pack(side=tkinter.LEFT)
    #     file_label = tkinter.Label(frame, textvariable=self.proj_file)
    #     file_label.pack(side=tkinter.LEFT)
    #
    #     frame2 = tkinter.Frame(proj_window, bg='#363535')
    #     frame2.pack(anchor='w', pady=2, padx=5)
    #     seq_label = tkinter.Label(frame2, text='Sequence to project to: ')
    #     seq_label.pack(side=tkinter.LEFT)
    #     proj_seq_entry = tkinter.Entry(frame2, fg='black', bg='white', width=33, highlightbackground='#363535')
    #     self.proj_seq_entry = proj_seq_entry
    #     proj_seq_entry.pack(side=tkinter.LEFT)
    #
    #     frame3 = tkinter.Frame(proj_window, bg='#363535')
    #     frame3.pack(anchor='w', pady=2, padx=5)
    #     num_label = tkinter.Label(frame3, text='Number of lines in projection: ')
    #     num_label.pack(side=tkinter.LEFT)
    #     vcmd = frame3.register(_validate_numeric_input)
    #     num_lines_entry = tkinter.Entry(frame3, fg='black', bg='white', width=5, highlightbackground='#363535',
    #                                     validate="key", validatecommand=(vcmd, '%P'))
    #     self.num_lines_entry = num_lines_entry
    #     num_lines_entry.pack(side=tkinter.LEFT)
    #
    #     frame4 = tkinter.Frame(proj_window, bg='#363535')
    #     frame4.pack(anchor='w', pady=2, padx=5)
    #     run_button = tkinter.Button(frame4, text='Create projection', width=10, height=1, bg='#616161', fg='black',
    #                                    highlightbackground='#363535', command=self.do_custom_projection)
    #     run_button.pack(side=tkinter.LEFT)
    #
    # def _select_proj_file(self):
    #     self.proj_file.set(fd.askopenfilename())

    def run_program(self) -> None:
        """Collects user input from window fields and passes it to the controller."""
        self.line_count = 1.0
        self.output_field.delete(1.0, tkinter.END)

        seq_str = self.get_seq()
        uni_str = self.get_uniprot()
        low_threshold = self.get_low_threshold()
        high_threshold = self.get_high_threshold()
        len_threshold = self.get_len_threshold()
        color_scheme = self.get_color_scheme()

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
                self._run_program(seq_lst[i], uni_lst[i], low_threshold, high_threshold, len_threshold, run_name, color_scheme)
            else:
                self._run_program(seq_lst[i], '', low_threshold, high_threshold, len_threshold, run_name, color_scheme)

    def _run_program(self, seq_str: str, uni_str: str, low_threshold: int, high_threshold: int, len_threshold: int, run_name: str,
                     color_scheme: str) -> None:
        self.aln_cont.clear()

        self.printout("Running BLAST searches...\n")
        error = self.aln_cont.run_blast(seq_str, self.file_tree, low_threshold, high_threshold, len_threshold)
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

        output_name = self.aln_cont.assemble_output(run_name, color_scheme)

        if output_name != '':
            self.printout('Results recorded in ' + output_name + '\n')
        else:
            self.printout('There has been an error, output could not be parsed. See console window for more details.\n')

    # def do_custom_projection(self):
    #     filename = self.get_proj_file()
    #     output_file = base_path + '/output' + filename[filename.rindex('/'):filename.rindex('.')] + '_projection.txt'
    #     proj_seq = self.get_proj_file()
    #     color = self.get_color_scheme()
    #
    #     AlignController.do_custom_projection(filename, output_file, proj_seq, color)

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

    def get_low_threshold(self) -> int:
        """Return the current setting of the threshold slider."""
        return int(self.low_threshold_slider.get())

    def get_high_threshold(self) -> int:
        """Return the current setting of the threshold slider."""
        return int(self.high_threshold_slider.get())

    def get_len_threshold(self) -> int:
        """Return the current setting of the length threshold slider."""
        return int(self.len_threshold_slider.get())

    def get_uniprot(self) -> str:
        """Return the contents of the Uniprot code entry field."""
        return self.uniprot_entry.get(1.0, tkinter.END)

    def get_filename(self) -> str:
        return self.output_name_entry.get()

    def get_color_scheme(self) -> str:
        return self.color.get()

    def get_proj_file(self):
        return self.proj_file.get()

    def get_proj_seq(self):
        return self.proj_seq_entry.get()

    def get_num_seqs(self):
        return int(self.num_lines_entry.get())

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
