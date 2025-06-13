"""
Contains functions for creating, manipulating, and interacting with database taxonomy trees.
"""

import os
import tkinter


class TaxTreeNode:
    """A node in a taxonomy tree"""

    def __init__(self, path: str):
        self.path = path
        self.name = os.path.basename(path)
        self.children = []
        self.is_leaf = False
        self.include = tkinter.IntVar()
        self.checkbox = None

    def get_size(self) -> int:
        """Returns the number of nodes in the tree."""
        size = 1
        if not self.is_leaf:
            for child in self.children:
                size += child.get_size()
        return size


def make_tax_tree(directory: str) -> TaxTreeNode:
    """
    Creates a taxonomy tree data structure from the directory structure of the given directory.

    Internal node directories should only contain directories. Leaf node directories should only contain files.
    It is assumed that a valid and properly structured directory pathway is passed to the function.

    Parameters:
        directory: str, Path to the root database directory

    Returns:
        Root node of the newly generated taxonomy tree structure
    """

    root_node = TaxTreeNode(directory)
    no_dirs = True
    for file in os.scandir(directory):
        if file.is_dir():
            # Recursively evaluate child node
            no_dirs = False
            root_node.children.append(make_tax_tree(file.path))

    root_node.is_leaf = no_dirs

    return root_node


def tax_tree_print(node: TaxTreeNode, level: int = 0):
    """Prints out a text representation of a taxonomy tree. For testing purposes."""
    print(' ' * 2 * level + node.path)
    for child in node.children:
        tax_tree_print(child, level + 1)


def get_path_list(node: TaxTreeNode) -> (list[str], list[str]):
    """
    Returns a list of file paths corresponding to leaf nodes in the tax tree that have a path from the root node such
    that each node in the path has an include value of 1, and a list of nodes for which the include value is 0.
    """
    if node.include.get() == 1:
        if node.is_leaf:
            return [node.path], []
        else:
            path_list = []
            exclude_list = []
            for child in node.children:
                p_lst, e_lst = get_path_list(child)
                path_list = path_list + p_lst
                exclude_list = exclude_list + e_lst
            return path_list, exclude_list
    else:
        return [], [node.name]
