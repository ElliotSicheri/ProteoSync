"""Contains functions for creating, manipulating, and interacting with database taxonomy trees."""

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
    """Creates a taxonomy tree from the file structure of the given directory. Returns the root node.
    Internal node directories should only contain directories. Leaf node directories should only contain files.
    It is assumed that a valid and properly structured directory pathway is passed to the function."""

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


def get_path_list(node: TaxTreeNode) -> list[str]:
    """Returns a list of file paths for leaf nodes in the tax tree that have a path from the root node such that each
    node in the path has an include value of 1."""
    if node.include.get() == 1:
        if node.is_leaf:
            return [node.path]
        else:
            path_list = []
            for child in node.children:
                path_list = path_list + get_path_list(child)
            return path_list
    else:
        return []
