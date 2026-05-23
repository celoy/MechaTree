#cython: language_level=3
#cython: wraparound=False
#cython: boundscheck=False
#cython: cdivision=True

import sys
import os
import cython
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""
Function declaration: plot3D_tree(PyTree tree)
- This function plots a 3D-tree. It adds green leaves at the extremities of the
branches who have no children.
"""

def plot3D_tree(tree):

    """
    Variable declaration
    """
    cdef int i
    cdef str num
    cdef double x1
    cdef double y1
    cdef double z1
    cdef double x2
    cdef double y2
    cdef double z2
    cdef double xleaf
    cdef double yleaf
    cdef double zleaf
    cdef double leaf_length=0.2

    """
    If you want to add the axis on the figure set this option to 'on'
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.axis('off')

    """
    We will traverse the tree plotting each branch at a time
    """
    i=0
    for i in range(tree.get_number_of_branches()):

        """
        These are the coordinates of the branch's origin
        """
        x1 = tree.get_property(i,"x")
        y1 = tree.get_property(i,"y")
        z1 = tree.get_property(i,"z")

        """
        These are the coordinates of branch's tip
        """
        x2 = x1 + np.sin(tree.get_property(i,"theta")) * np.cos(tree.get_property(i,"phi"))
        y2 = y1 + np.sin(tree.get_property(i,"theta")) * np.sin(tree.get_property(i,"phi"))
        z2 = z1 + np.cos(tree.get_property(i,"theta"))

        """
        The figure is created and a straight line representing the branch is plotted
        in red
        """
        ax.plot([x1,x2],[y1,y2],[z1,z2],'k',linewidth=2.5)

        if tree.get_number_of_children(i)<1:
            xleaf = x2 + leaf_length * np.sin(tree.get_property(i,"theta")) * np.cos(tree.get_property(i,"phi"))
            yleaf = y2 + leaf_length * np.sin(tree.get_property(i,"theta")) * np.sin(tree.get_property(i,"phi"))
            zleaf = z2 + leaf_length * np.cos(tree.get_property(i,"theta"))

            ax.plot([x2,xleaf],[y2,yleaf],[z2,zleaf],'g',linewidth=5)

    plt.show()
