#cython: language_level=3
#cython: wraparound=False
#cython: boundscheck=False
#cython: cdivision=True

import sys
import os
import cython
import math
from matplotlib import pyplot as plt

"""
Function declaration: plot_tree(PyTree tree, int iteration)
- This function plots the tree at the iteration 'iteration' of the growth simulation.
The figures are all stocked in the directory 'working_directory/selfAvoid'. If
'selfAvoid' doesn't exist the function creates it. The result is a directory
containing snapshots of the ramified structure at every stage of its growths which can
be used to make a video.
"""

def plot_tree(tree, int iteration):
    """
    Variable declaration
    """
    cdef int i
    cdef str num
    cdef double x1
    cdef double z1
    cdef double x2
    cdef double z2

    """
    If you want to add the axis on the figure set this option to 'on'
    """
    plt.axis('off')

    """
    We will traverse the tree plotting each branch at a time
    """
    i=0
    for i in range(tree.get_number_of_branches()):

        """
        These are the coordinates of the branch's origin
        """
        x1=tree.get_property(i,"x")
        z1=tree.get_property(i,"y")

        """
        These are the coordinates of branch's tip
        """
        x2=x1+(tree.get_property(i,"L")*math.sin(tree.get_property(i,"theta")))
        z2=z1+(tree.get_property(i,"L")*math.cos(tree.get_property(i,"theta")))

        """
        The figure is created and a straight line representing the branch is plotted
        in red
        """
        plt.figure(1)
        plt.plot([x1,x2],[z1,z2],'red',linewidth=1)

    """
    Instructions to save the figure in the directory: 'working_directory/selfAvoid'.
    This function allows you to create as much as 9999 figures. Figure names are:
      'fig0001.png', 'fig0002.png', ... , 'fig9999.png'.
    """
    working_dir = os.getcwd()
    results_dir = os.path.join(working_dir, 'SimResults/selfAvoidImages/')

    """
    If the directory 'working_directory/selfAvoid' doesn't exist, this function
    creates it
    """
    if os.path.exists(results_dir) == False:
        os.makedirs(results_dir)

    if iteration<10:
        num=str(0)*3+str(iteration)
    elif iteration>=10 and iteration<100:
        num=str(0)*2+str(iteration)
    elif iteration>=100 and iteration<1000:
        num=str(0)*1+str(iteration)
    else:
        num=str(iteration)

    title="fig"
    form=".png"
    figname=title+num+form
    plt.savefig(results_dir+figname)
    plt.close()
