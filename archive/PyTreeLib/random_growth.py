#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program simulates a completely random growth
"""

import sys
import os
import random
import numpy as np

"""
Adding new dependencies
"""
script_dir = os.path.dirname(__file__)

lib_dir = os.path.join(script_dir, 'Lib/')
sys.path.append(lib_dir)

mods_dir = os.path.join(script_dir, 'Mod/')
mod_dir_list = [f.path for f in os.scandir(mods_dir) if f.is_dir() ]
for mod_dir in mod_dir_list:
    sys.path.append(mod_dir)

"""
Importing the PyTree class and two external modules to plot the tree and to make a
distance test to simulate the self avoiding process.
"""

from pytree import PyTree
from mod_3Dplot import plot3D_tree

"""
Parameters
"""

NbGen=300 # number of generations
p_b = 0.2 # branching probability

maxChildren = 4

"""
Tree initialisation:
  - x,y,z: coordinates of trunk's origin
  - theta,phi: direction of the trunk

  These are the properties of the trunk. Each time that a new branch is created,
  it will be initialised with the same properties but with different values.
"""

trunk = {"x":0,"y":0,"z":0,"theta":0,"phi":0}
tree = PyTree(trunk)

"""
Growth simulation
"""

for generations in range(NbGen):

    """
    Randomly choosing were to produce a new ramification and the number of
    branches that will born from this ramification
    """

    mother = int(random.random() * (tree.get_number_of_branches() - 1))
    nbchild = int(random.random() * (maxChildren - tree.get_number_of_children(mother)))

    print("Branch "+repr(mother)+" is having "+repr(nbchild)+" children.")

    """
    Extracting the mother's properties to calculate the position of the origin
    of the children
    """
    xmom = tree.get_property(mother,"x")
    ymom = tree.get_property(mother,"y")
    zmom = tree.get_property(mother,"z")
    thetamom = tree.get_property(mother,"theta")
    phimom = tree.get_property(mother,"phi")

    """
    Calculating the coordinates of the children origin
    """
    x = xmom + np.sin(thetamom) * np.cos(phimom)
    y = ymom + np.sin(thetamom) * np.sin(phimom)
    z = zmom + np.cos(thetamom)

    if z>=1:

        for k in range(nbchild):

            """
            Randomly choosing the direction of each children
            """
            theta = np.pi * 2.0 * random.random()
            phi = np.pi * 0.5 * random.random()

            """
            Creating the properties of the new branch
            """
            properties={"x":x,"y":y,"z":z,"theta":theta,"phi":phi}

            """
            Adding the new branch
            """
            tree.add_branch(mother,properties)


    """
    Randomly removing a branch
    """
    if (random.random()>0.95 and tree.get_number_of_branches()>1):

        """
        Choosing which branch to remove
        """
        branch2rmv=int(random.random() * (tree.get_number_of_branches() - 2))+1

        print("Branch "+repr(branch2rmv)+" is going to be removed. Its last descendant is "+repr(tree.get_last_descendant_index(branch2rmv))+".")

        """
        Removing the branch
        """
        tree.remove_branch(branch2rmv)

print("Final tree size: "+repr(tree.get_number_of_branches()))

"""
Plotting the tree
"""
plot3D_tree(tree)
