#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 10:28:35 2017

@author: bengobengo
"""

import sys
import time
import math
import random

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



from pytree import PyTree
from matplotlib import pyplot as plt
#------------------------------------------------------------------------------

def increment_volume(tree,index):
    if index<tree.get_number_of_branches():
        tree.set_property(index,"volume",tree.get_property(index,"volume")+1)
        if tree.has_parent(index):
            parent_index=tree.get_parent_index(index)
            increment_volume(tree,parent_index)
    else:
        print('Trying to increment the volume of a non existing branch')



#------------------------------------------------------------------------------

cost=1
alpha=2

#------------------------------------------------------------------------------

trunk_init={"x":10**-6,"y":10**-6,"theta":10**-6,"volume":1,"reserve":0,
            "production":10,"leaf":1}
tree=PyTree(trunk_init)

#------------------------------------------------------------------------------

SimTime=100
for t in range(SimTime):
    print("t is"+repr(t))
    # From leaves to trunk
    for i in reversed(range(tree.get_number_of_branches())):
        # Ressources harvesting
        if tree.get_property(i,"leaf")==1:
            tree.set_property(i,"reserve",tree.get_property(i,"production"))

        # Branch maintenance
        maintenance=alpha*tree.get_property(i,"volume")
        tree.set_property(i,"reserve",tree.get_property(i,"reserve")-maintenance)
        print("here1")
        if tree.get_property(i,"reserve")<0:
            print("here2")
            tree.remove_branch(i) # A branch is removed if maintenance>reserve
            if tree.get_number_of_branches()==0:
                print('The tree died after '+repr(i)+' generations, it had '+repr(tree.get_number_of_branches())+'branches.')
                sys.exit()

        # Ressources downflow
        elif tree.has_parent(i)==1:
            parent_index=tree.get_parent_index(i)
            tree.set_property(parent_index,"reserve",tree.get_property(i,"reserve"))
            tree.set_property(i,"reserve",0)

#------------------------------------------------------------------------------

    # From trunk to leaves
    for  i in range(tree.get_number_of_branches()):

        # Creation of new branches if there are enough ressources
        excedent=int(math.floor(tree.get_property(i,"reserve")-cost))

        if excedent>0:
            nc=tree.get_number_of_children(i)
            gauss=int(round(random.gauss(2.5,0.5)))

            if nc<gauss:
                x=tree.get_property(i,"x")+math.sin(tree.get_property(i,"theta"))
                y=tree.get_property(i,"y")+math.cos(tree.get_property(i,"theta"))
                theta=random.uniform(-math.pi*0.5,0.5*math.pi)

                branch_init={"x":x,"y":y,"theta":theta,"volume":1,
                            "reserve":0,"production":1,"leaf":1}
                tree.add_branch(i,branch_init)

                if tree.get_property(i,"leaf")==1:
                    tree.set_property(i,"leaf",0)

                # Increments the volume of each ancestor of the new parent
                increment_volume(tree,i)

print(tree.get_number_of_branches())
