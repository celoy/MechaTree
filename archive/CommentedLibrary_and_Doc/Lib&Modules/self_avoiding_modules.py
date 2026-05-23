#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program simulates the growth of a ramified form with self avoiding branches.
"""

import sys
import os
import math
import random
import numpy as np

"""
Importing the PyTree class and two external modules to plot the tree and to make a
distance test to simulate the self avoiding process.
"""
from pytree import PyTree
from mod_plot import plot_tree
from mod_dist import distance_test

"""
Parameters
"""

T=380 # total simulation time

dL=0.1 # length added to growing branches at each time step
mean_angle=math.pi/3 # mean branching angle
std_angle=math.pi/10 # standard deviation for the branching angle
influenceRadius=5*dL # thickness of the influence zone


p_b=0.2 # branching probability
L_b=2*influenceRadius

"""
Coral initialisation:
  - x,y: coordinates of trunk's origin
  - theta: direction of the trunk. Angle between the trunk and vertical axis
  - L: initial length of the trunk
  - grow: if grow==1 the trunk can grow, if grow==0 it can't

  These are the properties of the trunk. Each time that a new branch is created,
  it will be initialised with the same properties but with different values.
"""

trunk={"x":0,"y":0,"theta":10**-6,"L":2,"grow":1}
coral=PyTree(trunk)

"""
Growth simulation
"""

for i in range(T):

    plot_tree(coral,i)

    """
    Traversing the tree to simulate growth of every branch: growing loop
    """
    nb=coral.get_number_of_branches()
    for n in range(nb):

        """
        If the branch can't grow we pass to the next branch
        """
        if coral.get_property(n,"grow")==1:

            """
            Avoiding branches to traverse the ground
            """
            y_0=coral.get_property(n,"y")
            L=coral.get_property(n,"L")
            theta=coral.get_property(n,"theta")

            y_tip=y_0+L*math.cos(theta)

            if y_tip<0.5:
                coral.set_property(n,"grow",0)

            """
            Realizing distance test for branch n. If it can grow, dL is added
            to its actual length, if not we pass to the next branch. If the
            branch doesn't pass the test, then it means it can't continue
            growing: the property "grow" is set to 0.
            """
            if distance_test(coral,n,influenceRadius)==1:
                L=coral.get_property(n,"L")
                coral.set_property(n,"L",L+dL)
            else:
                coral.set_property(n,"grow",0)

    #--------------------------------------------------------------------------
    """
    Traversing the tree to simulate branching for every branch: branching loop
    """
    nb=coral.get_number_of_branches()
    n=0
    while n<nb:
        """
        If branch isn't allowed to grow, then it isn't allowed to ramify, then
        we pass to the next branch
        """
        if coral.get_property(n,"grow")==1:

            """
            If random number is smaller than branching probability and if branch
            n is long enough to ramify
            """
            if random.random()<=p_b and coral.get_property(n,"L")>L_b:

                """
                Angle between the two new branches is chosen randomly
                """
                angle=random.gauss(mean_angle,std_angle)
                """
                Ramification is symmetrical
                """
                thetaplus=coral.get_property(n,"theta")+angle/2
                thetaminu=coral.get_property(n,"theta")-angle/2

                """
                The origin of the new branches has the same coordinates
                (x_0,y_0) than the parent's tip
                """
                x=coral.get_property(n,"x")
                y=coral.get_property(n,"y")
                L=coral.get_property(n,"L")
                theta=coral.get_property(n,"theta")

                x_0=x+L*math.sin(theta)
                y_0=y+L*math.cos(theta)

                """
                New branches initialisation
                """
                branch_init1={"x":x_0,"y":y_0,"theta":thetaplus,"L":0,"grow":1}
                branch_init2={"x":x_0,"y":y_0,"theta":thetaminu,"L":0,"grow":1}
                """
                New branches are added to branch n with the "add_branch"
                function
                """
                coral.add_branch(n,branch_init1)
                coral.add_branch(n,branch_init2)

                """
                The new parent can't grow anymore
                """
                coral.set_property(n,"grow",0)
                """
                Skipping the newly created branches in the loop
                """
                n=n+2
                """
                Actualizing the number of branches in the coral
                """
                nb=coral.get_number_of_branches()

        """
        Passing to the next branch
        """
        n=n+1
