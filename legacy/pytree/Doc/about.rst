==================
About this Library
==================

PyTreeLib aims to be an useful tool for the simulation of theoretical models describing the growth of branching structures. The core of the library is written in C++ but the user interacts with the core using Python. The interface between C++ and Python was done using Cython. Our intention is to provide a plateform where even the the less familiar with the numeric tools, are able to simulate the growth of ramified structures in a simple and straightforward way, without loosing the characteristic performance of low-level programming languages.  

One of purposes of PyTreeLib is to be generic enough to allow the user to simulate any kind of growth process involving ramified structures, such as trees, corals or neurons for example. This is achieved by letting the user decide which are the properties of the structure which interests him, and how they evolve. This plateforms allows the user to interact in a straightforward way with a structure (the tree) consisting on a set of elements (the branches) linked by family relations and with user-defined characteristics.

If you are interested in contributing to the developpement of this library we advice you to read the whole documentation. If you are just interested in using the library and you don't need or want to understand its functioning in deep, you can stop your reading after "Users guide" section. 
