================
Getting started
================


Overview of the library
=======================

PyTreeLib consists on two C++ classes: the **Branch** class and the **Tree** class. Each instance of the Branch class is charatacterized by a set of family relations with other instances of the class (pointers to parent and children) and a set of user-defined properties of type **double**. An instance of the Tree class is a container (a C++ **vector**) of pointers to Branch objects. 

The Tree class is wrapped by the **PyTree** class using Cython. This means that PyTree is a "Python version" of the C++ class Tree. Many of the methods of the Tree class have an equivalent in the PyTree class, which call the C++ methods when called. 

The user can write a Python script where PyTree objects are initialized and manipulated via the set of wrapped methods of the PyTree class. These methods, which can be called from Python, are members of the PyTree class and act by calling the equivalent methods of the C++ Tree class. Hence, using a Python script we interact with the C++ core of the ibrary. 

Installing and using PyTreeLib
==============================

In order to use PyTreeLib you will need Python3 and a C++ compiler (gcc, g++). We advice you to install **anaconda** to have Python3 and all of its packages installed at once. If you choose to install **anaconda** you must install the package ``libgcc`` with the following command in the terminal:

``conda install libgcc``.

If you don't, you will probably have the following compilation error:

``Import error undefined symbol (C++ module in python) ZTINSt8ios_base7failureB5cxx11E``.

Before building the library, be sure that all the files of the library are in the same directory. Then, set your working directory to the directory containing the library and enter the following command at the terminal:

``python3 setup.py build_ext --inplace``.

You have built the library and now are able to use it! To do this you must write the following line in your Python code: 

``from pytree import PyTree``.

Be sure to write your Python code in the directory containing the files of the library.


Instruction to download.



