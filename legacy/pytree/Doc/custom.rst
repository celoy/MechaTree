===========================================
Customize the library: add external modules
===========================================

Because of its generic characteristics, the library supports the importation of external modules containing the specificities of different kinds of growth. For example, if you want to simulate the growth of a tree under mechanical constraint, , you may want that all branches have already defined properties such as its deformation or its constraints, or to dispose of an arsenal of adapted functions. In order to take advantage of these features without overcharging the main simulation programs, it may interest you to acces them from external files or modules. 

Add existing modules
====================

The library comes with some existing modules that you can add from the ``setup_mods.py`` file. To do so you have to add the extension files names to the sources command:

``sources=[mod1, mod2, ..., modN]``

where the ``modI`` are the names of the files containing the external modules. The modules are written in Cython ".pyx" files. 

Create your own modules
=======================

You may also want to create your own modules, because they might be more adapted to your interests. To do so you should write Cython ".pyx" files where the custom functions you want to add are declared. In order to make this as performant as you can, be sure you know how to take advantage of Cython features. A good start is to declare all the variables that you use in the function, just like if you were writing C code. Be sure that the modules you write are in the same directory as the ``setup_mods.py`` file, then just compile them as explained above. You should now be able to use your own extension files from Python.




