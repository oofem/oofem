************************************************************
Introduction
************************************************************

The aim of OOFEM project is to develop efficient and robust tool for
FEM computations as well as to provide modular and extensible
environment for future development.

The aim of this document is to provide the introduction to understand the OOFEM
principles and internal structure, and provides commented examples how to implement new components, such as custom elements, material models, or even new problem formulations.

The OOFEM is divided into several modules. The core module, called OOFEMlib, introduces the fundamental, core classes. These classes are problem independent
and they provide the common definitions and support to remaining modules. 
The other modules typically implement problem specific parts, such as elements, materials, solvers, etc. 
The OOFEMlib module is the compulsory module, that is always included in all builds. The other modules are optional and can be included or excluded from the build, depending on the user configuration.

The OOFEM package is written in C++. This 
document contains many examples and listings of a parts of the source
code. Therefore, the ideal reader should be familiar with C++ programming
language. However, any reader with an object-oriented background should
follow this document, since the examples are written in a form, which is hopefully
easy to read and understand.
