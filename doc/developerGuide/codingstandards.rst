Coding Standards
================

Naming Conventions
------------------

The names of classes, attributes, services, variables, and functions
in a program serve as comments of a sort. So don't choose terse names instead, look for names that give useful information about the meaning of the variable or function. In a OOFEM
program, names should be English, like other comments. 

Local variable names can be shorter, because they are used only within one context, where (presumably) comments explain
their purpose. 

Try to limit your use of abbreviations in symbol names. It is ok to make a few abbreviations, explain what they mean, and
then use them frequently, but don't use lots of obscure abbreviations. 

Please use capital letters to separate words in a name. Stick to
lower case; reserve upper case for macros, and for name-prefixes that
follow a uniform convention. The function or service name should always begin with lovercase
letter, the first uppercase letter in function name indicates, that
function is returning newly allocated pointer, which has to be dealocated.
For example, you should use names like ignoreSpaceChangeFlag; 

When you want to define names with constant integer values, use enum rather than `#define`. GDB knows about
enumeration constants. 

Use descriptive file names. The class declarations should be placed in \*.h files
and class implementation in corresponding \*.C files. For each class, create a separate files.


Parameters and Return Values
-----------------------------

The prefered argunent passing method for objects is by reference. 
Try to avoid rerturning pointers to arrays or matrices (or generally
to any component), since it is not clear, whether to delocate the
returned pointer or not. The most prefered way is to create local
variable of vector or matrix type, pass it (using reference) to called
function. The calling function is responsible to properly resize the
(output) parameter and set values accordingly. The point is, that
destructors are called for local variables automatically by compiler,
so there is no possibility for memory leaks and the local variable can
be reused for multiple calls to target function (inside loop) and
therefore there is no need for repeating memory allocation and
dealocation. Sometimes it may be reasonable to return pointer to
constant float array or matrix (representing for example
nodal coordinates), since passing output array as parametr will require 
array copying. 


	