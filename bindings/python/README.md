# Welcome to OOFEM Python Bindings

OOFEM is parallel, object-oriented finite element code for solving mechanical, transport and fluid mechanics problems. 

OOFEM is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

Copyright (C) 1993 - 2021   Borek Patzak


This directory contains sources to generate OOFEM bindings to Python programming language using Pybind11.
To generate the binding, the pybind11 is required to be installed. Python module pytest is required as well.

## Prerequisites
```
pip3 install pybind11
```

Sometimes it happens that running any example raises an error:
expected constructor, destructor, or type conversion before ‘(’ token PYBIND11_MODULE(oofempy, m) {

In such a case, install pybind11 rather from source:
```
git clone https://github.com/pybind/pybind11.git; cd pybind11; pip3 install
```

For examples, please go to tests or examples subdirectories.

## Configuration
To build python bindings, one has to configure OOFEM. In particular,
set USE_PYTHON_BINDINGS to ON


## Running the tests
```
export PYTHONPATH=.
python3 -m pytest /path/to/oofem.git/bindings/python/tests
```

## Generating documentation
```
cd oofem.git/bindings/python2/docs; make html
```