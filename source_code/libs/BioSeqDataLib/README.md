BioSeqDataLib
=============

Library to deal with biological data. Main focus is currently on sequence data,
with special focus on protein domains. This library will be the core library
for the different domain related programs developed in the group.

Requirements
------------

We try to keep the dependencies as little as possible. Current dependencies are:
- cmake (https://cmake.org)
- compiler supporting c++11 - e.g. g++ 4.8 or higher
- boost modules: system filesystem (http://www.boost.org/)

Optional:
- boost modules: unit_test - only needed when unit_tests should be compiled

In most Linux distributions (e.g. Ubuntu, Arch Linux) it should be possible to install these dependencies using the package manager.

Installation
------------

Change into the BioSeqDataLib and run the following commands:

```bash
mkdir build
cd build
cmake ..
make
```


For install run:

```bash
sudo make install
```

Documentation
------------

Each function is documented using Doxygen. To make the documentation change into the doc folder and type:

```bash
doxygen Doxyfile
```

Problems, Bugs & Suggestions
----------------------------

We try our best not to have any bugs in the code, unfortunately some will probably avoid detection. If you encounter one or if you have suggestions or questions, please write an email to domainworld[ at ]uni-muenster.de.

Credits
-------

Code of AlgorithmPack has been used.
