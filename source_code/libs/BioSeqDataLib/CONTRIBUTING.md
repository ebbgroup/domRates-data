Contribution guidelines:
========================

Best have a look at existing code to see how you should format your code. Below are a couple of guidelines:

Formatting:
-----------

variable and method names:
- use lower camelCase: myFunction
- short but descriptive: calcIdentity 

indent style:
- a single tab for indentation

methods:
- return values on a seperate line above the function name
- otherwise: formatting as BSD/Allman style in eclipse

```C++
/**
* \brief Calculates a value.
* \param value A parameter
* \param str A second parameter
* \return The product of the two arguments.
*/
int
myFunction(int value, string str)
{
	return value * std::stoi(str);
}
```


Documentation
-------------

Please write a careful documentation using doxygen (http://www.stack.nl/~dimitri/doxygen/):
- document every function, class, file, etc. with doxygen
- add normal comments inside functions/classes where necessary to understand the code


Testing
-------

The functions in this library might be used in many different programs/projects. It is therefore important to carefully test each function:
- each method needs a test written in the boost unit testing framework
- try to think of the different scenarios that need to be handled and write a test for each case
- to run the tests: 

```bash
mkdir build
cd build
cmake -DWITH_UNIT_TEST=1 ..
make
```

- it is advisable to run valgrind to test your test for memory leaks and memory access problems
- add -DCMAKE_BUILD_TYPE=Debug to the cmake command to compile with debug options

Merge Requests
--------------

Once you finished your code and before you ask the maintainer to merge it into the master branch do the following:
- make sure that you followed all guidelines above
- merge the master branch into your development branch to solve potential conflicts
- check if really all files are included that are needed, including e.g. input files for tests

