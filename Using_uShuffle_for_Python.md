## How to use uShuffle in Python

Follow instructions found here for writing a wrapper module: 

http://digital.cs.usu.edu/~mjiang/ushuffle/python.html

For building the module, instead of their commands enter:

```
gcc -c -I/usr/include/python2.6 -fPIC ushuffle.c ushufflemodule.c
gcc -shared ushuffle.o ushufflemodule.o -o ushuffle.so

```
