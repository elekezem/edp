# EDP

## Introduction
Electron density plotter is a tiny tool for generating electron density plots from VASP output.

## Compilation
EDP depends on a couple of libraries, which are normally directly available
by your favorite package manager.

* Cairo
* PCRECPP
* TCLAP

For instance, on Debian you can install these packages by:
```
sudo apt-get install libcairo2-dev libpcre3-dev libtclap-dev
```

To install the program, simply run:
```
make
```

## Usage
A short tutorial on using the program is provided in this [blog post](http://www.ivofilot.nl/posts/view/27/Visualising+the+electron+density+of+the+binding+orbitals+of+the+CO+molecule+using+VASP).
