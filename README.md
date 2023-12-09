The Thermal Energy Algorithms Library (TEAL)
=====

![Green-winged Teal](teal-image.jpg)

**NOTE**: This project is very early in development. 

# Getting Started 

TEAL is a MOOSE-based application. Users will need to first follow the MOOSE installation instructions
[here](https://mooseframework.inl.gov/getting_started/installation/index.html). 

## Basic MOOSE instructions 

 - Basic Instructions for the moose conda/mamba environment 
 
```
mamba init
exit
```
**NOTE**: You will need to exit a terminal to initialize conda/mamba

```
conda config --add channels https://conda.software.inl.gov/public
mamba create -n moose moose-dev
mamba activate moose
```

 - Basic Instructions for keeping moose environment updated 
 
```
mamba update --all
```

 - Basic Instructions for cloning, building, and testing MOOSE
 
You should clone MOOSE into a `~/projects` directory and place all MOOSE apps in that same directory. 

```
cd ~/projects
git clone https://github.com/idaholab/moose.git
cd moose 
git checkout master
cd ~/projects/moose/test
make -j4
./run_tests -j4
```

**NOTE**: All MOOSE tests and framework modules must be built with the MOOSE conda environment active 

## Basic TEAL instructions

Once MOOSE is installed and tested, clone this repository into the `projects` directory. 

```
cd ~/projects
git clone https://github.com/aladshaw3/teal.git
```

OR

```
cd ~/projects
git clone git@github.com:aladshaw3/teal.git
```

After you have cloned the project, go into the `teal` directory to build the library from source.

```
cd ~/projects/teal
make -j4
./run_tests -j4
```

**NOTE**: All MOOSE-based applications must be built with the MOOSE conda environment active 


# Committing Code 

This repo has a pre-commit style code formatting requirement. You cannot commit without first running:

```
git clang-format
```

# More Information on Methods and Solvers

This project is a spin-off of a project called [CATS](https://github.com/aladshaw3/cats) (Catalysis and Treatment Simulations). 
You can visit the CATS GitHub page for more detailed technical information on how to use MOOSE-based applications. 

Solver options and kernel mathematically details are available in the [CATS User Guide](https://github.com/aladshaw3/cats/blob/master/CATS-UserGuide-09172022.pdf).

# More MOOSE Information

Fork "teal" to create a new MOOSE-based application.

For more information see: [https://mooseframework.org/getting_started/new_users.html#create-an-app](https://mooseframework.org/getting_started/new_users.html#create-an-app)
