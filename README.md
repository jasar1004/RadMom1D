# RadMom1D

SOFTWARE
--------

Welcome to RadMom1D.  RadMom1D  is a computational  framework  for
solving the radiative transfer equation in  one spatial dimension.
It includes a  package of C++ subroutines for  solving radiative
transfer problems  using  numerical methods (i.e., Computational  Fluid
Dynamics  or  CFD). In particular, a Godunov-type finite volume method with
piecewise linear reconstruction of the solution within the computational
domain is used to solve the system of partial differential equations of
interest. All required files are in  this directory and  the corresponding
subdirectories. In this framework, systems of PDEs, resulting from the
application of the method of moments to the radiative transfer equation, are
solved. The current implementation involves three moment closure techniques:

- M1: the first-order maximum entropy moment closure, which involves solving
      angular moments up to first order. The unclosed second-order moment in
      the transport equation for the first-order moment is expressed in terms
      of the lower-order moments by assuming a distribution that maximizes the
      radiative entropy and reproduces the lower-order moments;

- P1: the first-order spherical harmonic moment closure, which involves solving
      angular moments up to first order. The unclosed second-order moment in
      the transport equation for the first-order moment is expressed in terms
      of the lower-order moments by expanding the distribution as a series expansion
      in terms of orthogonal legendre polynomials (1D equivalent of 3D spherical
      harmonics) with the constraint that the lower-order moments are angular moments
      of the assumed distribution;

- P3: the third-order spherical harmonic moment closure, which involves solving
      angular moments up to third order. The unclosed fourth-order moment in
      the transport equation for the third-order moment is expressed in terms
      of the lower-order moments by expanding the distribution as a series expansion
      in terms of orthogonal legendre polynomials (1D equivalent of 3D spherical
      harmonics) with the constraint that the lower-order moments are angular moments
      of the assumed distribution.

The RadMom1D code is templated such that it can be used to solve the systems of PDEs
arising from moment closure techniques of varying order, which involve different number
of independent variables. As an example, the M1 and P1 closures involve only two independent
variables in 1D, whereas the P3 closure involves 4 independent variables in one spatial
dimension. The type of moment closure technique to be used must be specified prior to
compiling the RadMom1D code so as to ensure the proper specialization of the templated
classes are compiled.

---------
COPYRIGHT
---------

Copyright (C) Joachim Sarr.


-------
LICENSE
-------

----------
DISCLAIMER
----------


------------
CONTRIBUTORS
------------

Contributors to the development of RadMom1D include:
Joachim Sarr

----------------------------------------------------------------------- 

------------------------------
Classes:
------------------------------

RadMom1D_Input_Parameters: This class includes definition and manipulation of 1D RadMom input variables.

RadMom1D_UniformMesh: contains definitions for RadMom1D computational cell (solution and geometry).

RadMom1D_pState_First_Order: contains primitive variable solution state class definition for first-order moment closures (i.e, M1 and P1).

RadMom1D_cState_First_Order: contains conserved variable solution state class definition for first-order moment closures (i.e, M1 and P1).

RadMom1D_pState_Third_Order: contains primitive variable solution state class definition for third-order moment closures (i.e, P3).

RadMom1D_cState_Third_Order: contains conserved variable solution state class definition for third-order moment closures (i.e, P3).

Medium1D_State: contains gas state class definition for an absorbing, scattering participating medium. In particular, it stores radiative properties of the participating medium such as the absorption coefficient and the scattering coefficient. For now, the radiative properties implemented are based on the assumption of gray media (not spectrally dependent radiative properties) and of isotropic scattering. More complex/realistic gas radiative property models will be implemented in the future.

Cell1D_Uniform: contains class definitions for the 1D uniform cell making up the computational domain.

------------------------------
How to run the RadMom1D code:
------------------------------

------
Step 1:
------
In the RadMom1D.cc file,  set the variable "Closure_Type_ptr" to the desired moment closure technique to be used for providing approximate solution to the radiative transfer equation.
The possible options are:
 - Closure_Type_ptr = "M1": for the M1 moment closure
 - Closure_Type_ptr = "P1": for the P1 moment closure
 - Closure_Type_ptr = "P3": for the P3 moment closure

-------
Step 2:
-------
In the command line, execute the command "make clean" to clean up all the object files generated from previous compilations of the RadMom1D code.

-------
Step 3:
-------
In the command line, execute the command "make radMom1D" to compile the RadMom1D code. This will generate an executable named "radMom1D".

-------
Step 4:
-------
Once the code compiled, it can be used along with a given set of inputs to solve the resulting set of PDEs and write the solution to an output file (in *.dat format with can be visualized using Tecplot). The command used to run a test case with input file name "test.in" is the following:

./radmom1D -f test1D.in


-------
Tests:
-------
In the "Test" directory, some test samples are provided for use with the radMom1D exectutable.
These tests involves soving the radiative transfer equation using both the M1, P1, and P3 moment closures (see subdirectories "M1", "P1", and "P3" in "Test" directory).
The input file name associated with each of these test is "test1D.in". The scripts for executing the RadMom1D code can be found in these subdirectories, i.e., M1_Closure_script.sh for M1, P1_Closure_script.sh for P1, and P3_Closure_script.sh for P3 and can be exectuted in the command line using:
./M1_Closure_script.sh, or ./P1_Closure_script.sh, or ./P3_Closure_script.sh
Upon completion of the numerical simulations, an output file "test1D.dat" will be generated.

A Tecplot macro is provided in the "Test" subdirectory which extracts solution profiles for the M1, P1, and P3 closure (Comparisons_E.mcr for the radiative energy density or zeroth-order moment, Comparisons_F.mcr for the radiative heat flux or first-order moment, and Comparisons_Sr.mcr for the radiative source term) and compares them on the same figure.

