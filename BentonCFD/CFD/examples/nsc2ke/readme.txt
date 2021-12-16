------------------------------------------------------------
NSC2KE: release 1.0 bijan mohammadi, bijan.mohamadi@inria.fr
____________________________________________________________

NSC2KE: Euler, Navier-Stokes, K-Epsilon, Two-Layer or Wall-Laws.

The fortran code is in NSC2KE.F, make "fsplit NSC2KE.F" to get files ".f". 

you should have a file 'param2d.h' containing declarations.

make make (after making the makefile compatible with your system).

The computation configuration (Mach, Reynolds,...) is in file 'DATA'.

NSC2KE is the .EXE.

"MESH" contains a mesh for a NACA 0012 with  801 points.

Try to run NSC2KE.

Visualize the results by typing 'gnuplot<GNU.DATA' (before have a look at GNU.DATA).

You see the pressure distribution CP, the friction at the wall CF.
The mesh, The pressure distribution in the flow and the Mach number
distribution and finally the L2 norm of the residual.

If you run the code for this configuration on this mesh the residual should decrease
by more than 4 order of magnitude and to give an idea it takes 2 min 46 sec on my
workstation HP 715 running at around 10 mega flops. 

---------------------------------------------------------------------------------------

You can use the mesher FreeFem to built unstructured meshes for NSC2KE. You can get FreeFEm
in the same place than NSC2KE.

Otherwise, if you have a given mesher (structured or unstructured), you can  write an easy interface
between your mesher and NSC2KE. see routine mailla.f for the format used by NSC2KE.
If your mesher gives quadrangle, cut them just by two. Just do it avoiding degenerated triangles (3 dirichlet
boundary conditions). 

For instance if your local element numbering for a quadrangle is 1,2,3,4 in the positive sens
you can have 2 triangles (1,2,3) and (1,3,4).

NSC2KE asks for 
1) number of points, number of triangles
2) x,y,and an integer giving the kind of b.c. (0 internal, 3 body, 2 symmetry, 4-5 inflow-outflow, 6 for a given profile.)
3) the connectivity table. (see mailla.f)

-------------------------------------------------------------------------------
In principle, you have received by Email the following files:

READ_ME
NSC2KE.F 
makefile 
DATA
param2d.h 
MESH 
GNU.DATA 

-------------------------------------------------------------------------------

please send bugs and comments to bijan.mohammadi@inria.fr

good luck

