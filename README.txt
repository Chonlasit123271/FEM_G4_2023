# 3D 8-Noded Hexahedral FINITE ELEMENT SOLVER

This project is part of finite element method course that aims to apply the understanding
of FEA on structure under force and displacement boundaries.(This project doesn't consider body force.)

The material consitute relationship is isotropic.


## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)



## Installation
1. Clone the repository:

git clone https://github.com/


2. Everything should be in a single folder, for the project to run correctly.

3. This preject require Matlab to compile the script.



## Usage
To use the script,



 	###[ load_apply.mlx ] 

This script is for initializing cyclic target displacement in every cycle.


---->LOAD.amp <----
Simply input an array of target amplitude in ---->LOAD.amp <----

----> step <----
To input the number of step to reach the target amplitude,
input an integer in variable ----> step <----


This will create a file call 'u_input.csv' that contains array of load increment in each step


	###[ main_main_FEM3D_8H_GROUP_4 ] 

This script is the main function that house the whole calculation.

***Change the file name to match with the mesh file FEM_COL_XX or FEM_3D_COLUMN_XXxX under
 %Calling the exported gmsh file

the node coordinate and node connectivity will be automatically read.

----> Basic Input <----


E :  Modulus N/m^2
nu : Poisson Ratio


----> support <----
this matrix will define support boundaries condition
The arguments are as follows
1 = fix 
0 = not fix
 [NODE, fix_x, fix_y, fix_z]


----> CONTROL <----

For Force control          use 'd'
For Displacement control   use 'f'
CONTROL = 'd'; is used as default due to the main problem will be displacement control


----> input_disp <---
input_disp = [ NODE uDisp_x uDisp_y uDisp_z]

for x   uDisp( #node_number# *3-2)+u_input(step)
for y   uDisp( #node_number# *3-1)+u_input(step)
for z   uDisp( #node_number# *3)+u_input(step)

for example, the input in the x direction should be

 input_disp=[ 7 uDisp( 7 *3-2)+u_input(step) 0 0;
            5  uDisp( 5 *3-2)+u_input(step) 0 0];


----> nodal_force  <----

Applying nodal force [NODE fx fy fz]
the input can be ignored if the CONTROL = 'd', it won't affect the calculation

for example

nodal_force = [1 3000/4 0 0;
            3 3000/4 0 0;
            5 3000/4 0 0;
            7 3000/4 0 0];

If everything input is done correctly, the script would start calculating stiffness matrix of each element.
Then the global stiffness matrix will be mapped.
The result solved by {u} = [K]{f} will be recorded for each load applying step.

The result of solved displacement vector and sovled force vector can be access in these two files
in the same folder:

fem_u_result.csv
fem_f_result.csv

where each column shows the results displacement/force for each load step.

in 'fem_u_result.csv' the results are shown as follows

firstvalue	-------> resultant nodal displacement vector of node 1 in [#column] step (x-direction)
secondvalue	-------> resultant nodal displacement vector of node 1 in [#column] step (y-direction)
thirdvalue	-------> resultant nodal displacement vector of node 1 in [#column] step (z-direction)
forthvalue	-------> resultant nodal displacement vector of node 2 in [#column] step (x-direction)
fifthvalue	-------> resultant nodal displacement vector of node 2 in [#column] step (y-direction)
sixthvalue	-------> resultant nodal displacement vector of node 2 in [#column] step (z-direction)
seventhvalue	-------> resultant nodal displacement vector of node 3 in [#column] step (x-direction)
.....		------->					    ...			  ...

similarly, in 'fem_v_result.csv' the results are shown as follows

firstvalue	-------> resultant nodal force vector of node 1 in [#column] step (x-direction)
secondvalue	-------> resultant nodal force vector of node 1 in [#column] step (y-direction)
thirdvalue	-------> resultant nodal force vector of node 1 in [#column] step (z-direction)
forthvalue	-------> resultant nodal force vector of node 2 in [#column] step (x-direction)
fifthvalue	-------> resultant nodal force vector of node 2 in [#column] step (y-direction)
sixthvalue	-------> resultant nodal force vector of node 2 in [#column] step (z-direction)
seventhvalue	-------> resultant nodal force vector of node 3 in [#column] step (x-direction)
.....		------->					    ...			  ...







##Authors
1) Paweeta Wintachai
2) Phoowanun Tunnikorn
3) Kyi Noo Khin
4) Suthep Phosri
5) Abzal Ilahi

##Licences



	
	










































