# Project structure

## main.jl
This is the file that should be run to the the experiment. It is divided into:

- A synthetic data generation step at the start where cylinder of random transform (or identity for testing) and radiuses are generated
- Then a fictitious camera is created
- The relative quadrics, dual quadrics, singular points, planes, points are projected into conics, dual conics, etc
- A function to find a cylinder silhouette is declared and is used to find the lines defining it for each cylinder. Those are used both for visualization and for solving the problem
- A visualization step where the 3D scene is show on the left and the projected scene seen throug the camera is shown on the right
- The solver part:
	- We pick 3 lines (the minimum) from the ones that we found for the silhouettes
	- A parametrization of the rotation matrix using only linear equations in R is defined and solved for using (11)
	- Once the rotation is found the same thing is done for the translation using the original dual quadric

asserts are placed inside a begin - end block so they can be collapsed in the editor when not used

## cylinder.jl
Definitions of factory methods for quadrics creation from full parameter definitions, or randomly generatered ones

## camera.jl
Factory methods to generate camera matrices

## space.jl
Type definitions and factory methods for low level primitives related to 3D objects

## plotting.jl
Functions to plot points, lines, cylinders, conics onto the figure of Makie.jl

## fromMatrixToFormula.jl
A tool for converting from matrix or vector form to equations in x, y, z to be plotted on external online tools for debugging

## debug.jl
Functions used only during debugging like extra plotting or extra steps

## utils.jl
Definition of functions for random generation of data and notably of the ≃ "almost equal" sign used in a similar way of the ≈ approximate sign but with a precision that does not the depend of the machine epsilon but is manually choosed depending on the problem at hand

## docker-compose.yml
A docker compose for a julia environment that stays up and to which you can connect to to run and test the code

## camera.stl / camera.blend
A 3d model of a camera used for visualziation