# MES

Console app that uses FEM to simulate the process of heating up an infinitly long bar with rectangular cross section.

Console gives you a temperature in every element after every step but the app also saves the min and max temperatures in an iteration in a csv file.

To use it clone it and play with the variables in the main function:
- H - hight of the cross section
- B - width of the cross section
- nH - number of nodes on hight
- nB - number of nodes on width
- numberOfPoints - is a number of integration points in Gaussian quadrature
- tempOS - temperature of surroundings
- temp0 - starting temperature
- dt - time step
- tend - time of simulation
- ro - density of material
- cp - specific heat of material
- k - heat conductivity of material
- alph - heat transfer coefficient
