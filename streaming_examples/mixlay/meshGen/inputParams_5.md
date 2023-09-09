#----------------------------------------------------------
#   Input parameters for the mesh around an ellipse
#----------------------------------------------------------
#NOTE: An input starting with 'N' is considered integer.
#NOTE: Do NOT leave an empty line in this file.
#----------------------------------------------------------
#
##Mesh dimension, 2D or 3D (extruded)
meshDim = 2D   
#
##Main Ellipse
L = 0.4
D = 1.25
sF1 = 0.12
sB1 = 0.12
##Auxiliary Ellipse (x38: its center)
L2 = 3.
D2 = 3
x38 = 0.2
sF2 = .85
sB2 = 0.85
##Layer around the ellipses
H1 = 4.0
H2 = 12.0
M1 = 1.0
M2 = 4.
M3 = 15.
##x-coor of the Centers of the Auxiliary Circles (x values)
x39 =  0
x40 = -0.5
#
#**** Number of Divisions & Grid Compressions (if <1: stretching)
##Chord of the ellipse
Nc1 = 16
Nc2 = 12
sc2 = .99
Nc3 = 12
##Normal to the ellipse (between the two ellipses)
Nx1 = 12
sx1 = 0.93  
##Normal to the bigger ellipse
Nx2 = 12
sx2 = 0.92
##Normal to the middle layer
Nx3 = 18
sx3 = .92
##After the ellipse
Nx4 = 11
sx4 = 1.08
Nx5 = 30
sx5 = 1.02
#Only for 3D (extruded) ellipse (length and num of divisions in z-direction)
Lz = 4  
Nz = 20  
#Adhoc values for a few points in the mesh (Change with Care)
adHoc14 = -14.0
adHoc15 = 4.0
adHoc27 = -25
adHoc16x =  1.4
adHoc16y = -1.5
