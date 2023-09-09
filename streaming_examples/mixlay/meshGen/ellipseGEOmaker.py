################################################################
#   Script to create ellipse.geo from input arg inputParams.md
#      which is used by gmsh to generate structured mesh
#      around a 2D or 3D (extruded) ellipse.
################################################################
# Saleh Rezaeiravesh, salehr@kth.se
#-------------------------------------------------------
#
import sys
import math as mt
import numpy as np

#/////////////////////////////
def inParams_reader(fileName):
    """
        Read input parameters and make a dictionary based on them
    """
    F=open(fileName,'r')
    #Define list of quantities whose values are to be read from the file
    #NOTE: If a key starts with 'N', it is assumed to be an integer
    qList=['meshDim','L','D','sF1','sB1','L2','D2','x38','sF2','sB2','H1',
           'H2','M1','M2','M3','x39','x40','Nc1','Nc2','sc2','Nc3','Nx1','sx1',
           'Nx2','sx2','Nx3','sx3','Nx4','sx4','Nx5','sx5','Lz','Nz',
           'adHoc14','adHoc15','adHoc27','adHoc16x','adHoc16y']
    #read in the lines
    ain=F.readlines()
    ain_sep=[]
    for i in range(len(ain)):
        ain_sep.append(ain[i].split())
    #search in the file for the values associted to the qList
    params={}
    for q_ in qList:
        for i in range(len(ain_sep)):
            key_=ain_sep[i][0]
            if q_==key_:               
               val_=ain_sep[i] [-1]
               if key_!='meshDim':
                  if key_[0]=='N':  #if starts with 'N'
                     val_=int(val_) 
                  else:
                     val_=float(val_) 
               params.update({key_:val_}) 
    return params               
#
#////////////////////////
def set_inputs(fileName):
    """
       Read input parameters, construct a dictionary from them, and convert the dictionay keys to the code variables.
    """
    params=inParams_reader(fileName)
    #convert the dictionary entries to code variables
    for par in params.keys():
        globals()[par]=params[par]
#
#/////////////////////////
def yEllipse(x,xc,yc,a,b):
    """  
        Returns y associated to x, such that the point (x,y) is located on an ellipse with 
        center (xc,yc) and width a and height b. 
    """
    a_=a/2.
    b_=b/2.
    return(yc+b_*np.sqrt(1.-(x-xc)**2./a_**2.))
#    
#/////////////////
def geo_header(F):
    """
       Write the header for the .geo file
    """
    F.write('/*\n')
    F.write('***** gmsh Script for generating 2D/3D mesh for flow over an ellipse *****\n')
    F.write('\tFor nomenclature, see the attached figure.\n')
    F.write('\tChoose either 2D or 3D mesh in GENERAL SETTING\n')
    F.write('\tSet the mesh parameters in SETTINGS\n')
    F.write('\tTo generate 3D mesh: gmsh ellipse.geo -3 -order 2\n')
    F.write('\t>>> To generate 2D mesh: gmsh ellipse.geo -2 -order 2\n')
    F.write('\tSaleh Rezaeiravesh, salehr@kth.se *****\n')
    F.write('*/\n\n\n')
#
#/////////////////
def geo_points(F):
    """
       Define the coordinates of the points constructing the geometery
    """
    #----- ASSIGNMent
    idLastMainPts=40   #The ID of the last main points in the mesh block. From this point on we have auxiliary points on the ellipses
    #----------------------------
    def point(i_,x_,y_,z_,F):
        F.write("Point(%d)={%g,%g,%g,1.0};\n" %(i_,x_,y_,z_))
        return {'x':x_,'y':y_,'z':z_}
    def auxEllipsePts(pStrt,pEnd,idStrt,idEnd,lineID,pCenter,a,b,iSgn):
        """ 
           Create nAuxEll auxiliary points between pStrt & pEnd so that all
                  points reside on an ellipse with center pCenter, length a, and width b.
           lineID is the ID of the line between pStrt and pEnd
           iSgn: +/-, sign of y for a given x s.t. the point (x,y) resides on the ellipse
                   
        """
        dx=(pEnd['x']-pStrt['x'])/(nAuxEll+1)
        x=np.linspace(pStrt['x']+dx,pEnd['x']-dx,nAuxEll)
        y=iSgn*yEllipse(x,pCenter['x'],pCenter['y'],a,b)
        auxIDStrt=(lineID-1)*nAuxEll+idLastMainPts
        auxID=[]
        for i in range(nAuxEll):
            auxID_=auxIDStrt+i+1
            auxID.append(auxID_) 
            point(auxID_,x[i],y[i],0.0,F)
        aux_={'mainPt.start':str(idStrt),'mainPt.end':str(idEnd),'nAuxEll':nAuxEll,'auxPts.x':x,'auxPts.y':y,'auxID':auxID}
        return aux_

    #>>> Main Ellipse
    # Main points
    P1 =point(1,0,0,0,F)     #Ellipse Leading edge=Coordinates origin
    P4 =point(4,L,0,0,F)     #Ellipse Trailing edge 
    P37=point(37,L/2.,0,0,F) #Center of main ellipse (body)
    P2 =point(2,sF1, yEllipse(sF1,P37['x'],P37['y'],L,D),0.,F)
    P6 =point(6,sF1,-yEllipse(sF1,P37['x'],P37['y'],L,D),0.,F)
    P3 =point(3,L-sB1, yEllipse(L-sB1,P37['x'],P37['y'],L,D),0.,F)
    P5 =point(5,L-sB1,-yEllipse(L-sB1,P37['x'],P37['y'],L,D),0.,F)
    # Auxiliary points (between main points to define a spline)
    auxList=[]
    auxList.append(auxEllipsePts(P1,P2,1,2,1,P37,L,D,1))
    auxList.append(auxEllipsePts(P2,P3,2,3,2,P37,L,D,1))
    auxList.append(auxEllipsePts(P3,P4,3,4,3,P37,L,D,1))
    auxList.append(auxEllipsePts(P4,P5,4,5,4,P37,L,D,-1))
    auxList.append(auxEllipsePts(P5,P6,5,6,5,P37,L,D,-1))
    auxList.append(auxEllipsePts(P6,P1,6,1,6,P37,L,D,-1))
    #>>> Bigger Ellipse
    # Main points
    P38=point(38,x38,0,0,F)    #Center of the bigger ellipse
    P7 =point(7 ,P38['x']-L2/2,0.0,0.0,F)   #Leading edge of the bigger ellipse
    xTmp=P7['x']+sF2
    P8 =point(8 ,xTmp, yEllipse(xTmp,P38['x'],P38['y'],L2,D2),0.0,F)
    P12=point(12,xTmp,-yEllipse(xTmp,P38['x'],P38['y'],L2,D2),0.0,F)
    P10=point(10,P7['x']+L2,0.0,0.0,F)      #trailing edge of the bigger ellipse
    xTmp=P10['x']-sB2
    P9 =point(9 ,xTmp, yEllipse(xTmp,P38['x'],P38['y'],L2,D2),0.0,F)
    P11=point(11,xTmp,-yEllipse(xTmp,P38['x'],P38['y'],L2,D2),0.0,F)
    # Auxiliary points (between main points to define a spline)
    auxList.append(auxEllipsePts(P7,P8 ,7,8 ,7,P38,L2,D2,1))
    auxList.append(auxEllipsePts(P8,P9 ,8,9 ,8,P38,L2,D2,1))
    auxList.append(auxEllipsePts(P9,P10,9,10,9,P38,L2,D2,1))
    auxList.append(auxEllipsePts(P10,P11 ,10,11 ,10,P38,L2,D2,-1))
    auxList.append(auxEllipsePts(P11,P12 ,11,12 ,11,P38,L2,D2,-1))
    auxList.append(auxEllipsePts(P12,P7 ,12,7 ,12,P38,L2,D2,-1))
    #>>> Layer around the ellipses
    P39=point(39,x39,0.0,0.0,F)    #center of the smaller circlular arc
    #upper
    adHoc14_=adHoc14*sF1  #adhoc
    adHoc15_=adHoc15*sF1  #adhoc
    P14=point(14,P8['x']+adHoc14_ ,H1,0.0,F)  
    P15=point(15,P9['x']+adHoc15_ ,H1,0.0,F)  
    R1=mt.sqrt((P14['x']-P39['x'])**2.+(P14['y']-P39['y'])**2.)
    P13=point(13,P39['x']-R1,P39['y'],0.0,F)    
    P16=point(16,P15['x']+M1+adHoc16x,P14['y']+adHoc16y,0.0,F)
    P17=point(17,P16['x']+M2-adHoc16x,P16['y'],0.0,F)
    P18=point(18,P17['x']+M3,P16['y'],0.0,F)
    #lower
    P19=point(19,P14['x'],-P14['y'],0.0,F)
    P20=point(20,P15['x'],-P15['y'],0.0,F)
    P21=point(21,P16['x'],-P16['y'],0.0,F)
    P22=point(22,P17['x'],-P17['y'],0.0,F)
    P23=point(23,P18['x'],-P18['y'],0.0,F)
    #center
    P24=point(24,P17['x']-adHoc16x,P10['y'],0.0,F)   
    P25=point(25,P18['x'],P10['y'],0.0,F)   

    #>>> The outermost layer of the mesh
    P40=point(40,x40,0.0,0.0,F)    #center of the larger circlular arc
    #upper
    adHoc27_=adHoc27*sF1   #adhoc
    P27=point(27,P14['x']+adHoc27_,P14['y']+H2,0.0,F)   
    R2=mt.sqrt((P27['x']-P40['x'])**2.+(P27['y']-P40['y'])**2.)
    P26=point(26,P40['x']-R2,P40['y'],0.0,F)    
    P28=point(28,P15['x'],P15['y']+H2,0.0,F)   
    #P29=point(29,P28['x']+M1,P15['y']+H2,0.0,F)   
    P29=point(29,P16['x'],P15['y']+H2,0.0,F)   
    P30=point(30,P17['x'],P15['y']+H2,0.0,F)   
    P31=point(31,P18['x'],P15['y']+H2,0.0,F)   
    #lower
    P32=point(32,P27['x'],-P27['y'],0.0,F)
    P33=point(33,P28['x'],-P28['y'],0.0,F)
    P34=point(34,P29['x'],-P29['y'],0.0,F)
    P35=point(35,P30['x'],-P30['y'],0.0,F)
    P36=point(36,P31['x'],-P31['y'],0.0,F)
    return auxList
#
#////////////////////////////////
def geo_ellipseLines(F,auxList):
    """
       Define lines of the main and bigger ellipses
    """
    F.write("\n\n// Lines\n")
    #Define auxiliary points between main points on the main ellipse and then fit a spline
    for iLn in range(len(auxList)):   #line number
        ptIDstart=auxList[iLn]['mainPt.start']
       
        nAuxEll  =auxList[iLn]['nAuxEll']
        ptIDend  =auxList[iLn]['mainPt.end']
        F.write("Line(%d) = {%d, " %(iLn+1,int(ptIDstart))) 
        for j in range(nAuxEll):
            F.write("%d, " %auxList[iLn]['auxID'][j])
        F.write("%d};\n" %(int(ptIDend))) 
    #lines between the two ellipses
    F.write("Line(13) = {7, 1};\n")
    F.write("Line(14) = {8, 2};\n")
    F.write("Line(15) = {9, 3};\n")
    F.write("Line(16) = {10, 4};\n")
    F.write("Line(17) = {11, 5};\n")
    F.write("Line(18) = {12, 6};\n")
#
#/////////////////////
def geo_otherLines(F):
    """
       Define lines of all the blocks both the two ellipses
    """        
    def line(iLn,iPtStart,iPtEnd,F):
        F.write("Line(%d) = {%d, %d};\n" %(iLn,iPtStart,iPtEnd))
    def circle(iLn,iPtStart,iCenter,iPtEnd,F):
        F.write("Circle(%d) = {%d, %d, %d};\n" %(iLn,iPtStart,iCenter,iPtEnd))
    
    circle(19,13,39,14,F)
    line(20,14,15,F)
    line(21,15,16,F)
    line(22,16,17,F)
    line(23,17,18,F)
    circle(24,13,39,19,F)
    line(25,19,20,F)
    line(26,20,21,F)
    line(27,21,22,F)
    line(28,22,23,F)
    line(40,10,24,F)
    line(41,24,25,F)
    #radial lines
    line(29,13,7,F)
    line(30,14,8,F)
    line(31,15,9,F)
    line(32,16,10,F)
    line(33,17,24,F)
    line(34,18,25,F)
    line(35,19,12,F)
    line(36,20,11,F)
    line(37,21,10,F)
    line(38,22,24,F)
    line(39,23,25,F)
    #outermost layer of the mesh
    circle(42,26,40,27,F)
    line(43,27,28,F)
    line(44,28,29,F)
    line(45,29,30,F)
    line(46,30,31,F)
    circle(47,26,40,32,F)
    line(48,32,33,F)
    line(49,33,34,F)
    line(50,34,35,F)
    line(51,35,36,F)
    line(52,26,13,F)
    line(53,27,14,F)
    line(54,28,15,F)
    line(55,29,16,F)
    line(56,30,17,F)
    line(57,31,18,F)
    line(58,32,19,F)
    line(59,33,20,F)
    line(60,34,21,F)
    line(61,35,22,F)
    line(62,36,23,F)
#
#///////////////////
def geo_surfaces(F):
    """
       Define block surfaces 
    """
    def surf(id_,i1,i2,i3,i4,F):
        F.write("Line Loop(%d) = {%d,%d,%d,%d}; Plane Surface(%d) = {%d};\n" %(id_,i1,i2,i3,i4,id_,id_))
    F.write("\n//Surfaces\n")
    surf(1,7,14,-1,-13,F)
    surf(2,8,15,-2,-14,F)
    surf(3,9,16,-3,-15,F)
    surf(4,10,17,-4,-16,F)
    surf(5,11,18,-5,-17,F)
    surf(6,12,13,-6,-18,F)
    surf(7,19,30,-7,-29,F)
    surf(8,20,31,-8,-30,F)
    surf(9,21,32,-9,-31,F)
    surf(10,22,33,-40,-32,F)
    surf(11,23,34,-41,-33,F)
    surf(12,29,-12,-35,-24,F)
    surf(13,-11,-36,-25,35,F)
    surf(14,-10,-37,-26,36,F)
    surf(15,40,-38,-27,37,F)
    surf(16,41,-39,-28,38,F)
    surf(17,42,53,-19,-52,F)
    surf(18,43,54,-20,-53,F)
    surf(19,44,55,-21,-54,F)
    surf(20,45,56,-22,-55,F)
    surf(21,46,57,-23,-56,F)
    surf(22,52,24,-58,-47,F)
    surf(23,25,-59,-48,58,F)
    surf(24,26,-60,-49,59,F)
    surf(25,27,-61,-50,60,F)
    surf(26,28,-62,-51,61,F)
#
#////////////////
def geo_trans(F):
    """
       Construct the 1D mesh
    """
    F.write("\n// Assigne divisions of line segments \n")
    def trans(lineIDlist,n):
        F.write("Transfinite Line {")
        for i in range(len(lineIDlist)):
            id_=lineIDlist[i]
            tmp_=', '
            if i==len(lineIDlist)-1:
               tmp_=' '
            F.write(("%d"+tmp_) %id_)
        F.write("} = %d;\n" %n)
    def transPB(lineIDlist,n,v_,task_):
        """ task_='Progression' or 'Bump' """
        F.write("Transfinite Line {")
        for i in range(len(lineIDlist)):
            id_=lineIDlist[i]
            tmp_=', '
            if i==len(lineIDlist)-1:
               tmp_=' '
            F.write(("%d"+tmp_) %id_)
        F.write("} = %d Using %s %g;\n" %(n,task_,v_))
    #ellipse chordwise       
    trans([1,7,19,42],Nc1)
#    trans([2,8,20,43],Nc2)
    transPB([2,8,20,43],Nc2,sc2,'Bump')
    trans([3,9,21,44],Nc3)
    trans([4,10,26,49],Nc3)
    transPB([5,11,25,48],Nc2,sc2,'Bump')
    trans([6,12,24,47],Nc1)
    #normal to the ellipse
    transPB([13,14,15,16,17,18],Nx1,sx1,'Progression')
    #layer above the bigger ellipse
    transPB([29,30,31,32,33,34,35,36,37,38,39],Nx2,sx2,'Progression')
    transPB([50,27,40,22,45],Nx4,sx4,'Progression')
    transPB([51,28,41,23,46]  ,Nx5,sx5,'Progression')
    #the outermost layer of the mesh
    transPB([52,53,54,55,56,57,58,59,60,61,62],Nx3,sx3,'Progression')
#
#////////////////
def geo_setBC(F):
    F.write("\n// Set Boundary conditions\n")        
    if meshDim=='2D':
       F.write("Physical Line(\"inlet\") = {42,47};\n")
       F.write("Physical Line(\"outlet\") = {57,34,39,62};\n")
       F.write("Physical Line(\"wall\") = {1,2,3,4,5,6};\n")
       F.write("Physical Line(\"freestreamUp\") = {43,44,45,46};\n")
       F.write("Physical Line(\"freestreamLo\") = {48,49,50,51};\n")
       F.write("Physical Surface(\"wholeDomain\")={1:26};\n\n")
       F.write('Recombine Surface "*";\n')
       F.write('Transfinite Surface "*";\n')
    elif meshDim=='3D':
       F.write('Recombine Surface "*";\n')
       F.write('Transfinite Surface "*";\n\n')
       F.write("//make a 3d mesh by extrusion in z-direction\n")       
       F.write("mesh3D[]=Extrude {0,0,%g}\n" %Lz )
       F.write("{\n")
       F.write("\tSurface{1:26};\n")
       F.write("\tLayers{%d};\n" %Nz) 
       F.write("\tRecombine;\n")
       F.write("};\n")
       F.write("Physical Surface(\"inlet\") = {423,545};\n")
       F.write("Physical Surface(\"outlet\") = {515, 295, 405, 625};\n")
       F.write("Physical Surface(\"wall\") = {167, 145, 123, 189, 79, 101};\n")
       F.write("Physical Surface(\"freestreamUp\") = {445, 467, 489, 511};\n")
       F.write("Physical Surface(\"freestreamLo\") = {563, 585, 607, 629};\n")
       F.write("Physical Surface(\"periodicL\")={1:26};\n")
       F.write("Physical Surface(\"periodicR\")={84, 106, 128, 150, 172, 194, 216, 238, 260, 282, 304, 326, 348, 370, 392, 414, 436, 458, 480, 502, 524, 546, 568, 590, 612, 634};\n\n")
       F.write("Physical Volume(\"flowDomain\") = {1:26};\n")
       F.write("Recombine Volume \"*\";\n")
#
#
###########
#MAIN
###########
#Read in and set input parameters
inFile=sys.argv[1]
set_inputs(inFile)
#---------------------------
#   SETTINGS
#---------------------------
nAuxEll = 120   #Number of auxiliary points between the main points on the ellipses in order to construct the ellipse smooth splines
#---------------------------
#Output
geoFile=open('./ellipse.geo','w')
#
#Write in the geoFile
geo_header(geoFile)
auxEllList=geo_points(geoFile)
geo_ellipseLines(geoFile,auxEllList)
geo_otherLines(geoFile)
geo_trans(geoFile)
geo_surfaces(geoFile)
geo_setBC(geoFile)
geoFile.write('Coherence;\n')
geoFile.write('Mesh.MshFileVersion = 2.2;')   #to write .msh compatible with gmsh version2 
print('... ellipse.geo is created!')
