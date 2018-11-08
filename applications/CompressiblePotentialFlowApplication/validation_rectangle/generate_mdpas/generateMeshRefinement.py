# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.4.0 with dump python functionality
###
import os
'''
run this script with:
salome -t python generateMeshRefinement.py
'''
Number_Of_Refinements = TBD

Initial_Number_Of_Segments = TBD
Refinement_Factor = TBD

path = os.getcwd()
output_path = '/home/inigo/simulations/naca0012/07_salome/06_Rectangle/output_salome'

case = 0
NumberOfSegments_tmp = Initial_Number_Of_Segments

for i in range(Number_Of_Refinements):
    NumberOfSegments = int(NumberOfSegments_tmp)
    FarField_MeshSize = 100.0/NumberOfSegments
    print 'FarField_MeshSize = ', FarField_MeshSize
    print 'NumberOfSegments = ', NumberOfSegments

    import sys
    import salome

    salome.salome_init()
    theStudy = salome.myStudy

    import salome_notebook
    notebook = salome_notebook.NoteBook(theStudy)
    sys.path.insert( 0, r'/'+path)

    ###
    ### GEOM component
    ###

    import GEOM
    from salome.geom import geomBuilder
    import math
    import SALOMEDS


    geompy = geomBuilder.New(theStudy)

    #Create origin and axis
    O = geompy.MakeVertex(0, 0, 0)
    OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
    OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
    OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

    #Create domain
    Face_1 = geompy.MakeFaceHW(100, 100, 1)

    #Explode edges
    [Edge_1,Edge_2,Edge_3,Edge_4] = geompy.ExtractShapes(Face_1, geompy.ShapeType["EDGE"], True)

    #Generate groups
    Parts_Parts_Auto1 = geompy.CreateGroup(Face_1, geompy.ShapeType["FACE"])
    geompy.UnionIDs(Parts_Parts_Auto1, [1])
    PotentialWallCondition2D_Far_field_Auto1 = geompy.CreateGroup(Face_1, geompy.ShapeType["EDGE"])
    geompy.UnionIDs(PotentialWallCondition2D_Far_field_Auto1, [3, 6, 8, 10])

    #Add to study
    geompy.addToStudy( O, 'O' )
    geompy.addToStudy( OX, 'OX' )
    geompy.addToStudy( OY, 'OY' )
    geompy.addToStudy( OZ, 'OZ' )
    geompy.addToStudy( Face_1, 'Face_1' )
    geompy.addToStudyInFather( Face_1, Edge_1, 'Edge_1' )
    geompy.addToStudyInFather( Face_1, Edge_2, 'Edge_2' )
    geompy.addToStudyInFather( Face_1, Edge_3, 'Edge_3' )
    geompy.addToStudyInFather( Face_1, Edge_4, 'Edge_4' )
    geompy.addToStudyInFather( Face_1, Parts_Parts_Auto1, 'Parts_Parts_Auto1' )
    geompy.addToStudyInFather( Face_1, PotentialWallCondition2D_Far_field_Auto1, 'PotentialWallCondition2D_Far_field_Auto1' )

    ###
    ### SMESH component
    ###

    import  SMESH, SALOMEDS
    from salome.smesh import smeshBuilder

    smesh = smeshBuilder.New(theStudy)
    Parts_Parts_Auto1_1 = smesh.Mesh(Parts_Parts_Auto1)

    #Set mesh
    Regular_1D = Parts_Parts_Auto1_1.Segment(geom=Face_1)

    Number_of_Segments_1 = Regular_1D.NumberOfSegments(NumberOfSegments)
    Quadrangle_2D = Parts_Parts_Auto1_1.Quadrangle(algo=smeshBuilder.QUADRANGLE,geom=Face_1)

    Number_of_Segments_1.SetNumberOfSegments( 10 )

    isDone = Parts_Parts_Auto1_1.Compute()

    Regular_1D_1 = Parts_Parts_Auto1_1.Segment(geom=PotentialWallCondition2D_Far_field_Auto1)
    PotentialWallCondition2D_Far_field_Auto1_1 = Regular_1D_1.GetSubMesh()

    #Set the number of segments
    Number_of_Segments_2 = Regular_1D_1.NumberOfSegments(NumberOfSegments)

    #Mesh
    isDone = Parts_Parts_Auto1_1.Compute()

    #Split into triangles
    isDone = Parts_Parts_Auto1_1.SplitQuadObject( Parts_Parts_Auto1_1, 1 )
    Parts_Parts_Auto1_1.SetAutoColor( 1 )

    fluid_path = output_path + '/Parts_Parts_Auto1_Case_' + str(case) + '_NumberOfSegments_' + str(NumberOfSegments) + '.dat'

    far_field_path = output_path + '/PotentialWallCondition2D_Far_field_Auto1_Case_' + str(case) + '_NumberOfSegments_' + str(NumberOfSegments) + '.dat'

    #Export meshes
    try:
      Parts_Parts_Auto1_1.ExportDAT( r'/' + fluid_path)
      pass
    except:
      print 'ExportDAT() failed. Invalid file name?'
    try:
      Parts_Parts_Auto1_1.ExportDAT( r'/' + far_field_path, PotentialWallCondition2D_Far_field_Auto1_1 )
      pass
    except:
      print 'ExportPartToDAT() failed. Invalid file name?'


    ## Set names of Mesh objects
    smesh.SetName(PotentialWallCondition2D_Far_field_Auto1_1, 'PotentialWallCondition2D_Far_field_Auto1')
    smesh.SetName(Number_of_Segments_2, 'Number of Segments_2')
    smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
    smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
    smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
    smesh.SetName(Parts_Parts_Auto1_1.GetMesh(), 'Parts_Parts_Auto1')

    NumberOfSegments_tmp *= Refinement_Factor
    case +=1

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
