Work_Dir=/home/inigo/simulations/naca0012/07_salome/06_Rectangle

Generate_Mesh_File_Path=$PWD/generate_mdpas/generateMeshRefinement.py
Salome_Converter_File_Path=$PWD/generate_mdpas/use_converter.py
Mesh_Refinement_File_Path=$PWD/MeshRefinement.py

sed 's|'"Number_Of_Refinements = TBD"'|'"Number_Of_Refinements = $Number_Of_Refinements"'|g' -i /$Generate_Mesh_File_Path \
                                                            /$Salome_Converter_File_Path /$Mesh_Refinement_File_Path


sed 's|'"Initial_Number_Of_Segments = TBD"'|'"Initial_Number_Of_Segments = $Initial_Number_Of_Segments"'|g' -i /$Generate_Mesh_File_Path \
                                                            /$Salome_Converter_File_Path /$Mesh_Refinement_File_Path

sed 's|'"Refinement_Factor = TBD"'|'"Refinement_Factor = $Refinement_Factor"'|g' -i /$Generate_Mesh_File_Path \
                                                            /$Salome_Converter_File_Path /$Mesh_Refinement_File_Path