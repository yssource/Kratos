sed 's|'"Number_Of_Refinements = $Number_Of_Refinements"'|'"Number_Of_Refinements = TBD"'|g' -i /$Generate_Mesh_File_Path \
                                                                                                /$Salome_Converter_File_Path /$Mesh_Refinement_File_Path


sed 's|'"Initial_Number_Of_Segments = $Initial_Number_Of_Segments"'|'"Initial_Number_Of_Segments = TBD"'|g' -i /$Generate_Mesh_File_Path \
                                                                        /$Salome_Converter_File_Path /$Mesh_Refinement_File_Path

sed 's|'"Refinement_Factor = $Refinement_Factor"'|'"Refinement_Factor = TBD"'|g' -i /$Generate_Mesh_File_Path \
                                                                        /$Salome_Converter_File_Path /$Mesh_Refinement_File_Path