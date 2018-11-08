'''
This is a brief example on how mesh files from Salome in *.dat format can be read with the Salome-kratos-converter
It works with a simple 2D cantilever example, fixed on one side and a load applied on the other side
'''

# Note that this has to be on the path in order to work, or you manually specify the path
import kratos_io_utilities as kratos_utils
import global_utilities as global_utils
import os

Number_Of_Refinements = TBD

Initial_Number_Of_Segments = TBD
Refinement_Factor = TBD


output_salome_path = "/home/inigo/simulations/naca0012/07_salome/06_Rectangle/output_salome/"
output_mdpa_path = "/home/inigo/simulations/naca0012/07_salome/06_Rectangle/mdpas/"

case = 0
print('Writing mdpa...')
NumberOfSegments_tmp = Initial_Number_Of_Segments

for i in range(Number_Of_Refinements):
    NumberOfSegments = int(NumberOfSegments_tmp)
    model = kratos_utils.MainModelPart() # Main mesh object to which we will add the submeshes (Kratos Name: ModelPart)

    # Specifying the names of the submeshes (Kratos Name: SubModelPart)
    smp_dict_fluid           = {"smp_name": "Parts_Parts_Auto1"}
    smp_dict_far_field       = {"smp_name": "PotentialWallCondition2D_Far_field_Auto1"}

    file_name_fluid         = output_salome_path +  'Parts_Parts_Auto1_Case_' + str(case) + '_NumberOfSegments_' + str(NumberOfSegments) + '.dat'
    file_name_far_field     = output_salome_path +  'PotentialWallCondition2D_Far_field_Auto1_Case_' + str(case) + '_NumberOfSegments_' + str(NumberOfSegments) + '.dat'
   

    def ReadDatFile(file_name):
        valid_file, nodes, geom_entities = global_utils.ReadAndParseSalomeDatFile(os.path.join(os.getcwd(),file_name))
        if not valid_file:
            raise Exception("Invalid File!\n" + file_name)
        return nodes, geom_entities

    nodes_fluid,            geom_entities_fluid         = ReadDatFile(file_name_fluid)
    nodes_far_field,        geom_entities_far_field     = ReadDatFile(file_name_far_field)
    
    # Here we specify which Kratos-entities will be created from the general geometric entities
    mesh_dict_fluid         = {'write_smp': 1,
                           'entity_creation': {203: {'Element': {'Element2D3N': '1'}}}}
    mesh_dict_far_field     = {'write_smp': 1,
                           'entity_creation': {102: {'Condition': {'LineCondition2D2N': '0'}}}}
   
    model.AddMesh(smp_dict_fluid,           mesh_dict_fluid,            nodes_fluid,            geom_entities_fluid)
    model.AddMesh(smp_dict_far_field,       mesh_dict_far_field,        nodes_far_field,        geom_entities_far_field)
    
    mdpa_info = "mdpa for demonstration purposes"
    mdpa_file_name = output_mdpa_path + 'naca0012_Case_' + str(case) + '_NumberOfSegments_' + str(NumberOfSegments)
    
    model.WriteMesh(mdpa_file_name, mdpa_info)
    NumberOfSegments_tmp *= Refinement_Factor
    case += 1