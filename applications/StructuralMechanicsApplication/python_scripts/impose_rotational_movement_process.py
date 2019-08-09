from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeRotationalMovementProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ImposeRotationalMovementProcess(KratosMultiphysics.Process):
    """This class is used in order to impose a rigid body movement in a certain region of the problem

    This class constructs the model parts containing the constrains that enforce the rigid body movement
    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                        : "This process uses LinearMasterSlaveConstraint in order to impose an unified rtotational  movement in the given submodelpart",
            "computing_model_part_name"   : "computing_domain",
            "model_part_name"             : "please_specify_model_part_name",
            "new_model_part_name"         : "",
            "interval"                    : [0.0, 1e30],
            "transformation_settings"     : {
                "rotation_settings"         : {
                    "automatic_center"         : false,
                    "master_node_id"           : 0,
                    "center"                   : [0,0,0],
                    "axis_of_rotation"         : [0.0,0.0,0.0],
                    "constrained"              : false
                    "angle_degree"             : 0.0
                }
            }
        }
        """)

        # Detect "End" as a tag and replace it by a large number
        if settings.Has("interval"):
            if settings["interval"][1].IsString():
                if settings["interval"][1].GetString() == "End":
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        # Transforming value to string
        if settings["transformation_settings"]["rotation_settings"]["angle_degree"].IsNumber():
            angle_degree = settings["transformation_settings"]["rotation_settings"]["angle_degree"].GetDouble()
            string_angle_degree = str(angle_degree)
            #settings["transformation_settings"]["rotation_settings"].RemoveValue("angle_degree")
            #settings["transformation_settings"]["rotation_settings"].AddEmptyValue("angle_degree")
            settings["transformation_settings"]["rotation_settings"]["angle_degree"].SetString(string_angle_degree)

        # The computing model part
        computing_model_part_name = settings["computing_model_part_name"].GetString()
        self.computing_model_part = Model["Structure"].GetSubModelPart(computing_model_part_name)

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # We get the corresponding model parts
        self.model_part = Model[settings["model_part_name"].GetString()]
        new_model_part_name = settings["new_model_part_name"].GetString()
        if new_model_part_name != "":
            if self.model_part.HasSubModelPart(new_model_part_name):
                self.rotational_model_part = self.model_part.GetSubModelPart(new_model_part_name)
            else:
                self.rotational_model_part = self.model_part.CreateSubModelPart(new_model_part_name)
        else:
            settings["new_model_part_name"].SetString(settings["model_part_name"].GetString())
            self.rotational_model_part = self.model_part

        # Create the process
        rotational_parameters = KratosMultiphysics.Parameters("""{}""")
        rotational_parameters.AddValue("model_part_name", settings["model_part_name"])
        rotational_parameters.AddValue("new_model_part_name", settings["new_model_part_name"])
        rotational_parameters.AddValue("transformation_settings", settings["transformation_settings"])
        self.rotational_movement_process = StructuralMechanicsApplication.ImposeRotationalMovementProcess(self.computing_model_part, rotational_parameters)

        # Trasfering the entities
        if new_model_part_name != "":
            transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.rotational_model_part, self.model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODES)
            transfer_process.Execute()

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.rotational_movement_process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        # We activate/deactivate conditions dependeding of interval
        if self.interval.IsInInterval(current_time):
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, self.rotational_model_part.MasterSlaveConstraints)
            self.rotational_movement_process.ExecuteInitializeSolutionStep()
        else:
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, False, self.rotational_model_part.MasterSlaveConstraints)
