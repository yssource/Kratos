import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication
from numpy import *
import itertools
import matplotlib.pyplot as plt

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeLiftProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ComputeLiftProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        default_parameters = KratosMultiphysics.Parameters("""
            {
                
                "model_part_name" : "please specify the main model part",
                "mesh_id": 0,
                "velocity_infinity": [1.0,0.0,0],
                "reference_area": 1,
                "problem_name": "embedded_potential_flow"
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)
        self.problem_name=settings["problem_name"].GetString()
        self.fluid_model_part=Model.GetModelPart(settings["model_part_name"].GetString()).GetRootModelPart() 
        self.result_force=KratosMultiphysics.Vector(3)  
        self.process=KratosMultiphysics.CompressiblePotentialFlowApplication.ComputeLiftLevelSetProcess(self.fluid_model_part,self.result_force)


    def ExecuteFinalizeSolutionStep(self):
        print("wip_compute_lift_level_set_process")
        self.process.Execute()
        x_upper=[]
        cp_upper=[]
        x_lower=[]
        cp_lower=[]
        for element in self.fluid_model_part.Elements:
            if element.Is(KratosMultiphysics.BOUNDARY) and element.IsNot(KratosMultiphysics.STRUCTURE):
                gp=element.GetValue(KratosMultiphysics.BODY_FORCE)#provisional name
                pressure=element.GetValue(KratosMultiphysics.PRESSURE)
                normal=element.GetValue(KratosMultiphysics.NORMAL)
                if pressure==0.0:
                    print(element.Id)
                if normal[1]<=0:
                    x_upper.append(gp[0])
                    cp_upper.append(pressure)
                else:
                    x_lower.append(gp[0])
                    cp_lower.append(pressure)
        max_x=max(max(x_upper),max(x_lower))
        min_x=min(min(x_upper),min(x_lower))
        for i in range(0,len(x_upper)):
            x_upper[i]=(x_upper[i]-min_x)/abs(max_x-min_x)
        for i in range(0,len(x_lower)):
            x_lower[i]=(x_lower[i]-min_x)/abs(max_x-min_x)
        # for i in range(len(x_upper)):
            # if x_upper[i]>0.4 and x_upper[i]<0.6:
                # print(x_upper[i],cp_upper[i])
        print("Cl:",self.result_force[1]) 
        print("Cd:",self.result_force[0])   
        # self.fluid_model_part.SetValue(KratosMultiphysics.FRICTION_COEFFICIENT,self.result_force[1])
        # for node in self.fluid_model_part.GetSubModelPart("KuttaLS").Nodes:
        #     Cl_jump=node.GetValue(KratosMultiphysics.TEMPERATURE)
        # print("Cl_jump:",Cl_jump)
        with open('cp_distribution_naca0012.dat') as cp_ref:
            lines=cp_ref.readlines()
            x_ref=[]
            cp_ref=[]
            for line in lines:
                x_ref.append(float(line.split(' ')[0]))
                cp_ref.append(float(line.split(' ')[1]))
        max_x=max(x_ref)
        min_x=min(x_ref)
        for i in range(0,len(x_ref)):
            x_ref[i]=(x_ref[i]-min_x)/abs(max_x-min_x)
        plt.plot(x_upper,cp_upper,'.',label='Upper surface')
        plt.plot(x_lower,cp_lower,'r.',label='Lower surface')
        plt.plot(x_ref, cp_ref ,'k.')
        title="Cl: %.5f, Cd: %.5f" % (self.result_force[1],self.result_force[0])
        plt.title(title)
        plt.legend()
        # plt.ylim(-1.8, 1.2)
        # plt.xlim(-0.2, 1.2)
        plt.gca().invert_yaxis()
        # plt.show()
        plt.savefig(self.problem_name+'.png', bbox_inches='tight')
        plt.close('all')
    
