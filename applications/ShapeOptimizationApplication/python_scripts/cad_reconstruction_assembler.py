# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Oberbichler Thomas, https://github.com/oberbichler
#
# ==============================================================================

# Import ANurbs library
import ANurbs as an
import numpy as np

# Additional imports
import time

# ==============================================================================
class Assembler():
    # --------------------------------------------------------------------------
    def __init__(self, cad_model):
        self.cad_model = cad_model

        self.conditions = None
        self.dof_ids ={}
        self.dofs = []
        self.eqs = []
        self.lhs = np.zeros((0, 0))
        self.rhs = np.zeros((0, 0))

    # --------------------------------------------------------------------------
    def Initialize(self, conditions):
        self.conditions = conditions

        # Assign dof ids
        for conditions_face_i in self.conditions.values():
            for condition in conditions_face_i:
                for dof in condition.dof_list:
                    _ = self.__GetDofId(dof)

        num_dofs = len(self.dofs)

        # Initilize equation system
        self.lhs = np.zeros((num_dofs, num_dofs))
        self.rhs = np.zeros(num_dofs)

    # --------------------------------------------------------------------------
    def AssembleSystem(self):
        print("\n> Starting assembly....")
        start_time = time.time()

        self.lhs.fill(0)
        self.rhs.fill(0)

        for face in self.cad_model.GetByType('BrepFace'):
            face_key = face.Key()

            print("Processing face", face_key, "with", len(self.conditions[face_key]), "conditions.")

            for condition in self.conditions[face_key]:

                condition_lhs, condition_rhs, condition_dof_list = condition.CalculateLocalSystem()

                global_dof_ids = []
                for dof in condition_dof_list:
                    global_dof_ids.append(self.__GetDofId(dof))

                self.lhs[np.ix_(global_dof_ids, global_dof_ids)] += condition_lhs
                self.rhs[global_dof_ids] += condition_rhs

        print("> Finished assembly in" ,round( time.time()-start_time, 3 ), " s.")

        return self.lhs, self.rhs

    # --------------------------------------------------------------------------
    def AssembleRHS(self):
        print("\n> Starting to assemble RHS....")
        start_time = time.time()

        self.rhs.fill(0)

        for face in self.cad_model.GetByType('BrepFace'):
            face_key = face.Key()

            print("Processing face", face_key, "with", len(self.conditions[face_key]), "conditions.")

            for condition in self.conditions[face_key]:

                local_rhs, condition_dof_list = condition.CalculateRHS()

                global_dof_ids = []
                for dof in condition_dof_list:
                    global_dof_ids.append(self.__GetDofId(dof))

                self.rhs[global_dof_ids] += local_rhs

        print("> Finished assembling RHS in" ,round( time.time()-start_time, 3 ), " s.")

        return self.rhs

    # --------------------------------------------------------------------------
    def GetDofs(self):
        if self.dofs == None:
            raise Exception("> Assembler:: No dofs specified yet! First initialize assembler.")
        return self.dofs

    # --------------------------------------------------------------------------
    def GetDofIds(self):
        if self.dofs == None:
            raise Exception("> Assembler:: No dof ids specified yet! First initialize assembler.")
        return self.dof_ids

    # --------------------------------------------------------------------------
    def __GetDofId(self, dof):
        if dof not in self.dof_ids:
            self.dof_ids[dof] = len(self.dof_ids)
            self.dofs.append(dof)
        return self.dof_ids[dof]

    # --------------------------------------------------------------------------
    def __GetDof(self, index):
        return self.dofs[index]

# ==============================================================================