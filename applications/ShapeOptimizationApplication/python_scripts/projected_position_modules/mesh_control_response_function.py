# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division


from projected_position_modules.matrix import *
from projected_position_modules.utils import *

class MeshControlResponseFunction:

	def __init__(self,model_part,response_settings):
		self.power = response_settings["power"].GetDouble()
		self.value = None
		self.gradient = None

	def Initialize(self):
		pass

	def GetValue(self):
		return self.value

	def GetGradient(self):
		return self.gradient

	def CalculateValue(self,p):
		self.value = 1
		if p.isEnforcingFinalFeasibility:
			self.value = -1

	# the gradient is the gradient of the p-norm of the area and of the p-norm of the ratio of the length of the triangles
	def CalculateGradient(self,p):
		gradPowArea = self.__DiscreteGradient(self.__Pow(self.__AreaTriangle),p)
		gradPowRatio = self.__DiscreteGradient(self.__Pow(self.__RatioTriangle),p)
		gradArea = self.__DiscreteGradient(self.__AreaTriangle,p)

		gradient = plus( scalprod(1/norm2(gradPowArea),gradPowArea) , scalprod(1/norm2(gradPowRatio),gradPowRatio) )

		self.gradient = scalprod(1/1000,gradient) # low norm to be always far from the feasible domain


	def __AreaTriangle(self,nodes):
		if not isinstance(nodes,list): # kratos_nodes
			nodes = KratosNodeToPyNode(nodes)
		x1,x2,x3 = nodes
		v1 = minus(x2,x1)
		v2 = minus(x3,x1)
		v3 = minus(x3,x2)
		return 1/2*norm2(cross(v1,v2))

	def __RatioTriangle(self,nodes):
		if not isinstance(nodes,list): # kratos_nodes
			nodes = KratosNodeToPyNode(nodes)
		x1,x2,x3 = nodes
		v1 = minus(x2,x1)
		v2 = minus(x3,x1)
		v3 = minus(x3,x2)
		norms = [norm2(v1),norm2(v2),norm2(v3)]
		return max(norms)-min(norms)

	# return handle
	def __Pow(self,funhandle):
		return lambda x: funhandle(x)**self.power

	# compute the discrete gradient of the function handle
	def __DiscreteGradient(self,valuefun,p):
		pertubation = 0.0001
		gradient = [zeros(3) for _ in range(len(p.nodeIds))]
		for cdn in p.DesignSurface.Conditions: # conditions are triangles
			krnodes = cdn.GetNodes()
			nodes = KratosNodeToPyNode(krnodes)
			initialValue = valuefun(nodes)
			for num in range(3):
				nId = krnodes[num].Id
				nIndex = p.idIndex[nId]
				for dim in range(3):
					nodes[num][dim] += pertubation
					modifiedValue = valuefun(nodes)
					gradient[nIndex][dim] += (modifiedValue - initialValue)/pertubation # cumulative
					nodes[num][dim] -= pertubation
		return b2x(gradient)
