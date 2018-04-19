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

# NOW USELESS, USE CONSISTENT DAMPING INSTEAD
class FixPointResponseFunction:

	def __init__(self,model_part,response_settings):
		self.response_settings = response_settings
		self.model_part = model_part

		self.nodeId = response_settings["node_id"].GetInt()
		self.initialPosition = None
		self.value = None
		self.gradient = None


	def Initialize(self):
		pass

	def GetValue(self):
		return self.value

	def GetGradient(self):
		return self.gradient

	def CalculateValue(self,p):
		if self.initialPosition is None:
			print(p.DesignSurface.Nodes[0])
			self.initialPosition = KratosNodeToPyNode( p.DesignSurface.Nodes[self.nodeId] )

		currentPosition = KratosNodeToPyNode( p.DesignSurface.Nodes[self.nodeId] )
		self.value = norm2(minus(currentPosition,self.initialPosition))

	# the gradient is the gradient of the p-norm of the area and of the p-norm of the ratio of the length of the triangles
	def CalculateGradient(self,p):
		if self.initialPosition is None:
			self.initialPosition = KratosNodeToPyNode( p.DesignSurface.Nodes[self.nodeId] )

		currentPosition = KratosNodeToPyNode( p.DesignSurface.Nodes[self.nodeId] )
		dx = smv(minus(currentPosition,self.initialPosition))
		if norm2(dx)==0:
			gradnode = [10,10,10]
		else:
			gradnode = scalprod(1/norm2(dx),dx)

		self.gradient = b2x([ gradnode if i==self.nodeId-1 else zeros(3) for i in range(p.DesignSurface.NumberOfNodes())  ])