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

class RosenbrockLineResponseFunction:

	def __init__(self,model_part,response_settings):
		pos1kratos = response_settings["position_1"].GetVector()
		pos2kratos = response_settings["position_2"].GetVector()
		self.pos1 = [pos1kratos[0],pos1kratos[1]]
		self.pos2 = [pos2kratos[0],pos2kratos[1]]
		delta = [u-v for u,v in zip(self.pos2,self.pos1)]
		dir = [-delta[1],delta[0]]
		self.dir = [u/norm2(dir) for u in dir]

		self.value = None
		self.gradient = None

	def Initialize(self):
		pass

	def CalculateValue(self,p):
		self.value = -dot(minus(p.xRosenbrock,self.pos1),self.dir)

	def CalculateGradient(self,p):
		self.gradient = [-u for u in self.dir]

	def GetValue(self):
		return self.value

	def GetGradient(self):
		return self.gradient
