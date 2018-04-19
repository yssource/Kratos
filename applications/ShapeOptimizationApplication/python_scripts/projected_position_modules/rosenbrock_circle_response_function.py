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

class RosenbrockCircleResponseFunction:

	def __init__(self,model_part,response_settings):
		centerkratos = response_settings["center"].GetVector()
		self.center = [centerkratos[0],centerkratos[1]]
		self.radius = response_settings["radius"].GetDouble()
		self.value = None
		self.gradient = None

	def Initialize(self):
		pass

	def CalculateValue(self,p):
		dx = [u-v for u,v in zip(p.xRosenbrock,self.center)]
		r = norm2(dx)
		dr = r-self.radius
		self.value = -dr

	def CalculateGradient(self,p):
		dx = [u-v for u,v in zip(p.xRosenbrock,self.center)]
		r = norm2(dx)
		dr = r-self.radius
		self.gradient = [-u/norm2(dx) for u in dx]

	def GetValue(self):
		return self.value

	def GetGradient(self):
		return self.gradient