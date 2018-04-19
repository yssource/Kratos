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

class RosenbrockObjectiveResponseFunction:

	def __init__(self,model_part,response_settings):
		self.a = response_settings["factor_a"].GetDouble()
		self.b = response_settings["factor_b"].GetDouble()
		self.value = None
		self.gradient = None

	def Initialize(self):
		pass

	def CalculateValue(self,p):
		a,b = self.a, self.b
		x,y = p.xRosenbrock
		self.value = (a-x)**2+b*(y-x**2)**2

	def CalculateGradient(self,p):
		a,b = self.a, self.b
		x,y = p.xRosenbrock
		self.gradient = [2*(x-a)+b*2*(x**2-y)*2*x,  b*2*(y-x**2)]

	def GetValue(self):
		return self.value

	def GetGradient(self):
		return self.gradient
