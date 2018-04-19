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

# linear function of the position
class TranslationResponseFunction:

	def __init__(self,model_part,response_settings):
		self.direction = response_settings["direction"].GetInt() # -1 -2 -3 1 2 3
		self.value = None
		self.gradient = None

	def Initialize(self):
		pass

	def CalculateValue(self,p):
		x = p.ReadDesignSurfaceToList()
		X = x2b(x)
		m = len(X)

		# mean of coord
		sv = 0
		for col in X:
			sv += col[abs(self.direction)-1]
		sv /= m

		if self.direction>0:
			sv = -sv

		self.value = sv

	def CalculateGradient(self,p):
		x = p.ReadDesignSurfaceToList()
		X = x2b(x)
		m = len(X)

		gradient = [ [ (i==abs(self.direction)-1)*(-1 if self.direction>0 else 1)/m for i in range(3)] for i in range(m)]
		self.gradient = b2x(gradient)

	def GetValue(self):
		return self.value

	def GetGradient(self):
		return self.gradient
