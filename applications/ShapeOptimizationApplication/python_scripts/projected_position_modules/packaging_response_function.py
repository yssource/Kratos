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
from projected_position_modules.projection import projHs,projHp

# for one packaging domain, there are 2 constraints functions, the inner and the outer one, both implemented here
class PackagingResponseFunction:

	def __init__(self,model_part,response_settings,strInnerOrOuter,identifier):
		self.model_part = model_part
		self.response_settings = response_settings
		self.strInnerOrOuter = strInnerOrOuter
		self.identifier = identifier # the identifier is of the form "package_inner_geomname"
		self.identifierRaw = self.identifier.replace("inner","").replace("outer","").replace("packaging","").replace("_","")

		self.value = None
		self.gradient = None
		self.lateInitDone = False
		self.lastIter = -1
		# see also LateInit

	# late initialization ot have access to p, the instance of the algorithm class
	def __LateInit(self,p):
		# the response parameters need to be defined only for the inner or the outer response functions, not both
		self.response_settings = self.__ShareResponseSettingsBetweenInstance(p)
		paramResponse = [
			["fileName","file_name","string",False],
			["buffer","buffer","double",True,0]
		]
		ReadParameters(self,paramResponse,self.response_settings)
		# set geometry name
		if self.fileName.endswith(".obj"):
			self.geometryName = "mesh"
			self.objFileName = self.fileName
		else:
			self.geometryName = self.fileName
		self.center,self.diameter,self.rotation = None,None,None

		# read geometry specific parameters
		if self.geometryName == "mesh":
			paramGeom = [
				["center","center","vector",True,None],
				["diameter","diameter","double",True,None],
				["rotation","rotation","vector",True,None]
			]
			ReadParameters(self,paramGeom,self.response_settings)
		elif self.geometryName == "sphere_analytic":
			paramGeom = [
				["center","center","vector",False],
				["diameter","diameter","double",False],
			]
			ReadParameters(self,paramGeom,self.response_settings)
			self.objFileName = "sphere.obj"
			raise ValueError("SHOULD NOT BE USED ANYMORE, use .obj file instead")
		elif self.geometryName == "halfspace":
			paramGeom = [
				["position","position","vector",False],
				["direction","direction","vector",False],
			]
			ReadParameters(self,paramGeom,self.response_settings)
			self.direction = scalprod(1/norm2(self.direction),self.direction)
			self.center = self.position # for scaled obj
			self.objFileName = "square.obj"
		else:
			raise NameError("geometryName not recognized: for"+self.fileName+", or give proper mesh file ")

		self.nodes,self.triangles,self.normals = getTransformedObjFile(self.objFileName,self.center,self.diameter,self.rotation)
		scaledFilePath = p.resultsFolder+"/scaled_"+self.objFileName
		writeObjFile(self.nodes,self.triangles,self.normals,scaledFilePath)

	def __ShareResponseSettingsBetweenInstance(self,p):
		if not hasattr(p,"packagingData"):
			p.packagingData = {}
		if not self.identifierRaw in p.packagingData:
			p.packagingData[self.identifierRaw] = {"last_iter":-1, "dxs": None}

		if "response_settings" in p.packagingData[self.identifierRaw]:
			if self.response_settings.Has("file_name") and p.packagingData[self.identifierRaw]["response_settings"].Has("file_name"):
				raise ValueError("The parameters of packaging should be given only once either in packaging_inner or in packaging_outer. Here param defined twice")
			return p.packagingData[self.identifierRaw]["response_settings"]
		else:
			if not self.response_settings.Has("file_name"):
				raise ValueError("file_name missing")
			p.packagingData[self.identifierRaw]["response_settings"] = self.response_settings
			return self.response_settings

	def Initialize(self):
		pass

	def CalculateValue(self,p):
		self.value = self.__InequalityPackaging(p,"value")

	def CalculateGradient(self,p):
		self.gradient = self.__InequalityPackaging(p,"gradient")

	def GetValue(self):
		return self.value

	def GetGradient(self):
		return self.gradient

	def __InequalityPackaging(self,p,strValueOrGradient):
		# late initialization
		if not self.lateInitDone:
			self.__LateInit(p)
			self.lateInitDone = True

		if self.geometryName=="sphere_analytic":
			deltaXBoolPerNode = lambda node: self.__DeltaXBoolPerNodeSphere(node,self.center,self.radius,self.buffer,p.stepLength)
		elif self.geometryName == "halfspace":
			deltaXBoolPerNode = lambda node: self.__DeltaXBoolPerNodeHalfspace(node,self.position,self.direction,self.buffer,p.stepLength)
		elif self.geometryName == "mesh":
			deltaXBoolPerNode = lambda node: self.__DeltaXBoolPerNodeMesh(node,self.nodes,self.triangles,self.normals,self.buffer,p.stepLength)
		else:
			raise NameError("wrong name")

		# compute only once per iteration, and store result
		if p.packagingData[self.identifierRaw]["last_iter"]!=p.iterCounter:
			p.packagingData[self.identifierRaw]["last_iter"] = p.iterCounter
			x = p.ReadDesignSurfaceToList()
			dxInnerToBorder,dxOuterToBorder = self.__DeltaXInnerOuter(x,deltaXBoolPerNode,self.buffer,p.stepLength)
			p.packagingData[self.identifierRaw]["dxs"] = [dxInnerToBorder,dxOuterToBorder]
		else:
			dxInnerToBorder,dxOuterToBorder = p.packagingData[self.identifierRaw]["dxs"]

		if self.strInnerOrOuter=="inner":
			dxBorder = dxInnerToBorder
			isFeasible = False
		elif self.strInnerOrOuter=="outer" :
			dxBorder = dxOuterToBorder
			isFeasible = True
		else:
			raise NameError("wrong name")

		value,gradient = deltaXBoolToValueGradient(dxBorder,isFeasible)

		if strValueOrGradient=="gradient":
			return gradient
		elif strValueOrGradient=="value":
			return value
		else:
			raise NameError("wrong name")

	def __DeltaXInnerOuter(self,x,deltaXBoolPerNode,buffer,stepLength):
		nodes = x2b(x)
		distOuter = 1.2*stepLength
		dxInnerToBorder = []
		dxOuterToBorder = []

		for node in nodes:

			deltaNodeToBorder,isInside = deltaXBoolPerNode(node) # <0 if outside, >0 if inside
			distNodeToBorder = norm2(deltaNodeToBorder)

			if isInside:
				dxInnerToBorder.append(deltaNodeToBorder)
				dxOuterToBorder.append(zeros(3))
			elif distNodeToBorder<distOuter:
				dxInnerToBorder.append(zeros(3))
				dxOuterToBorder.append(deltaNodeToBorder)
			else:
				dxInnerToBorder.append(zeros(3))
				dxOuterToBorder.append(zeros(3))
		return b2x(dxInnerToBorder), b2x(dxOuterToBorder)

#########################################

	def __DeltaXBoolPerNodeSphere(self,node,center,radius,buffer,stepLength):
		radius = radius+buffer

		nodeRel = minus(node,center)
		rNode = norm2(nodeRel)
		if rNode<radius:
			isInside = True
		else:
			isInside = False
		eRadial = scalprod(1/rNode,nodeRel)
		dxToBorder = scalprod(radius-rNode,eRadial)
		return dxToBorder,isInside

	def __DeltaXBoolPerNodeHalfspace(self,node,position,direction,buffer,stepLength):
		position = smv(plus(position,scalprod(-buffer,direction)))

		xProj = projHp(node,[position],[direction]) # careful, projection onto hyperplane, not halfspace
		dxToBorder = smv(minus(xProj,node))
		if dot(direction,minus(node,position))>0 :
			isInside = True
		else:
			isInside = False
		return dxToBorder,isInside

	def __DeltaXBoolPerNodeMesh(self,node,meshNodes,triangles,normals,buffer,stepLength):
		#nearest neighbor
		nbrNode = meshNodes[0]
		nbrNodeId = 0
		nbrDist = norm2(minus(node,nbrNode))
		for i,mn in enumerate(meshNodes):
			if norm2(minus(node,mn))<nbrDist:
				nbrNode = mn
				nbrNodeId = i
				nbrDist = norm2(minus(node,mn))

		# get triangles that include the nearest neighbor. together they form a corner
		normalsCorner = []
		for i,tr in enumerate(triangles):
			if nbrNodeId in tr:
				normalsCorner.append(normals[i])

		# check if inside or outside corner
		isInside = True
		minDistHs = 1e9
		minDistNormal = normalsCorner[0]
		for i,nc in enumerate(normalsCorner):
			distHs = dot(minus(node,nbrNode),nc)
			if abs(distHs)<minDistHs:
				minDistHs = abs(distHs)
				minDistNormal = nc
			if distHs>buffer:
				isInside = False

		# project onto corner
		if isInside:
			nodeProj = projHs(node,[nbrNode],[minDistNormal])
		else:
			if nbrDist>buffer+1.2*stepLength:
				nodeProj = nbrNode # simpler formula
			else:
				pos = [nbrNode for _ in normalsCorner]
				dir = [ scalprod(-1,nc) for nc in normalsCorner]
				nodeProj = projHs(node,pos,dir)

		dxToBorder = smv(minus(nodeProj,node))

		if buffer!=0:
			if isInside:
				dxToBorder = smv(plus(dxToBorder, scalprod( buffer/norm2(dxToBorder) ,dxToBorder) ) )
			else:
				dxToBorder = smv(plus(dxToBorder, scalprod( -buffer/norm2(dxToBorder) ,dxToBorder) ) )

		return dxToBorder, isInside


