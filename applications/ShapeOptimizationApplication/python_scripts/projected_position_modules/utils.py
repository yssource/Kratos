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

from math import sqrt,cos,sin,pi
import csv, datetime
from projected_position_modules.matrix import *
import re


def ReadParameters(obj,paramTemplates,settings):
    iAttr,iJson,iType,iOpt,iDefVal = 0,1,2,3,4
    for t in paramTemplates:
        val = __ReadValueParamDefault(t,settings)
        setattr(obj,t[iAttr],val)

def ReadAlgoParameters(p,algoParamTemplates,osa):
    iAttr,iJson,iType,iOpt,iDefVal = 0,1,2,3,4
    for t in algoParamTemplates:
        val = __ReadValueParamDefault(t,osa)
        setattr(p,t[iAttr],val)


def ReadResponseParameters(p,objectiveParamTemplates,inequalityParamTemplates,equalityParamTemplates,OptimizationSettings):
    iAttr,iJson,iType,iOpt,iDefVal = 0,1,2,3,4

    # objective
    for t in objectiveParamTemplates:
        val = __ReadValueParamDefault(t,OptimizationSettings["objectives"][0])
        setattr(p,t[iAttr],val)

    # init
    for t in inequalityParamTemplates + equalityParamTemplates: # concat
        if not hasattr(p,t[iAttr]):
            setattr(p,t[iAttr],[]) # init empty array

    # constraints
    for i in range(OptimizationSettings["constraints"].size()):
        oc = OptimizationSettings["constraints"][i]
        ctype = oc["type"].GetString()

        if ctype == "=":
            tps = equalityParamTemplates
        else:
            tps = inequalityParamTemplates

        for t in tps:
            val = __ReadValueParamDefault(t,oc)
            getattr(p,t[iAttr]).append(val)

def __ReadValueParamDefault(t,param): # t: param Templates
    iAttr,iJson,iType,iOpt,iDefVal = 0,1,2,3,4

    if param.Has(t[iJson]):
        # read value in json
        if not isinstance(t[iType],list): # wrap in list for simplicity
            t[iType] = [t[iType]]
        val = "param_value_not_read"

        for ptype in t[iType]:
            if ptype == "string" and param[t[iJson]].IsString():
                val = param[t[iJson]].GetString()
                break
            elif ptype == "double" and (param[t[iJson]].IsDouble() or param[t[iJson]].IsInt() ):
                val = param[t[iJson]].GetDouble()
                break
            elif ptype == "int" and param[t[iJson]].IsInt():
                val = param[t[iJson]].GetInt()
                break
            elif ptype == "bool" and param[t[iJson]].IsBool():
                val = param[t[iJson]].GetBool()
                break
            elif ptype == "vector" and param[t[iJson]].IsVector():
                vect = param[t[iJson]].GetVector()
                val = [vect[i] for i in range(vect.Size())]
                break

        if val=="param_value_not_read": # value has not been read
            raise ValueError("value not read, probably wrong type: "+t[iJson]+": "+str(t[iType]))

    else:
        # check default value
        if not t[iOpt]:
            raise ValueError(t[iJson]+" not specified")
        else:
            val = t[iDefVal]
    return val

def BuildProjectNormalsIds(p):
    now = datetime.datetime.now()
    p.timestamp = str(now.year)+"{:02}".format(now.month)+"{:02}".format(now.day)+"_"+"{:02}".format(now.hour)+"{:02}".format(now.minute)+"{:02}".format(now.second)

    projIds = []
    if p.objectiveProjectNormals:
        projIds.append(p.objectiveId)
    for i,b in enumerate(p.inequalityProjectNormals):
        if b:
            projIds.append(p.inequalityIds[i])
    for i,b in enumerate(p.equalityProjectNormals):
        if b:
            projIds.append(p.equalityIds[i])
    return projIds

def BuildRunId(p):
    geom = p.OptimizationSettings["design_variables"]["optimization_model_part_name"].GetString()
    obj = "o-"+p.objectiveId
    ineq = "-".join(["i",*p.inequalityIds]) if len(p.inequalityIds)>0 else ""
    eq = "-".join(["e",*p.equalityIds]) if len(p.equalityIds)>0 else ""
    runid = "_".join([p.timestamp , geom, obj , ineq ,  eq ] )

    strToReplace = [
        ["strain_energy","strain"],
        ["mesh_control","mshctrl"],
        ["translation","transl"],
        ["packaging",""],
        ["inner",""],
        ["outer",""],
        ["__","_"]
    ]
    for s,r in strToReplace:
        runid = runid.replace(s,r)

    if "packaging_inner" in p.inequalityIds:
        for i in range(len(p.inequalityIds)):
            if "packaging_inner" == p.OptimizationSettings["constraints"][i]["identifier"].GetString():
                packagingName = p.OptimizationSettings["constraints"][i]["kratos_response_settings"]["file_name"].GetString()
                packagingName = packagingName.replace(".","")
                runid = runid.replace("packaging_inner-packaging_outer",packagingName)
                break

    if runid[-1]=="_":
        runid = runid[:-1]
    return runid

def ReadDesignSurfaceToList(DesignSurface):
    list_of_values = []
    for node in DesignSurface.Nodes:
        list_of_values.append(node.X)
        list_of_values.append(node.Y)
        list_of_values.append(node.Z)
    return list_of_values

def simpleDistanceToFeasibleDomain(lInequality0, lEquality0):
    return max( [max(l,0) for l in lInequality0] + [abs(l) for l in lEquality0] + [0] ) # add [0] for case with no constraint

def BuildNodeIdIndex(DesignSurface):
    nodeIds = []
    for node in DesignSurface.Nodes:
        nodeIds.append(node.Id)

    idIndex = zeros(max(nodeIds)+1)
    i = 0
    for node in DesignSurface.Nodes:
        idIndex[node.Id] = i
        i += 1

    return nodeIds,idIndex

def KratosNodeToPyNode(krnodes):
    if not isinstance(krnodes,list):
        krnodes = [krnodes]
    pynodes = []
    for krn in krnodes:
        pynodes.append([krn.X,krn.Y,krn.Z])
    if len(pynodes)==1:
        return pynodes[0]
    return pynodes

def deltaXBoolToValueGradient(dx,isFeasible):
    if norm2(dx)==0:
        value = -10
        gradient = [0.0001 for _ in dx] #avoid division by 0
        return value,gradient

    if isFeasible:
        value = min(-normmin3d(dx),-0.01)
    else:
        value = norminf3d(dx)

    gradient = scalprod( -abs(value)/norm2(dx)**2 ,dx)
    return value,gradient


# def valueGradientToDeltaXBool(value,gradient):
#     isFeasible = value<=0
#     dx = scalprod(-abs(value)/norm2(gradient)**2,gradient)
#     return dx,isFeasible

# return min 2-norm of the node that is not 0 (=distance to the infeasible domain)

def normmin3d(dx):
    if norm2(dx)==0:
        return 0
    nodes = x2b(dx)
    return min(norm2(node) for node in nodes if norm2(node)>0)


def rigidBodyTransform(nodes,center,diameter,rotation):

    maxNode = [ max(n[i] for n in nodes) for i in range(3)  ]
    minNode = [ min(n[i] for n in nodes) for i in range(3)  ]
    center0 = scalprod(1/2, plus(maxNode,minNode) )
    diameter0 = max(smv(minus(maxNode,minNode)))

    # nodes centered on 0
    nodes = minus(nodes, center0 )

    #scale
    if diameter is not None:
        nodes = [ scalprod( diameter/diameter0  ,n) for n in nodes]

    #rotate
    if rotation is not None:
        rx,ry,rz = scalprod(pi/180, rotation)
        Rx = [ # vector of COLUMNS!
            [1,0,0],
            [0,cos(rx),sin(rx)],
            [0,-sin(rx),cos(rx)]
        ]
        Ry = [
            [cos(ry),0,-sin(ry)],
            [0,1,0],
            [sin(ry),0,cos(ry)]
        ]
        Rz = [
            [cos(rz),sin(rz),0],
            [-sin(rz),cos(rz),0],
            [0,0,1]
        ]
        R = prod(Rx,prod(Ry,Rz))
        nodes = prod(R,nodes)

    # center
    if center is None:
        nodes = plus(nodes,center0)
    else:
        nodes = plus(nodes,center)

    return nodes

def readObjFile(fileName):
    groupToNode = lambda m: [float(m.group(i)) for i in range(1,4)]
    groupToTriangle = lambda m: [int(m.group(i)) for i in range(1,4)]

    nodes,triangles = [],[]
    patternNode = re.compile("(\S+) (\S+) (\S+)$")
    patternTriangle = re.compile("(\S+)//\S+ (\S+)//\S+ (\S+)//\S+")

    with open(fileName,"r") as file:
        for line in file:
            if line[0:2]=="v ": #vertexscaleCenterDiameter(self.nodes,self.center,self.diameter)
                m = patternNode.search(line)
                nodes.append( groupToNode(m) )
            elif line[0:2]=="f ":
                m = patternTriangle.search(line)
                tr = groupToTriangle(m)
                tr = [u-1 for u in tr] # begin at 0
                triangles.append( tr )
    return nodes,triangles

def computeNormals(nodes,triangles):
    normals = []
    for tr in triangles:
        n1,n2,n3 = nodes[tr[0]], nodes[tr[1]], nodes[tr[2]]
        v1 = minus(n2,n1)
        v2 = minus(n3,n2)
        cr = cross(v1,v2)
        normals.append(scalprod(1/norm2(cr),cr))
    return normals

def getTransformedObjFile(fileName,center,diameter,rotation):
    nodes,triangles = readObjFile(fileName)
    nodes = rigidBodyTransform(nodes,center,diameter,rotation)
    normals = computeNormals(nodes,triangles)
    return nodes,triangles,normals

def writeObjFile(nodes,triangles,normals,filePath,addNormalsToNode=False):
    if addNormalsToNode: # check if normals in right direction
        dn = zermat(nodes)
        for t,n in zip(triangles,normals):
            for idn in t:
                dn[idn] = smv(plus(dn[idn],n))
                dn[idn] = scalprod(1/norm2(dn[idn])*5,dn[idn])
        for i in range(len(nodes)):
            nodes[i] = smv(plus(nodes[i],dn[i]))

    with open(filePath, 'w') as file:
        file.write("o scaled_obj")
        for n in nodes:
            file.write("\nv {:.5f} {:.5f} {:.5f}".format(*n))
        for nml in normals:
            file.write("\nvn {:.5f} {:.5f} {:.5f}".format(*nml))
        for i,tr in enumerate(triangles):
            tr = [u+1 for u in tr]
            s = "\nf {}//i {}//i {}//i".format(*tr)
            s = s.replace("i",str(i+1))
            file.write(s)

def testObjFile(fileName):
    center = [0,0,-10]
    diameter = None
    rotation = [0,0,90]
    outputFileName = fileName.replace(".obj","_out.obj")
    normalsFileName = fileName.replace(".obj","_norm.obj")

    nodes,triangles,normals = getTransformedObjFile(fileName,center,diameter,rotation)
    writeObjFile(nodes,triangles,normals,outputFileName)
    writeObjFile(nodes,triangles,normals,normalsFileName,addNormalsToNode=True)