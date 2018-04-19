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

from math import sqrt,inf
from projected_position_modules.matrix import *
import csv
from copy import deepcopy

# called in the algorithm, p is the instance of the algorithm class
def plusSphere(eObjective,lInequality,eInequality,lEquality,eEquality,p):

    target = 1
    tolerance = 0.03
    nInterval = 5
    # check length init
    lObjective = -10 # too low value
    threshold = 10 # too large value
    sInit = normCoord(*filterShareThreshold(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,threshold,p),p)
    printListValues(["plusSphere> ",["li0:[{}]","{:.2f} ",lInequality],["le0:[{}]","{:.2f} ",lEquality],["                (sInit:{})","{:.2f}",sInit]])
    if sInit>1: # reduce with threshold
        funSearchStepLength = lambda threshold: normCoord(*filterShareThreshold(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,threshold,p),p)
        lMin = 0
        lMax = 1.2
        threshold,sDx = searchDichotomy(funSearchStepLength,target,tolerance,lMin,lMax,nInterval,p)
    else: # fill with objective
        if not p.isEnforcingFinalFeasibility:
            funSearchStepLength = lambda lObj: normCoord(*filterShareThreshold(lObj,eObjective,lInequality,eInequality,lEquality,eEquality,threshold,p),p)
            lMin = p.objectiveMinShare
            lMax = 1.3
            lObjective,sDx = searchDichotomy(funSearchStepLength,target,tolerance,lMin,lMax,nInterval,p)
        else:
            sDx = sInit # take step as it is

    # actual values:
    lObjective,eObjective,lInequality,eInequality,lEquality,eEquality = filterShareThreshold(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,threshold,p)

    if p.isEnforcingFinalFeasibility:
        print("plusSphere> enforce final feasibility iter")
        for i in range(len(lEquality)):
            lEquality[i] *= p.reduceOscillationFactor # prevent oscillation
        for i in range(len(lInequality)):
            lInequality[i] += 0.1 # speed up towards the feasible domain

    printListValues(["plusSphere> ",[" li:[{}] ","{:.2f} ",lInequality], [" le:[{}] ","{:.2f} ",lEquality], ["lo:{} ","{:.2f}",lObjective], ["       (sDx:{})","{:.2f}",sDx]])

    x = sumCoord(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,p)
    return x,lObjective,lInequality,lEquality,sInit,sDx

# FOR OLD ALGORITHM
def plusSphereOld(eObjective,lInequality,eInequality,lEquality,eEquality,p):
    # read function parameters
    if p.sphereSum == "projection":
        sumCoord = projCoord
    elif p.sphereSum == "addition":
        sumCoord = addCoord
    else:
        raise NameError("wrong name")
    if p.sphereNorm == "geom_norminf3d":
        normCoord = geomNormInf3d
    elif p.sphereNorm == "coord_norm1":
        normCoord = coordNorm1
    else:
        raise NameError("wrong name")

    lObjective = -10 # too low value
    # check length init
    sInit = normCoord(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,sumCoord)
    printListValues(["plusSphereOld> ",["li0:[{}]","{:.2f} ",lInequality],["le0:[{}]","{:.2f} ",lEquality],["                (sInit:{})","{:.2f}",sInit]])
    if sInit>1: # if the projection of 0 onto the hs and hp is outside of the ball, no solution exist. return
        xInit = sumCoord(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality)
        lObjective = dot(xInit,scalprod(1/norm2(eObjective)**2,eObjective))
        print("projection.projBallHpHs> sInit={:.2f}".format(sInit))
        return (xInit,lObjective,sInit,sInit,lInequality)

    target = 1
    tolerance = 0.03
    nInterval = 5
    lMin = -0.2
    lMax = 1.2
    funSearchStepLength = lambda lObj: normCoord(lObj,eObjective,lInequality,eInequality,lEquality,eEquality,sumCoord)
    lObjective,sDx = searchDichotomy(funSearchStepLength,target,tolerance,lMin,lMax,nInterval,p)

    printListValues(["plusSphere> ",[" li:[{}] ","{:.2f} ",lInequality], [" le:[{}] ","{:.2f} ",lEquality], ["lo:{} ","{:.2f}",lObjective], ["       (sDx:{})","{:.2f}",sDx]])

    x = sumCoord(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality)

    lXConstraint = zeros(len(lInequality))
    for i in range(len(lInequality)):
        lXConstraint[i] = max(lInequality[i],dot(x,scalprod(1/norm2safe(eInequality[i])**2,eInequality[i])))

    return x,lObjective,sInit,sDx,lXConstraint


def normCoord(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,p):
    if p.sphereNorm == "geom_norminf3d":
        return geomNormInf3d(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,p)
    elif p.sphereNorm == "coord_norm1":
        return coordNorm1(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,p)
    else:
        raise NameError("wrong name")

#NORMCOORD
def geomNormInf3d(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,p):
    x = sumCoord(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,p)
    return norminf3d(x)

#NORMCOORD
def coordNorm1(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,p):
    raise NameError("COORDNORM1 SHOULD USUALLY NOT BE USED, JUST FOR DEMONSTRATION PURPOSE")
    nm = lObjective if lObjective>0 else 0
    for i in range(len(lEquality)):
        nm += abs(lEquality[i])
    for i in range(len(lInequality)):
        if lInequality[i]>0:
            nm += lInequality[i]
    return nm

def sumCoord(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,p):
    if p.sphereSum == "projection":
        return projCoord(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality)
    elif p.sphereSum == "addition":
        return addCoord(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality)
    else:
        raise NameError("wrong name")

def addCoord(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality):
    print("ADDITION SHOULD NOT BE USED ANYMORE, USE PROJECTION INSTEAD")
    n = len(eObjective)
    x = zeros(n)
    if lObjective>0:
        x = scalprod(lObjective,eObjective)
    for i in range(len(lInequality)):
        if lInequality[i]>0:
            x = plus(x,scalprod(lInequality[i],eInequality[i]))
    for i in range(len(lEquality)):
        x = plus(x,scalprod(lEquality[i],eEquality[i]))
    return smv(x)

def projCoord(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality):
    try:
        lObjective,eObjective,lInequality,eInequality,lEquality,eEquality = deepcopy((lObjective,eObjective,lInequality,eInequality,lEquality,eEquality))
        hpPos,hpDir,hsPos,hsDir = coordToHpHs(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality)
        n = len(eObjective)
        return projHpHs(zeros(n),hpPos,hpDir,hsPos,hsDir)
    except:
        print("projHpHs fail, fallback on addition")
        return addCoord(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality)

# apply min max share and threshold
def filterShareThreshold(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality,threshold,p):
    #apply min/maxShare
    lo = max(lObjective,p.objectiveMinShare)
    li = [min(l,maxshare) for l,maxshare in zip(lInequality,p.inequalityMaxShare)]
    le = [inInterval(l,maxshare) for l,maxshare in zip(lEquality,p.equalityMaxShare)]

    #apply threshold
    lo = min(lo,threshold)
    li = [min(l,threshold) for l in li]
    le = [inInterval(l,max(threshold,0)) for l in le] # if threshold negative, lEquality = 0

    if p.isEnforcingFinalFeasibility: # remove objective contribution
        lo = -10 # too low value

    return lo,eObjective,li,eInequality,le,eEquality

# convert from length-direction format to position-direction format
def coordToHpHs(lObjective,eObjective,lInequality,eInequality,lEquality,eEquality):
    hsPos = []; hsDir = []; hpPos = []; hpDir = []
    for i in range(len(lInequality)):
        if norm2(eInequality[i])>0 and lInequality[i]>-3: # do not add halfspaces far away (l<-3)
            hsPos.append(scalprod(lInequality[i],eInequality[i]))
            hsDir.append(scalprod(1/norm2(eInequality[i]),eInequality[i]))
    for i in range(len(lEquality)):
        if norm2(eEquality[i])>0:
            hpPos.append(scalprod(lEquality[i],eEquality[i]))
            hpDir.append(scalprod(1/norm2(eEquality[i]),eEquality[i]))
    hsPos.append(scalprod(lObjective,eObjective))
    hsDir.append(scalprod(1/norm2(eObjective),eObjective))
    return hpPos,hpDir,hsPos,hsDir

#####################
# PROJECTION

# project orthogonally x on the halfspaces and hyperplanes
def projHpHs(xOriginal0,hpPos0,hpDir0,hsPos0,hsDir0): # hss (point go through, direction) and hps (point go through,direction)
    hpPos0, hpDir0 = filterPosDir(hpPos0,hpDir0)
    hsPos0, hsDir0 = filterPosDir(hsPos0,hsDir0)

    nHp = len(hpDir0)
    nHs = len(hsDir0)
    if nHp==0 and nHs==0:
        return xOriginal0

    # convert to reduced coordinates
    basis = basisGramSchmidt(horzcat(hpDir0,hsDir0))
    hpPos = xb(minus(hpPos0,xOriginal0),basis)
    hpDir = xb(hpDir0,basis)
    hsPos = xb(minus(hsPos0,xOriginal0),basis)
    hsDir = xb(hsDir0,basis)

    # project on the HP
    xHp = projHp(zeros(size2(basis)),hpPos,hpDir)

    # conversion of HS to the reduced coordinates:  project HSposition on the other HPs, along its own HP
    for i in range(nHs):
        hsPos[i] = projHp(hsPos[i],horzcat(hpPos,hsPos[i]),horzcat(hpDir,hsDir[i]))
        hsDir[i] = projHp(hsDir[i],zermat(hpDir),hpDir)

    # project on the HS in the space defined by the HP
    x = projHs(xHp,hsPos,hsDir)

    x = plus(bx(x,basis),xOriginal0)
    testHp(x,hpPos0,hpDir0)
    testHs(x,hsPos0,hsDir0)
    return  smv(x)

# project x0 onto the HP
def projHp(x0,pos,dir):
    pos,dir = filterPosDir(pos,dir)
    if len(dir)==0:
        return x0

    scalMat = prod(trans(dir),dir)
    l = [ dot(dir[j],minus(pos[j],x0)) for j in range(size2(dir)) ]

    lambdas = solveLinear(scalMat,l)

    x = plus(bx(lambdas,dir),x0)
    testHp(x,pos,dir)
    return smv(x)

# projetc x0 onto the HS, solve quadratic programming problem with the interior point alorithm
def projHs(x0,pos,dir):
    pos,dir = filterPosDir(pos,dir)
    if len(dir)==0:
        return x0

    A = trans(dir)
    b = [ dot(dir[j],minus(pos[j],x0)) for j in range(size1(A)) ]
    dx = quadSolve(A,b) # solve quadratic programming problem
    x = plus(x0,dx)
    testHs(x,pos,dir)
    return smv(x)

# build orthogonal basis of the space defined by V (not necessarily independent). Stabilized Gram-Schmidt algorithm (code copy-pasted from Wikipedia, and modified)
def basisGramSchmidt(V):
    n = size1(V)
    k = size2(V)
    U = []

    U.append( scalprod(1/norm2(V[0]),V[0]) )
    for v in V:
        for u in U:
            v = minus( v , scalprod( dot(v,u)/norm2(u)**2 , u ) )
        if norm2(v)>1e-8: # add only if vector is independent
            U.append( scalprod(1/norm2(v),v) )
    return U

def filterPosDir(pos,dir):
    if size1(pos)!=size1(dir) or size2(pos)!=size2(dir):
        disp(pos,dir)
        raise ValueError("wrong size")
    if isempty(dir):
        return pos,dir
    for i in range(len(dir)):
        if norm2(dir[i])==0:
            raise ValueError("direction is zero")
        dir[i] = scalprod(1/norm2(dir[i]),dir[i])
    return (pos,dir)

############################################################################################
## TESTS

def testLinear():
    print(">>test solveLinear()")
    n = 4
    A = [ [sqrt(i+j**2) for j in range(n)] for i in range(n)]
    b = [ sqrt(i) for i in range(n)]
    printmat(A)
    printmat(b)
    x = solveLinear(A,b)
    printmat(x)
    printmat(prod(A,x))
    printmat(minus(prod(A,x),b))
    print(">>end test solveLinear()")

def testProd():
    print("\n>>test prod:")
    n = 4
    A = [zeros(n) for i in range(n)]
    for i in range(n-1):
        A[i+1][i] = 1
    printmat(A)
    B = A
    for i in range(n-1):
        B = prod(B,A)
        printmat(B)
    print(">>end test prod")

def testTrans():
    print("\n>>test trans:")
    n = 4
    A = [zeros(n) for i in range(n)]
    for i in range(n-1):
        A[i+1][i] = 1
    printmat(A)
    printmat(trans(A))
    printmat(A)
    print(">>end test trans")

def testQuadSolve():
    print("\n>>test quadSolve:")
    n = 4
    A = [[1,5,1],[2,7,6],[3,9,8]]
    b = [0,6,2]
    solmatlab = [0.1935, 0.2710,0.3484]
    printmat(A)
    printmat(b)
    c = quadSolve(A,b)
    printmat(c)
    diff = minus(c,solmatlab)
    printmat(diff)
    if max(smv(diff))>1e-3:
        raise ValueError("test failed")
    print(">>end test quadSolve")

def testHp(x,pos,dir):
    n = size1(dir)
    nHp = size2(dir)
    hpDist = [dot(minus(x,pos[i]),dir[i])/n for i in range(nHp)]
    hpTest = [abs(hpDist[i])>1e-2 for i in range(nHp)]
    if sum(hpTest)>0:
        # printmat(hpDist)
        # printmat(hpTest)
        print("projection WARNING> HP not enforced, ",hpDist)
        # raise ValueError("HP not enforced")

def testHs(x,pos,dir):
    nHs = size2(dir)
    if nHs ==0:
        return
    n = size1(dir)
    dist = [dot(minus(x,pos[i]),dir[i])/n for i in range(nHs)]
    test = [dist[i]<-1e-2 for i in range(nHs)]
    if sum(test)>0:
        # dispv(dist)
        # dispv(test)
        print("projection WARNING> HS not enforced, ",dist)
        # raise ValueError("HS not enforced")

def testMatlabHs(): # compare with matlab value
    x0 = [1,2,3]
    pos = [[0,3,4]]
    dir = [[0,1,1]]
    res = projHs(x0,pos,dir)
    printmat(res)
    solmatlab = [1,3,4]

def testMatlabHp(): # compare with matlab value
    x0 = [1,2,3]
    pos = [[0,3,4]]
    dir = [[0,1,1]]
    res = projHp(x0,pos,dir)
    printmat(res)
    solmatlab = [1,3,4]

def testMatlabHpHs(): # compare with matlab value
    x0 = [1,2,3]
    pos = [[0,3,4]]
    dir = [[0,1,1]]
    res = projHpHs(x0,[],[],pos,dir)
    printmat(res)
    solmatlab = [1,3,4]

def printListValues(dataList,printConsole=True):
    ss = ""
    for v in dataList:
        if isinstance(v,str):
            ss += v
            continue

        if isinstance(v,list): # [str{},formatNbString,listvalues]
            strWrap = v[0]
            strFormat = v[1]
            values = v[2]
            if isinstance(values,(int,float)):
                values = [values]
            if len(values)==0:
                continue
            ss += strWrap.format("".join([strFormat.format(l) for l in values]))
            continue
        raise ValueError("wrong type")
    if printConsole:
        print(ss)
    return ss

# if __name__=="__main__":
    # testLinear()
    # testProd()
    # testTrans()
    # testQuadSolve()
    # testMatlabHs()
    # testMatlabHp()
    # testMatlabHpHs()