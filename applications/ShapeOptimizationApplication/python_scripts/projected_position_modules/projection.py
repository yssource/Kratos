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

from math import sqrt
from projected_position_modules.matrix import *


def plusBall(radiusBall,directionBall,lConstraint,eConstraint,lEquality,eEquality):
    hsPos = []; hsDir = []; hpPos = []; hpDir = []
    for i in range(len(eConstraint)):
        if norm2(eConstraint[i])>0 and lConstraint[i]>-radiusBall:
            hsPos.append(scalprod(lConstraint[i],eConstraint[i]))
            hsDir.append(scalprod(1/norm2(eConstraint[i]),eConstraint[i]))
    for i in range(len(eEquality)):
        if norm2(eEquality[i])>0:
            hpPos.append(scalprod(lEquality[i],eEquality[i]))
            hpDir.append(scalprod(1/norm2(eEquality[i]),eEquality[i]))
    directionBall = scalprod(1/norminf3d(directionBall) , directionBall)

    x,lBall,sInit,normInfX = projBallHpHs(radiusBall,directionBall,hpPos,hpDir,hsPos,hsDir)

    lXConstraint = zeros(len(lConstraint))
    for i in range(len(lConstraint)):
        lXConstraint[i] = max(lConstraint[i],dot(x,scalprod(1/norm2safe(eConstraint[i])**2,eConstraint[i])))
    x = smv(x)
    return (x,lBall,sInit,normInfX,lXConstraint)


def projBallHpHs(radiusBall,directionBall,hpPos,hpDir,hsPos,hsDir):
    n = len(directionBall)

    # check length init
    xHpHsInit = projHpHs(zeros(n),hpPos,hpDir,hsPos,hsDir)
    sInit = norminf3d(xHpHsInit)
    if sInit>radiusBall: # if the projection of 0 onto the hs and hp is outside of the ball, no solution exist. return
        lBall = dot(xHpHsInit,scalprod(1/norm2(directionBall)**2,directionBall))
        print("projection.projBallHpHs> sInit={:.2f}".format(sInit/radiusBall))
        return (xHpHsInit,lBall,sInit,sInit)

    nl = []
    lMin = -0.2
    lMax = 1.2
    nInterval = 5
    nTry = 0
    while nTry<10:
        for lBall in linspace(lMin,lMax,nInterval):
            x = projHpHs(zeros(n),hpPos,hpDir,horzcat(hsPos,scalprod(lBall*radiusBall,directionBall)),horzcat(hsDir,directionBall))
            normInfX = norminf3d(x)/radiusBall
            nl.append((lBall,normInfX))
            print("projection.projBallHpHs> norminf={:.2f}".format(normInfX),end="")
            if normInfX>=0.95 and normInfX<=1:
                print("\n",end="")
                return (x,lBall*radiusBall,sInit,normInfX*radiusBall)
            else:
                print("\r",end="")
                if normInfX>1:
                    break
        lMin,lMax = getIntervalBall(nl)
        nTry += 1
    raise ValueError("projball not found")

def getIntervalBall(nl):
    listLower = [(l,n) for l,n in nl if n<0.95]
    listUpper = [(l,n) for l,n in nl if n>1]
    lMin,_ = max(listLower, key = lambda t: t[1])
    if len(listUpper)==0:
        lMax,_ = max(nl, key = lambda t: t[0])
        lMax = 2*lMax
    else:
        lMax,_ = min(listUpper, key = lambda t: t[1])
    return (lMin,lMax)

# def projBallHpHs(radiusBall,directionBall,hpPos,hpDir,hsPos,hsDir):
#     # disp("projball",radiusBall,directionBall,hpPos,hpDir,hsPos,hsDir,norminf3d(directionBall),norm2(hsDir[0]))

#     n = len(directionBall)

#     # check length init
#     # disp("hphsinit")
#     xHpHsInit = projHpHs(zeros(n),hpPos,hpDir,hsPos,hsDir)
#     sInit = norminf3d(xHpHsInit)
#     if sInit>radiusBall: # if the projection of 0 onto the hs and hp is outside of the ball, no solution exist. return
#         lBall = dot(xHpHsInit,scalprod(1/norm2(directionBall)**2,directionBall))
#         print("projection.projBallHpHs> sInit={:.2f}".format(sInit/radiusBall))
#         return (xHpHsInit,lBall,sInit,sInit)

#     # search for x with 0.95 <= norminf(x)/radiusBall <= 1
#     nVal = 100
#     nf = zeros(nVal)
#     # disp("hphs split")
#     x = projHpHs(zeros(n),hpPos,hpDir,horzcat(hsPos,zeros(n)),horzcat(hsDir,directionBall))
#     if norminf3d(x)>radiusBall:
#         lBallRange = linspace(-radiusBall,0,nVal)
#     else:
#         lBallRange = linspace(0,2*radiusBall,nVal)

#     for lBall in lBallRange:
#         x = projHpHs(zeros(n),hpPos,hpDir,horzcat(hsPos,scalprod(lBall,directionBall)),horzcat(hsDir,directionBall))
#         normInfX = norminf3d(x)

#         if normInfX>0.95*radiusBall: # && nf(k)<=radiusBall
#             print("projection.projBallHpHs> norminf={:.2f}".format(normInfX/radiusBall))
#             return (x,lBall,sInit,normInfX)
#         else:
#             print("projection.projBallHpHs> norminf={:.2f}".format(normInfX/radiusBall), end="\r")
#     raise ValueError("projball not found")


# project orthogonally x on the hps and hss
def projHpHs(xOriginal0,hpPos0,hpDir0,hsPos0,hsDir0): # hss (point go through, direction) and hps (point go through,direction)
    hpPos0, hpDir0 = filterPosDir(hpPos0,hpDir0)
    hsPos0, hsDir0 = filterPosDir(hsPos0,hsDir0)

    nHp = len(hpDir0)
    nHs = len(hsDir0)
    iHp = range(nHp)
    iHs = range(nHp,nHp+nHs)
    if nHp==0 and nHs==0:
        return xOriginal0

    # disp("hphs0",hpPos0,hpDir0,hsPos0,hsDir0)

    # convert to reduced coordinates
    basis = basisGramSchmidt(horzcat(hpDir0,hsDir0))

    hpPos = xb(minus(hpPos0,xOriginal0),basis)
    hpDir = xb(hpDir0,basis)
    hsPos = xb(minus(hsPos0,xOriginal0),basis)
    hsDir = xb(hsDir0,basis)

    # project on the hps
    xHp = projHp(zeros(nHp),vertslice(hpPos,iHp),vertslice(hpDir,iHp))

    # conversion of HS to the reduced coordinates:  project HSposition on the other HPs, along its own HP
    for i in range(nHs):
        hsPos[i] = projHp(hsPos[i],horzcat(hpPos,hsPos[i]),horzcat(hpDir,hsDir[i]))

    # project on the hss in the space defined by the hps
    xHs = projHs(zeros(nHs),vertslice(hsPos,iHs),vertslice(hsDir,iHs))

    x = vertcat(xHp,xHs)
    x = plus(bx(x,basis),xOriginal0)
    testHp(x,hpPos0,hpDir0)
    testHs(x,hsPos0,hsDir0)
    return  smv(x)

# project x0 onto the hps
def projHp(x0,pos,dir):
    if len(dir)==0:
        return x0

    scalMat = prod(trans(dir),dir)
    l = [ dot(dir[j],minus(pos[j],x0)) for j in range(size2(dir)) ]

    lambdas = solveLinear(scalMat,l)

    x = plus(bx(lambdas,dir),x0)
    testHp(x,pos,dir)
    return smv(x)

# projetc x0 onto the hss, solve quadratic programming problem with the interior point alorithm
def projHs(x0,pos,dir):
    if len(dir)==0:
        return x0

    A = trans(dir)
    b = [ dot(dir[j],minus(pos[j],x0)) for j in range(size1(A)) ]
    dx = quadSolve(A,b) # solve quadratic programming problem
    x = plus(x0,dx)
    testHs(x,pos,dir)
    return smv(x)

# build orthogonal basis of the space defined by V (not necessarily independent). Stabilized Gram-Schmidt algorithm (code copy-pasted from Wikipedia, and modified)
# 'is' contains the indices of the vector of V that form an independent family
def basisGramSchmidt(V):
    n = size1(V)
    k = size2(V)
    U = [zeros(n) for j in range(k)]
    U[0] = scalprod(1/norm2(V[0]),V[0])
    for i in range(1,k):
        U[i] = V[i]
        for j in range(i-1):
            U[i] = minus( U[i] , scalprod( dot(U[i],U[j])/norm2(U[j])**2 , U[j] ) )
            U[i] = scalprod(1/norm2(U[i]),U[i])
    return U

def filterPosDir(pos,dir):
    if size1(pos)!=size1(dir) or size2(pos)!=size2(dir):
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


if __name__=="__main__":
    # testLinear()
    # testProd()
    # testTrans()
    # testQuadSolve()
    # testMatlabHs()
    # testMatlabHp()
    testMatlabHpHs()