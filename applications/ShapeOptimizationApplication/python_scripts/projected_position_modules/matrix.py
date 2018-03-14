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

##########################
## FUNCTIONS TO MANIPULATE MATRICES


def copy(x):
    if isvector(x):
        return x[:]
    return [[x[j][i] for i in range(size1(x))] for j in range(size2(x))]

def norm2(x):
    x = smv(x)
    return sqrt(sum(v**2 for v in x))

def norm2safe(x):
    nf = norm2(x)
    if nf==0:
        return 1
    return nf

# for vectors
def scalprod(t,x):
    x = smv(x)
    return [t*v for v in x]

def norminf3d(x):
    x = smv(x)
    ns = []
    for i in range(int(round(len(x)/3))):
        ns.append(norm2(x[3*i:3*(i+1)]))
    return max(ns)

def zeros(n):
    return [0 for i in range(n)]

def ones(n):
    return [1 for i in range(n)]

def isvector(x):
    if len(x)==0 or isinstance(x[0],(float,int)):
        return True
    else:
        return False

def ismatrix(x):
    if len(x)==0:
        return False
    if len(x[0])==0:
        return True
    if isinstance(x[0][0],(float,int)):
        return True
    else:
        return False

def isempty(x):
    if len(x)==0 or len(x[0])==0:
        return True
    return False

#safe vector to matrix
def svm(x):
    if isvector(x):
        return [x]
    return x

# safe matrix to vector
def smv(x):
    if isvector(x):
        return x
    if size2(x)>1:
        raise ValueError("wrong size")
    return x[0]

# for matrices or vector
def horzcat(a,b):
    a = svm(a)
    b = svm(b)
    if isempty(b):
        return a
    if isempty(a):
        return b
    return a+b

def vertcat(a,b):
    a = svm(a)
    b = svm(b)
    if isempty(b):
        return a
    if isempty(a):
        return b
    if size2(a)!=size2(b):
        raise ValueError("wrong size")
    c = []
    for i in range(size2(a)):
        c.append(a[i]+b[i])
    return c

def linspace(u,v,n):
    return [u+(v-u)*i/n for i in range(n)]

def minus(a,y):
    a = svm(a)
    y = smv(y)
    if isempty(a):
        return []
    if size1(a)!=len(y):
        printmat("error minus()",a,y)
        raise ValueError("wrong size")
    return [ [ a[j][i]-y[i] for i in range(size1(a))] for j in range(size2(a))]

def plus(a,y):
    a = svm(a)
    y = smv(y)
    if isempty(a):
        return []
    if size1(a)!=len(y):
        raise ValueError("wrong size")
    return [ [ a[j][i]+y[i] for i in range(size1(a))] for j in range(size2(a))]

# for vectors
def dot(x,y):
    x = smv(x)
    y = smv(y)
    return sum(u*v for u,v in zip(x,y))

# for matrices
def prod(a,b):
    a = svm(a)
    b = svm(b)
    if isempty(a) or isempty(b):
        return []
    if size2(a)!=size1(b):
        print("prod error:")
        disp(a,b)
        raise ValueError("wrong size")
    return [ [ sum(a[k][i]*b[j][k] for k in range(size2(a))) for i in range(size1(a)) ] for j in range(size2(b)) ]

# for matrices
def trans(a):
    if isempty(a):
        return []
    return [ [ a[j][i] for j in range(size2(a)) ] for i in range(size1(a))  ]

def eltprod(x,y):
    x = smv(x)
    y = smv(y)
    if len(x)!=len(y):
        raise ValueError("wrong size")
    return [x[i]*y[i] for i in range(len(x))]

def size1(a):
    if len(a)==0:
        return 0
    return len(a[0])

def size2(a):
    return len(a)

def vertslice(a,r):
    if isempty(a):
        return []
    return [ [ a[j][k] for k in r ] for j in range(size2(a)) ]

# get x coordinates in the basis
def xb(x,basis):
    return prod(trans(basis),x)

# get original x from coordinate in baiss
def bx(x,basis):
    return prod(basis,x)

def dispv(*xs):
    disp(*xs,printValue=True)

def disp(*xs,printValue=False):
    print("")
    for x in xs:
        if isinstance(x,str):
            print(x)
        elif isinstance(x,(float,int)):
            print("Number: ",round(x,2))
        elif len(x)==0:
            print("Empty list")
        elif isvector(x):
            print("Vector: ",len(x))
            if printValue:
                line = ""
                for i in range(len(x)):
                    line += str(round(x[i],2)) + "\t"
                print(line)
        elif ismatrix(x):
            print("Matrix: (",size1(x),",",size2(x),")")
            if printValue:
                for i in range(size1(x)):
                    line = ""
                    for j in range(size2(x)):
                        line += str(round(x[j][i],2)) + "\t"
                    print(line)
        else:
            print("unrecognized: ",type(x))
            if printValue:
                print(x)

def printmat(*As):
    print("")
    for A in As:
        if isinstance(A,str):
            print(A)
        elif isvector(A):
            print("Vector:",len(A))
            print(A)
        else:
            print("Matrix: (",size1(A),",",size2(A),")")
            for i in range(size1(A)):
                line = ""
                for j in range(size2(A)):
                    line += str(round(A[j][i],2)) + "\t"
                print(line)

def solveLinear(A,b):
    b = smv(b)
    A = trans(A)
    n = len(A)
    for i in range(n):
        A[i].append(b[i])
    for i in range(0, n):
        # Search for maximum in this column
        maxEl = abs(A[i][i])
        maxRow = i
        for k in range(i+1, n):
            if abs(A[k][i]) > maxEl:
                maxEl = abs(A[k][i])
                maxRow = k

        # Swap maximum row with current row (column by column)
        for k in range(i, n+1):
            tmp = A[maxRow][k]
            A[maxRow][k] = A[i][k]
            A[i][k] = tmp

        # Make all rows below this one 0 in current column
        for k in range(i+1, n):
            c = -A[k][i]/A[i][i]
            for j in range(i, n+1):
                if i == j:
                    A[k][j] = 0
                else:
                    A[k][j] += c * A[i][j]

    # Solve equation Ax=b for an upper triangular matrix A
    x = [0 for i in range(n)]
    for i in range(n-1, -1, -1):
        x[i] = A[i][n]/A[i][i]
        for k in range(i-1, -1, -1):
            A[k][n] -= A[k][i] * x[i]
    return x

# Interior point algorithm
# Nocedal & Wright, Numerical Optimization, p.506, Springer
# solve min x'*x   with A*x>=b, x is of size n and A of size m*n
def quadSolve(A,b):
    m = size1(A)
    n = size2(A)

    def rdrp(x,y,l,A,b):
        rd = minus(x,prod(trans(A),l))
        rp = minus(minus(prod(A,x),y),b)
        return (rd,rp)

    def gradResidu(y,l,A):
        m = size1(A)
        n = size2(A)
        grad = [zeros(n+2*m) for i in range(n+2*m)]
        for i in range(n):
            grad[i][i] = 1
        for i in range(n):
            for j in range(m):
                grad[n+m+j][i] = -A[i][j]
        for i in range(m):
            for j in range(n):
                grad[j][n+i] = A[j][i]
        for i in range(m):
            grad[n+i][n+i] = -1
        for i in range(m):
            grad[n+i][n+m+i] = l[i]
            grad[n+m+i][n+m+i] = y[i]
        return grad

    #init
    x = zeros(n)
    y = ones(m) # slack variables
    l = ones(m) # lagrange multipliers

    rd, rp = rdrp(x,y,l,A,b)
    k = 0
    while norm2(vertcat(rd,rp))>1e-6:
        # solve affine delta
        gradRes = gradResidu(y,l,A)
        deltaXYLAff = zeros(n+2*m)
        rhs = scalprod(-1, vertcat(rd,vertcat(rp,eltprod(y,l))))
        deltaXYLAff = solveLinear(gradRes,rhs)

        dxAff = deltaXYLAff[:n]
        dyAff = deltaXYLAff[n:n+m]
        dlAff = deltaXYLAff[n+m:]

        mu = dot(y,l)/m
        alphaS = 1
        alphaZ = 1
        for i in range(m):
            if dyAff[i]<0:
                alphaS = min(alphaS, abs(y[i]/dyAff[i]))
            if dlAff[i]<0:
                alphaZ = min(alphaZ, abs(l[i]/dlAff[i]))
        alphaAff = 0.8*min(alphaS,alphaZ)

        muAff = dot(plus(y,scalprod(alphaAff,dyAff)),plus(l,scalprod(alphaAff,dlAff))) / m
        sigma = (muAff/mu)**3

        # solve for delta
        rhs = scalprod(-1, vertcat(rd,vertcat(rp, plus( plus(eltprod(y,l),eltprod(dyAff,dlAff)) , scalprod(-sigma*mu,ones(m)) ))))
        deltaXYL = solveLinear(gradRes,rhs)

        dx = deltaXYL[:n]
        dy = deltaXYL[n:n+m]
        dl = deltaXYL[n+m:]

        alphaS = 1
        alphaZ = 1
        for i in range(m):
            if dy[i]<0:
                alphaS = min(alphaS, abs(y[i]/dy[i]))
            if dl[i]<0:
                alphaZ = min(alphaZ, abs(l[i]/dl[i]))
        alpha = 0.8*min(alphaS,alphaZ)

        # apply delta
        x = smv(plus(x,scalprod(alpha,dx)))
        y = smv(plus(y,scalprod(alpha,dy)))
        l = smv(plus(l,scalprod(alpha,dl)))
        rd,rp = rdrp(x,y,l,A,b)
        k = k+1
        if k>1000:
            raise ValueError("too many iterations")
    return x