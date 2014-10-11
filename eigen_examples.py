"""
eigen_example_2.py

Includes code for functions that do basic vector and
matrix arithmetic.  Most of these functions support
the matrix multiplication and norm calculation needed for
iterating to an eigenvector.  The iteration should converge
if the matrix has a dominant eigenvalue.

"""


def rows(mat):
    "return number of rows"
    return(len(mat))

def cols(mat):
    "return number of cols"
    return(len(mat[0]))
 
def zero(m,n):
    "Create zero matrix"
    new_mat = [[0 for col in range(n)] for row in range(m)]
    return new_mat
 
def transpose(mat):
    "return transpose of mat"
    new_mat = zero(cols(mat),rows(mat))
    for row in range(rows(mat)):
        for col in range(cols(mat)):
            new_mat[col][row] = mat[row][col]
    return(new_mat)

def dot(A,B):
    "vector dot product"
    if len(A) != len(B):
        print("dot: list lengths do not match")
        return()
    dot=0
    for i in range(len(A)):
        dot = dot + A[i]*B[i]
    return(dot)

def getCol(mat, col):
    "return column col from matrix mat"
    return([r[col] for r in mat])

def getRow(mat, row):
    "return row row from matrix mat"
    return(mat[row])

def matMult(mat1,mat2):
    "multiply two matrices"
    if cols(mat1) != rows(mat2):
        print("multiply: mismatched matrices")
        return()
    prod = zero(rows(mat1),cols(mat2))
    for row in range(rows(mat1)):
        for col in range(cols(mat2)):
            prod[row][col] = dot(mat1[row],getCol(mat2,col))
    return(prod)

def vectorQ(V):
    "mild test to see if V is a vector"
    if type(V) != type([1]):
        return(False)
    if type(V[0]) == type([1]):
        return(False)
    return(True)

def scalarMult(a,mat):
    "multiply a scalar times a matrix"
    if vectorQ(mat):
        return([a*m for m in mat])
    for row in range(rows(mat)):
        for col in range(cols(mat)):
            mat[row][col] = a*mat[row][col]
    return(mat)

def addVectors(A,B):
    "add two vectors"
    if len(A) != len(B):
        print("addVectors: different lengths")
        return()
    return([A[i]+B[i] for i in range(len(A))])


def show(mat):
    "Print out matrix"
    for row in mat:
        print(row)

### vectors vs rowVectors and colVectors
### the latter are matrices

def vec2rowVec(vec):
    "[a,b,c] -> [[a,b,c]]"
    return([vec])

def vec2colVec(vec):
    return(transpose(vec2rowVec(vec)))

def colVec2vec(mat):
    rowVec = transpose(mat)
    return(rowVec[0])


def norm_infinity(vec):
    """
    Given an augmented matrix it returns a list.  The [0]
    element is the infinity norm of vec, norm_vec[1]=p 
    the index of the first component whose absolute value
    is the norm, and norm_vec[2] is vec[p].
    
    """
    n = len(vec)
    x = [vec[i] for i in range(n)]
    norm = abs(x[0])
    norm_index = 0

    for i in range(n):
        if  norm < abs(x[i]):
            norm = abs(x[i])
            norm_index = i
    norm_vec = [norm,norm_index,vec[norm_index]]
    return(norm_vec)


def x_k(A,x,k):
    """
    Given an n by n matrix A, a nonzero n vector x, and a positve
    integer k it returns A^k x.  The idea is to show that
    the iterates become parallel, but blow up or converge to
    the zero vector in typical cases.
    
    """
    n = len(x)
    x_iterate = [x[i] for i in range(n)]  # a copy of x
    x_col_vec = vec2colVec(x_iterate)
    
    for i in range(k):
        x_col_vec = matMult(A,x_col_vec)
       
    return(x_col_vec)


def y_k(A,x,k):
    """
    Given an n by n matrix A, a nonzero n vector x, and a positve
    integer k it returns the kth iterate of an estimate for
    the dominant eigenvalue of A and a corresponding unit
    eigenvector.  This is an implementation of the Power Method.
    
    """
    eigenval_eigenvec = [0, [0,0,0]]

    norm_vec = norm_infinity(x)
    norm = norm_vec[0]
    p = norm_vec[1]
    x_iterate = [x[i] for i in range(len(x))]  # a copy of x
    y_iterate = scalarMult(1.0/norm, x_iterate)
    y = y_iterate

    for i in range(k-2):              
        y_col_vec = vec2colVec(y)
        y_iterate = matMult(A,y_col_vec)
        y = colVec2vec(y_iterate)
        norm_vec = norm_infinity(y)
        norm = norm_vec[0]
        p = norm_vec[1]
        y = scalarMult(1.0/norm, y)
        
    norm_vec = norm_infinity(y)
    y_p = norm_vec[2]
    
    y_col_vec = vec2colVec(y)
    y_iterate = matMult(A,y_col_vec)
    y = colVec2vec(y_iterate)
    norm_vec = norm_infinity(y)
    eigenval_eigenvec[0] = y[p]/y_p
    norm = norm_vec[0]
    y = scalarMult(1.0/norm, y)
    eigenval_eigenvec[1] = y

    return(eigenval_eigenvec)


### Some Testing data begins here.

A = [[2,1],
     [1,2]]


X0A = [1,0]


B = [[ 5,-2, 0],
    [-2, 4,-1],
    [ 0,-1, 1]]


X0B = [1,1,1]





