import math

def make_bezier():
    return [ [-1, 3,-3, 1],
             [ 3,-6, 3, 0],
             [-3, 3, 0, 0],
             [ 1, 0, 0, 0] ]

def make_hermite():
    return [ [ 2,-3, 0, 1],
             [-2, 3, 0, 0],
             [ 1,-2, 1, 0],
             [ 1,-1, 0, 0] ]

def generate_curve_coefs( p1, p2, p3, p4, t ):
    pmat = [ [p1,p2,p3,p4] ]
    return matrix_mult(t,pmat)

def make_translate( x, y, z ):
    rmat = new_matrix(4,4)
    rmat = ident(rmat)
    rmat[3][0] = x
    rmat[3][1] = y
    rmat[3][2] = z
    return rmat

def make_scale( x, y, z ):
    rmat = new_matrix(4,4)
    rmat[0][0] = x
    rmat[1][1] = y
    rmat[2][2] = z
    rmat[3][3] = 1
    return rmat
    
def make_rotX( theta ):
    rmat = new_matrix(4,4)
    rad = math.radians(theta)
    rmat[1][1] = math.cos(rad)
    rmat[2][2] = math.cos(rad)
    rmat[1][2] = math.sin(rad)
    rmat[2][1] = -math.sin(rad)
    rmat[0][0] = 1
    rmat[3][3] = 1
    return rmat
    
def make_rotY( theta ):
    rmat = new_matrix(4,4)
    rad = math.radians(theta)
    rmat[0][0] = math.cos(rad)
    rmat[2][2] = math.cos(rad)
    rmat[0][2] = math.sin(rad)
    rmat[2][0] = -math.sin(rad)
    rmat[1][1] = 1
    rmat[3][3] = 1
    return rmat
    
def make_rotZ( theta ):
    rmat = new_matrix(4,4)
    rad = math.radians(theta)
    rmat[0][0] = math.cos(rad)
    rmat[1][1] = math.cos(rad)
    rmat[0][1] = math.sin(rad)
    rmat[1][0] = -math.sin(rad)
    rmat[2][2] = 1
    rmat[3][3] = 1
    return rmat

def new_matrix(rows = 4, cols = 4):
    m = []
    for c in range( cols ):
        m.append( [] )
        for r in range( rows ):
            m[c].append( 0 )
    return m

def print_matrix( matrix ):
    s = ''
    for r in range( len( matrix[0] ) ):
        for c in range( len(matrix) ):
            s+= str(matrix[c][r]) + ' '
        s+= '\n'
    print s

def print_matrix_vert( matrix ):
    for c in matrix:
        print c

def ident( matrix ):
    idmat = new_matrix( len(matrix), len(matrix) )
    for i in range( len( idmat ) ):
        idmat[i][i] = 1
    return idmat

def matrix_copy(matrix):
    if type(matrix) is list:
        r = []
        for e in matrix:
            r.append(matrix_copy(e))
        return r
    else: 
        return matrix

def scalar_mult( matrix, x ):
    for c in range( len(matrix) ):
        for r in range( len( matrix[0] ) ):
            matrix[c][r] *= x

def matrix_mult( m1, m2 ):
    rmat = new_matrix(len(m1[0]),len(m2))
    for c in range( len(m2) ):
        for r in range( len(m1[0]) ):
            cell = 0
            for i in range( len(m1) ):
                cell += m1[i][r] * m2[c][i]
            rmat[c][r] = cell
    return rmat
