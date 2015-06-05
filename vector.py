def cross_prod(a,b):
    return [a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]]

def dot_prod(a,b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def vect_minus(a,b):
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]
