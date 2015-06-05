from display import *
from matrix import *
from vector import *
import math


def add_polygon( points, x0, y0, z0, x1, y1, z1, x2, y2, z2 ):
    add_point( points, x0, y0, z0 )
    add_point( points, x1, y1, z1 )
    add_point( points, x2, y2, z2 )
    
def add_polygon_p(points, p0, p1, p2):
    add_polygon(points, 
                p0[0], p0[1], p0[2],
                p1[0], p1[1], p1[2],
                p2[0], p2[1], p2[2])

def draw_polygons( points, screen, color ):
    def draw_polygon(p0,p1,p2):
        draw_line(screen, p0[0], p0[1], p1[0], p1[1], color)
        draw_line(screen, p1[0], p1[1], p2[0], p2[1], color)
        draw_line(screen, p2[0], p2[1], p0[0], p0[1], color)

    view_vect = [0, 0, -1]

    if len( points ) % 3 !=  0:
        print "Bad number of points to draw polygons: not div by 3?" 

    p = 0
    while p < len(points)-1:
        p0 = points[p]
        p1 = points[p+1]
        p2 = points[p+2]
        surf_norm = cross_prod(vect_minus(p1,p0),vect_minus(p2,p0))
                                        
        if dot_prod(surf_norm, view_vect) < 0:
            draw_polygon(points[p], points[p+1], points[p+2])
        p+= 3

def add_prism(points,x,y,z,w,h,d):
    def add_f(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3):
        add_polygon(points, x0,y0,z0,  x1,y1,z1,  x2,y2,z2)
        add_polygon(points, x0,y0,z0,  x2,y2,z2,  x3,y3,z3)
    add_f(x  ,y  ,z  ,  x+w,y  ,z  ,  x+w,y+h,z  ,  x  ,y+h,z  )#F
    add_f(x  ,y  ,z+d,  x  ,y  ,z  ,  x  ,y+h,z  ,  x  ,y+h,z+d)#L
    add_f(x+w,y  ,z+d,  x  ,y  ,z+d,  x  ,y+h,z+d,  x+w,y+h,z+d)#B
    add_f(x+w,y  ,z  ,  x+w,y  ,z+d,  x+w,y+h,z+d,  x+w,y+h,z  )#R
    add_f(x  ,y  ,z+d,  x+w,y  ,z+d,  x+w,y  ,z  ,  x  ,y  ,z  )#U
    add_f(x  ,y+h,z+d  ,x  ,y+h,z  ,  x+w,y+h,z  ,  x+w,y+h,z+d)#D

def add_sphere(points,cx,cy,cz,r,step):
    spts = []

    #generate points of sphere
    i = 0.0
    j = 0.0
    while i<=step:
        ti = i/step
        theta = math.pi * ti
        while j<=step:
            tj = j/step
            phi = 2 * math.pi * tj
            x = r*math.cos(theta) + cx
            y = r*math.sin(theta)*math.cos(phi) + cy
            z = r*math.sin(theta)*math.sin(phi) + cz
            spts += [[x,y,z]]
            j+=1
        i+=1
        j=1.0

    #connect points
    i = 0
    while i < len(spts)-step:
        a = (i)        % len(spts)
        b = (i+step)   % len(spts)
        if i%step<1:
            #special case for when switching to next layer
            c = (i+1)      % len(spts)
            d = (i-step+1) % len(spts)
        else:
            c = (i+step+1) % len(spts)
            d = (i+1)      % len(spts)
        add_polygon_p(points, spts[a], spts[b], spts[c])
        if i>step:
            add_polygon_p(points, spts[a], spts[c], spts[d])
        i+=1

def add_torus(points,cx,cy,cz,r1,r2,step):
    spts = []

    #generate points
    i = 1.0
    j = 1.0
    while i<=step:
        ti = i/step
        theta = 2 * math.pi * ti
        while j<=step:
            tj = j/step
            phi = 2 * math.pi * tj
            x = (r1*math.cos(theta) + r2) * math.cos(phi) + cx
            y = r1*math.sin(theta) + cy
            z = -1 * (r1*math.cos(theta) + r2) * math.sin(phi) + cz
            spts += [[x,y,z]]
            j+=1
        i+=1
        j=1.0
        
    #connect points
    i = 0
    while i < len(spts)-step:
        a = (i)        % len(spts)
        b = (i+step)   % len(spts)
        c = (i+1)      % len(spts)
        d = (i-step+1) % len(spts)
        add_polygon_p(points, spts[a], spts[b], spts[c])
        add_polygon_p(points, spts[a], spts[c], spts[d])
        i+=1

def add_parametric(points,param_x,param_y,param_z,step):
    x0 = param_x(0.0)
    y0 = param_y(0.0)
    z0 = param_z(0.0)
    i = 1.0
    while i<=step:
        t = i/step
        x1 = param_x(t)
        y1 = param_y(t)
        z1 = param_z(t)
        add_edge(points, x0, y0, z0, x1, y1, z1)
        x0 = x1
        y0 = y1
        z0 = z1
        i+=1
    
def add_circle( points, cx, cy, cz, r, step ):
    def xt(t):
        return cx + r*math.cos(2*math.pi*t)
    def yt(t):
        return cy + r*math.sin(2*math.pi*t)
    def zt(t):
        return 0
    add_parametric(points,xt,yt,zt,step)

def add_curve( points, x0, y0, x1, y1, x2, y2, x3, y3, step, curve_type ):
    if curve_type == 0 or curve_type == "hermite":
        hmat = make_hermite()
        dx0 = x1 - x0
        dy0 = y1 - y0
        dx2 = x3 - x2
        dy2 = y3 - y2
        xc = generate_curve_coefs(x0,x2,dx0,dx2,hmat)[0]
        yc = generate_curve_coefs(y0,y2,dy0,dy2,hmat)[0]
    elif curve_type == 1 or curve_type == "bezier":
        bmat = make_bezier()
        xc = generate_curve_coefs(x0,x1,x2,x3,bmat)[0]
        yc = generate_curve_coefs(y0,y1,y2,y3,bmat)[0]
    else:
        print "Invalid curve type somehow?"

    def xt(t):
        t2 = t*t
        t3 = t2*t
        return xc[0]*t3 + xc[1]*t2 + xc[2]*t + xc[3]
    def yt(t):
        t2 = t*t
        t3 = t2*t
        return yc[0]*t3 + yc[1]*t2 + yc[2]*t + yc[3]
    def zt(t):
        return 0
    add_parametric(points,xt,yt,zt,step)

def draw_lines( matrix, screen, color ):
    if len( matrix ) < 2:
        print "Need at least 2 points to draw a line"
        
    p = 0
    while p < len( matrix ) - 1:
        draw_line( screen, matrix[p][0], matrix[p][1],
                   matrix[p+1][0], matrix[p+1][1], color )
        p+= 2

def add_edge( matrix, x0, y0, z0, x1, y1, z1 ):
    add_point( matrix, x0, y0, z0 )
    add_point( matrix, x1, y1, z1 )

def add_point( matrix, x, y, z=0 ):
    matrix.append( [x, y, z, 1] )


def draw_line( screen, x0, y0, x1, y1, color ):
    dx = x1 - x0
    dy = y1 - y0
    if dx + dy < 0:
        dx = 0 - dx
        dy = 0 - dy
        tmp = x0
        x0 = x1
        x1 = tmp
        tmp = y0
        y0 = y1
        y1 = tmp
    
    if dx == 0:
        y = y0
        while y <= y1:
            plot(screen, color,  x0, y)
            y = y + 1
    elif dy == 0:
        x = x0
        while x <= x1:
            plot(screen, color, x, y0)
            x = x + 1
    elif dy < 0:
        d = 0
        x = x0
        y = y0
        while x <= x1:
            plot(screen, color, x, y)
            if d > 0:
                y = y - 1
                d = d - dx
            x = x + 1
            d = d - dy
    elif dx < 0:
        d = 0
        x = x0
        y = y0
        while y <= y1:
            plot(screen, color, x, y)
            if d > 0:
                x = x - 1
                d = d - dy
            y = y + 1
            d = d - dx
    elif dx > dy:
        d = 0
        x = x0
        y = y0
        while x <= x1:
            plot(screen, color, x, y)
            if d > 0:
                y = y + 1
                d = d - dx
            x = x + 1
            d = d + dy
    else:
        d = 0
        x = x0
        y = y0
        while y <= y1:
            plot(screen, color, x, y)
            if d > 0:
                x = x + 1
                d = d - dy
            y = y + 1
            d = d + dx

