from display import *
from matrix import *
from vector import *
from sys import maxint
import math

z = [[-500  for x in range(500)] for x in range(500)] 
#def new_z():
#    z = [[-500 for x in range(500)] for x in range(500)] 


def SameSide(p1,p2, a,b):
    cp1 = cross_prod(vect_minus(b,a), vect_minus(p1,a))
    cp2 = cross_prod(vect_minus(b,a), vect_minus(p2,a))
    if (dot_prod(cp1, cp2) >= 0):
        return True
    else:
        return False

def PointInTriangle(p, a,b,c):
    if ((SameSide(p,a, b,c) and SameSide(p,b, a,c)) and SameSide(p,c, a,b)):
        return True
    else:
        return False

def add_polygon( points, x0, y0, z0, x1, y1, z1, x2, y2, z2):
    add_point( points, x0, y0, z0 )
    add_point( points, x1, y1, z1 )
    add_point( points, x2, y2, z2 )

    p0 = [x0, y0, z0]
    p1 = [x1, y1, z1]
    p2 = [x2, y2, z2]
    
    zn = (p0[2] + p1[2] + p2[2]) / 3
   
    def z_buff(a, b, c):
        #print "dfgdsf"
        
        #print "a"
        #print a
        #print int(a)
        #print c
        #print z[int(a)][int(b)]
        if (z[int(a)][int(b)] < c):
            #print "change"
            #print "old"
            #print z[int(a)][int(b)]
            z[int(a)][int(b)] = c
            #print "new"
            #print z[int(a)][int(b)]

    
    for x in range ( int (min ([p0[0], p1[0], p2[0]])), (int (max ([p0[0], p1[0], p2[0]])))):
        for y in range ( int (min ([p0[1], p1[1], p2[1]])), (int (max ([p0[1], p1[1], p2[1]])))):
    #        p = [x, y, zn]
    #        if PointInTriangle(p, p0, p1, p2):
            z_buff(x, y, zn)
                #print "x"
                #print x
                #print "y"
                #print y
            
    z_buff (x0, y0, z0)
    z_buff (x1, y1, z1)
    z_buff (x2, y2, z2)
    
    #print "called"
    #print p0
    #print p1
    #print p2
    #add to z-buffer here???
   
def add_polygon_p(points, p0, p1, p2):
    add_polygon(points, 
                p0[0], p0[1], p0[2],
                p1[0], p1[1], p1[2],
                p2[0], p2[1], p2[2])

    

def draw_polygons(points, screen, color):
    
    def sortaequal(a,b,tol):
        return not abs(a-b)<tol
    def scanlines(p0,p1,p2):
        pts = sorted( (p0,p1,p2), key=lambda pt: pt[1])
        top = pts[0]; mid = pts[1]; bot = pts[2]

        yi = top[1]
        xi0 = top[0]
        xi1 = top[0]

        #YO
        #take into account when top/mid or mid/bot are same y coordinate
        #also, the edge case for drawing torus polygons appears to have
        #    vertices drawn in the wrong order
        dx0  = (bot[0]-top[0])/(bot[1]-top[1]) \
               if sortaequal(bot[1],top[1],0.00001) else 1
        dx1m = (mid[0]-top[0])/(mid[1]-top[1]) \
               if sortaequal(mid[1],top[1],0.00001) else 1
        dx1b = (bot[0]-mid[0])/(bot[1]-mid[1]) \
               if sortaequal(bot[1],mid[1],0.00001) else 1

        while yi < mid[1]:
            debug = False
            if abs(yi-mid[1])<0.1: debug = True
            xi0 += dx0
            xi1 += dx1m
            yi  += 1
            print mid[1]
            print xi0,yi,"!!",xi1,yi
            print "drawing"
            draw_line(screen, xi0,yi, xi1,yi, color)
        while yi < bot[1]:
            xi0 += dx0
            xi1 += dx1b
            yi  += 1
            draw_line(screen, xi0,yi, xi1,yi, color)

    def draw_polygon(p0,p1,p2, c):
        draw_line(screen, p0[0], p0[1], p1[0], p1[1], c)
        draw_line(screen, p1[0], p1[1], p2[0], p2[1], c)
        draw_line(screen, p2[0], p2[1], p0[0], p0[1], c)
        #scanlines(p0,p1,p2)

    view_vect = [0, 0, -1]

    if len( points ) % 3 !=  0:
        print "Bad number of points to draw polygons: not div by 3?" 

    p = 0
    while p < len(points)-1:
        p0 = points[p]
        p1 = points[p+1]
        p2 = points[p+2]
        surf_norm = cross_prod(vect_minus(p1,p0),vect_minus(p2,p0))

        def front(a, b, c):
            #if c >= (z[int(a)][int(b)]):
            #    print 'front'
            #else:
            #    print "back"
            print "front"
            print c
            print z[int(a)][int(b)]
            print c >= (z[int(a)][int(b)])
            return c >= (z[int(a)][int(b)])
        
        red = [255, 0, 0]
        #add z-buffer check here????
        if dot_prod(surf_norm, view_vect) < 0:
            if (front (p0[0], p0[1], p0[2]) and front (p1[0], p1[1], p1[2])) and front (p2[0], p2[1], p2[2]):
                draw_polygon(points[p], points[p+1], points[p+2], red)
                #pass
            else:
                draw_polygon(points[p], points[p+1], points[p+2], color)

        p+=3

    #print "tesT"
    #print z[40][40]
    for i in range (0, 500, 10):
        s = ""
        for k in range (0, 500, 10):
            s = s + "\t " + str(int(z[i][k]))
        print s + "\n"
    #new_z()
        

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

