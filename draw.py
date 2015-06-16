from display import *
from matrix import *
from vector import *
from sys import maxint
import math
import random

zbuffer = [[-5000 for x in range(500)] for x in range(500)] 
def clear_zbuffer():
    global zbuffer
    zbuffer = [[-5000 for x in range(500)] for x in range(500)] 

def add_polygon( points, x0, y0, z0, x1, y1, z1, x2, y2, z2):
    add_point( points, x0, y0, z0 )
    add_point( points, x1, y1, z1 )
    add_point( points, x2, y2, z2 )
   
def add_polygon_p(points, p0, p1, p2):
    add_polygon(points, 
                p0[0], p0[1], p0[2],
                p1[0], p1[1], p1[2],
                p2[0], p2[1], p2[2])

def draw_polygons(points, screen, env):
    color = [100,100,100] #placeholder
    def sortaequal(a,b,tol):
        return abs(a-b)<tol
    def light(x,y,z,ka,kd): #placeholder
        colortmp = []
        ia = env["ambient"]
        for i in range(3):
            ambient = ka[i] * ia[i]
            for light in env["lights"]:
                vector_l = vect_minus(light[3:],p0)
                diffuse = kd[i] * light[i] * \
                          max(0, dot_prod(normalize(surf_norm), 
                                          normalize(vector_l)))
                #print kd[i], light[i], dot_prod(normalize(surf_norm), normalize(vector_l))
            #specular = ks[i] * (dot_prod(specular_r_vector, vect_minus([0,0,-1],center))) ** sn
            
            colortmp.append(int(round(ambient+diffuse)))
        print colortmp
        return colortmp
            
    def scanlines(p0,p1,p2):
        if env["shading_mode"] == "flat":
            colortmp = light(p0[0],p0[1],p0[2],
                             (0.8,0.8,0.8), #k-ambient
                             (0.8,0.8,0.8)  #k-diffuse
                            )
        else:
            colortmp = random.sample(xrange(255),3)
    
        for i in range(3):
            p0[i] = math.floor(p0[i])
            p1[i] = math.floor(p1[i])
            p2[i] = math.floor(p2[i])
        
        pts = sorted( (p0,p1,p2), key=lambda pt: pt[1])
        top = pts[0]; mid = pts[1]; bot = pts[2]

        yi = top[1]
        x0 = top[0]
        x1 = top[0]
        z0 = top[2]
        z1 = top[2]

        #if sortaequal(bot[1],top[1], 0.0001):
        if bot[1] == top[1]:
            dx0 = 0
            dz0 = 0
        else:
            dx0 = (bot[0]-top[0])/(bot[1]-top[1])
            dz0 = (bot[2]-top[2])/(bot[1]-top[1])

        #if sortaequal(mid[1],top[1], 0.0001):
        if mid[1] == top[1]:
            dx1m = 0
            dz1m = 0
        else:
            dx1m = (mid[0]-top[0])/(mid[1]-top[1])
            dz1m = (mid[2]-top[2])/(mid[1]-top[1])

        #if sortaequal(bot[1],mid[1], 0.0001):
        if bot[1] == mid[1]:
            dx1b = 0
            dz1b = 0
        else:
            dx1b = (bot[0]-mid[0])/(bot[1]-mid[1])
            dz1b = (bot[2]-mid[2])/(bot[1]-mid[1])

        while yi < mid[1]:
            x1 += dx1m
            z1 += dz1m
            yi += 1
            x0 += dx0
            z0 += dz0
            draw_line(screen, x0,yi,z0, x1,yi,z1, colortmp)
        x1 = mid[0]
        yi = mid[1]
        z1 = mid[2]
        draw_line(screen, x0,yi,z0, x1,yi,z1, colortmp)
        while yi < bot[1]:
            x0 += dx0
            z0 += dz0
            x1 += dx1b
            z1 += dz1b
            yi += 1
            draw_line(screen, x0,yi,z0, x1,yi,z1, colortmp)
            
            

    def draw_polygon(p0,p1,p2, c):
        if env["shading_mode"]=="wireframe":
            draw_line(screen, p0[0],p0[1],p0[2], p1[0],p1[1],p1[1], c)
            draw_line(screen, p1[0],p1[1],p1[2], p2[0],p2[1],p2[2], c)
            draw_line(screen, p2[0],p2[1],p2[2], p0[0],p0[1],p0[2], c)
        else:
            scanlines(p0,p1,p2)

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
            draw_polygon(points[p], points[p+1], points[p+2], color)

        p+=3
        

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
    i = 1.0 + 1.0
    j = 1.0
    while i<=step:
        ti = i/step
        theta = 2 * math.pi * ti
        while j<=step:
            tj = j/step
            phi = 2 * math.pi * tj
            x = (r1*math.cos(theta) + r2) * math.cos(phi) + cx
            y = r1*math.sin(theta) + cy
            z = (r1*math.cos(theta) + r2) * math.sin(phi) + cz
            spts += [[x,y,z]]
            j+=1
        i+=1
        j=1.0
        
    #connect points
    i = -step
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

def draw_point(screen, color, x,y,z):
    x = int(x)
    y = int(y)
    if x >= len(zbuffer) or y>= len(zbuffer[0]) or x<0 or y<0: return
    if z > zbuffer[x][y]:
        zbuffer[x][y] = z
        plot(screen, color, x, y)

def draw_line( screen, x0, y0, z0, x1, y1, z1, color ):
    dx = x1 - x0
    dy = y1 - y0
    dz = z1 - z0
    
    if dx + dy < 0:
        dx = 0 - dx
        dy = 0 - dy
        dz = 0 - dz
        tmp = x0
        x0 = x1
        x1 = tmp
        tmp = y0
        y0 = y1
        y1 = tmp
        
    dzy = dz/dy if dy!=0 else 0
    dzx = dz/dx if dx!=0 else 0
    
    if dx == 0:
        y = y0
        z = z0
        while y <= y1:
            draw_point(screen,color,  x0,y,z)
            y += 1
            z += dzy
    elif dy == 0:
        x = x0
        z = z0
        while x <= x1:
            draw_point(screen,color, x,y0,z)
            x += 1
            z += dzx
    elif dy < 0:
        d = 0
        x = x0
        y = y0
        z = z0
        while x <= x1:
            draw_point(screen,color, x,y,z)
            if d > 0:
                y -= 1
                d -= dx
            x += 1
            z += dzx
            d -= dy
    elif dx < 0:
        d = 0
        x = x0
        y = y0
        z = z0
        while y <= y1:
            draw_point(screen,color, x,y,z)
            if d > 0:
                x -= 1
                d -= dy
            y += 1
            z += dzy
            d -= dx
    elif dx > dy:
        d = 0
        x = x0
        y = y0
        z = z0
        while x <= x1:
            draw_point(screen,color, x,y,z)
            if d > 0:
                y += 1
                d -= dx
            x += 1
            z += dzx
            d += dy
    else:
        d = 0
        x = x0
        y = y0
        z = z0
        while y <= y1:
            draw_point(screen,color, x,y,z)
            if d > 0:
                x += 1
                d -= dy
            y += 1
            z += dzy
            d += dx

