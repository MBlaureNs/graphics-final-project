from display import *
from matrix import *
from draw import *

def parse_file( fname, points, transform ):
    trans = [transform] #this is the only time i wish we
    pts = [points]      #used python 3 instead of python 2

    def do_l(args): #add line
        def inner(x0,y0,z0,x1,y1,z1):
            add_edge(pts[0],x0,y0,z0,x1,y1,z1)
        inner(int(args[0]),int(args[1]),int(args[2]),
              int(args[3]),int(args[4]),int(args[5]))

    def do_c(args): #add circle
        def inner(cx,cy,r):
            add_circle(pts[0],cx,cy,0,r,1000)
        inner(int(args[0]),int(args[1]),int(args[2]))

    def do_h(args): #add hermite
        def inner(x0,y0,x1,y1,x2,y2,x3,y3):
            add_curve(pts[0],x0,y0,x1,y1,x2,y2,x3,y3,1000,0)
        inner(int(args[0]),int(args[1]),int(args[2]),int(args[3]),
              int(args[4]),int(args[5]),int(args[6]),int(args[7]))

    def do_b(args): #add bezier
        def inner(x0,y0,x1,y1,x2,y2,x3,y3):
            add_curve(pts[0],x0,y0,x1,y1,x2,y2,x3,y3,1000,1)
        inner(int(args[0]),int(args[1]),int(args[2]),int(args[3]),
              int(args[4]),int(args[5]),int(args[6]),int(args[7]))

    def do_p(args): #add prism
        def inner(x,y,z,width,height,depth):
            add_prism(pts[0],x,y,z,width,height,depth)
        inner(int(args[0]),int(args[1]),int(args[2]),
              int(args[3]),int(args[4]),int(args[5]))

    def do_m(args): #add sphere
        def inner(x,y,radius):
            add_sphere(pts[0],x,y,0,radius,10)
        inner(int(args[0]),int(args[1]),int(args[2]))

    def do_d(args): #add torus
        def inner(x,y,radius1,radius2):
            #r1 - circle rad
            #r2 - rad from torus center to circle center
            add_torus(pts[0],x,y,0,radius1,radius2,12)
        inner(int(args[0]),int(args[1]),int(args[2]),int(args[3]))

    def do_i(args): #identity transform
        def inner():
            trans[0] = ident(new_matrix(4,4))
        inner()

    def do_s(args): #scale transform
        def inner(sx,sy,sz):
            trans[0] = matrix_mult(make_scale(sx,sy,sz),trans[0])
        inner(float(args[0]),float(args[1]),float(args[2]))

    def do_t(args): #translate transform
        def inner(tx,ty,tz):
            trans[0] = matrix_mult(make_translate(tx,ty,tz),trans[0])
        inner(float(args[0]),float(args[1]),float(args[2]))

    def do_x(args): #x rotation transform
        def inner(theta):
            trans[0] = matrix_mult(make_rotX(theta),trans[0])
        inner(float(args[0]))

    def do_y(args): #y rotation transform
        def inner(theta):
            trans[0] = matrix_mult(make_rotY(theta),trans[0])
        inner(float(args[0]))

    def do_z(args): #z rotation transform
        def inner(theta):
            trans[0] = matrix_mult(make_rotZ(theta),trans[0])
        inner(float(args[0]))

    def do_a(args): #apply transform
        def inner():
            pts[0] = matrix_mult(trans[0], pts[0])
        inner()

    def do_w(args): #clear points
        def inner():
            pts = [[]]
        inner()

    def do_v(args): #view image
        def inner():
            screen = new_screen()
            color = [0,255,0]
            draw_polygons(pts[0], screen, color)
            display(screen)
        inner()

    def do_g(args): #save image
        def inner(fname):
            screen = new_screen()
            color = [0,255,0]
            draw_polygons(pts[0], screen, color)
            save_extension(screen, fname)
        inner(args[0])

    #           {cmd:(argc,f), #description
    
    validcmds = {"l":(6,do_l), #add line to points
                 "c":(3,do_c), #add circle to points
                 "h":(8,do_h), #add hermite curve to points
                 "b":(8,do_b), #add bezier curve to points
                 "p":(6,do_p), #add rectangular prism to points
                 "m":(3,do_m), #add sphere munchkin to points
                 "d":(4,do_d), #add torus donut to points
                 
                 "i":(0,do_i), #set transform to identity matrix
                 "s":(3,do_s), #scale transform matrix
                 "t":(3,do_t), #translate transform matrix
                 "x":(1,do_x), #rotate transform about x-axis
                 "y":(1,do_y), #rotate transform about y-axis
                 "z":(1,do_z), #rotate transform about z-axis
                 "a":(0,do_a), #apply transform matrix to points
                 
                 "w":(0,do_w), #clear all edges from points
                 "v":(0,do_v), #draw point matrix lines to screen
                 "g":(1,do_g)} #save screen to file

    def isvalid(cmd):
        return cmd in validcmds.keys()
    def argcount(cmd):
        return validcmds[cmd][0]
    def parseerror():
        print "REKT"
    def execcmd(args):
        validcmds[args[0]][1](args[1:])
    def execcmds(cmds, points, transform):
        i = 0
        args = []
        while i<len(cmds):
            if len(args) == 0: 
                if not isvalid(cmds[i]): parseerror()
            args+=[cmds[i]]
            if len(args)>argcount(args[0]):
                execcmd(args)
                args=[]
            i+=1

    f = open(fname)
    ftxt = f.read()
    f.close()
    ctxt = ""
    for ln in ftxt.split("\n"):
        if "#" in ln: ln = ln[:ln.find("#")]
        ctxt += ln + "\n"
    cmds = ctxt.split()
    execcmds(cmds, points, transform)

