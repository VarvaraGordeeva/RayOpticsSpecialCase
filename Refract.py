# Module for calculation of coordinates and incidence angle of a ray refracted on a bubble

import matplotlib.pyplot as plt
import math
import numpy as np
import sys
import Primitives as prim

# consts:
n1 = 1.333 # refraction index of water
n2 = 1.000293 # refraction index of air
# Zd = 6 # sensor coordinate

def intPoint(ve:prim.Vector,sp:prim.Sphere,b:bool):
    # intersection point of a vector (ve) and a sphere (sp)
    # b is boolean for first or second point (before or after refraction)
    inP = prim.MyPoint(0,0,-1)
    if isinstance(ve,prim.Vector)and isinstance(sp,prim.Sphere):
        xc = sp.center.X
        yc = sp.center.Y
        zc = sp.center.Z
        R = sp.radius
        A0 = ve.X**2+ve.Y**2+ve.Z**2
        B0 = 2*(ve.X*(ve.start.X-xc)+ve.Y*(ve.start.Y-yc)+ve.Z*(ve.start.Z-zc))
        C0 = (ve.start.X-xc)**2+(ve.start.Y-yc)**2+(ve.start.Z-zc)**2-R*R
        D0 = B0*B0-4*A0*C0
        if D0>0:
            t1 = (-B0-math.sqrt(D0))/2/A0
            t2 = (-B0+math.sqrt(D0))/2/A0
            if not b:
                x1 = ve.X*t1+ve.start.X
                y1 = ve.Y*t1+ve.start.Y
                z1 = ve.Z*t1+ve.start.Z
                l1 = math.sqrt((x1-ve.start.X)**2+(y1-ve.start.Y)**2+(z1-ve.start.Z)**2)
                x2 = ve.X*t2+ve.start.X
                y2 = ve.Y*t2+ve.start.Y
                z2 = ve.Z*t2+ve.start.Z
                l2 = math.sqrt((x2-ve.start.X)**2+(y2-ve.start.Y)**2+(z2-ve.start.Z)**2)
                if l1<l2:
                    inP.X = x1
                    inP.Y = y1
                    inP.Z = z1
                else:
                    inP.X = x2
                    inP.Y = y2
                    inP.Z = z2
            else:
                t = max(t1,t2)
                inP.X = ve.X*t+ve.start.X
                inP.Y = ve.Y*t+ve.start.Y
                inP.Z = ve.Z*t+ve.start.Z
        else:
            return inP
    else:
        print('The objects should be Vector and Sphere')
    return inP

def incCos(ve1:prim.Vector,ve2:prim.Vector):
    # cosine of the incidence ray (between two vectors)
    dotmul = ve1.X*ve2.X+ve1.Y*ve2.Y+ve1.Z*ve2.Z
    l1 = ve1.length()
    l2 = ve2.length()
    cos = dotmul/(l1*l2)
    return cos

def refCos(al1,n01,n02):
    # cosine of refraction angle
    # al1 is cosine of incidence angle
    cos = np.sqrt(1-(n01/n02)**2*(1-al1**2))
    return cos

def refVec(p:prim.Vector,ms:prim.MyPoint,pl:prim.MyPoint,cos,ob:prim.Sphere):
    # search of refracted ray
    # p is the vectors of perpendicular starting in intersection point,
    # ms is the incidence vector m start
    # pl is incidence-refraction plane (x,y,z) for (A,B,C)
    # cos is cosine of the refraction angle
    # ob is the sphere (bubble)
    P1y = p.start.Y
    P1z = p.start.Z
    K1 = (pl.Z*p.X-pl.X*p.Z)/(pl.X*p.Y-pl.Y*p.X)
    K2 = ob.radius*pl.X*cos/(pl.X*p.Y-pl.Y*p.X) 
    K3 = (pl.Y*K1-pl.Z)/pl.X
    K4 = pl.Y*K2/pl.X
    aa = K3*K3+K1*K1+1
    bb = 2*(K1*K2+K3*K4)
    cc = K2*K2+K4*K4-1
    discr = bb**2-4*aa*cc
    Ap = 1/(ob.center.Y-P1y)
    Bp = -1/(ob.center.Z-P1z)
    Cp = P1z/(ob.center.Z-P1z)-P1y/(ob.center.Y-P1y)
    inc = np.sign(Ap*ms.Y+Cp)
    m = prim.Vector(p.start,0,0,0)
    if discr>0:
        m.Z=(-bb+np.sqrt(discr))/(2*aa)
        m.X=-m.Z*K3-K4
        m.Y=m.Z*K1+K2
        ref = np.sign(Ap*(m.Y+P1y)+Bp*(m.Z+P1z)+Cp)
        # refracted ray choise
        if inc == ref:
            m.Z=(-bb-np.sqrt(discr))/(2*aa)
            m.X=-m.Z*K3-K4
            m.Y=m.Z*K1+K2
            # print('Refracted ray:')
            # print(m)
    else:
        # print('Unable to build refracted ray')
        return(m)
    return(m)

def GetAngle(ve:prim.Vector,ob:prim.Sphere,Zd):
    # ve is initial (incident) ray vector
    # ob is an illuminated object (sphere)
    # Zd is a sensor Z-coordinate
    
    # initialization
    xc = ob.center.X
    yc = ob.center.Y
    zc = ob.center.Z
    x0 = ve.start.X
    y0 = ve.start.Y
    # result
    res = np.zeros([3])
    # error
    eps = ob.radius*1e-3
    
    # first intersection point
    P1 = prim.MyPoint(0,0,-1)
    P1 = intPoint(ve,ob,0)
    if P1.Z == -1:
        res[0] = x0
        res[1] = y0
        res[2] = 1
        # print('There is no intersection point')
        return(res)
    # print('Coordinates of the first intersection point: x1 = ',P1.X,'y1 =',P1.Y,'z1 = ',P1.Z)

    # incidence angle
    p1 = prim.Vector(ob.center,P1.X-xc,P1.Y-yc,P1.Z-zc)
    m01 = prim.Vector(P1,-ve.X,-ve.Y,-ve.Z)
    cos_0 = incCos(p1,m01)
    # print('Length of m0 =',m01.length())
    # print('incidence cos(alpha0)=',cos_0)
    # check for pure reflection
    if cos_0>np.sqrt(1-(n2/n1)**2):
        # print('incidence cos(alpha0) = ',cos_0)
        cos_1 = refCos(cos_0,n1,n2)
        # print('Refracted cos(alpha1)=',cos_1)
    else:
        # print('Pure reflection!')
        res[0] = x0
        res[1] = y0
        res[2] = 0
        return(res)
    
    # incidence-refraction plane
    sig = prim.MyPoint(0,0,-1) # actually (x,y,z) for (A,B,C) of the plane
    D = np.ones([3])
    Xx = np.array([[x0,y0,0],[P1.X,P1.Y,P1.Z],[xc,yc,zc]])
    delta = np.linalg.det(Xx)
    if delta !=0:
        coef = np.linalg.solve(Xx,D)
        sig.X=coef[0]
        sig.Y=coef[1]
        sig.Z=coef[2]
    else:
        # print('delta = ',delta,'. The matrix is degenerate')
        res[0] = ve.X*Zd/ve.Z+ve.start.X
        res[1] = ve.Y*Zd/ve.Z+ve.start.Y
        res[2] = ve.Z
        return(res)
    # print('Incidence-refraction plane: (',sig.X,sig.Y,sig.Z,')')
    # refracted ray vector
    m1 = prim.Vector(P1,0,0,0)
    p11 = prim.Vector(P1,-p1.X,-p1.Y,-p1.Z)
    m1 = refVec(p11,ve.start,sig,cos_1,ob)
    if m1.length() ==0:
        # print('Unable to build 1-refracted ray m1')
        return(res)
    # print('Refracted ray: (',m1.X,m1.Y,m1.Z,')')

    # second incidence angle
    P2 = intPoint(m1,ob,1)
    # print('Coordinates of the second intersection point: ',P2)
    p2 = prim.Vector(P2,P2.X-xc,P2.Y-yc,P2.Z-zc)
    # print('m1 = ',m1.length())
    # print('p2 = ',p2.length())
    cos_2 = incCos(p2,m1)
    # print('Second incidence cos(alpha2)=',cos_2)
    # second refraction angle
    cos_3 = refCos(cos_2,n2,n1)
    # print('Refracted cos(alpha3)=',cos_3)
    
    # second refraction vector m2
    m2 = refVec(p2,P2,sig,cos_3,ob)
    # if m2.length == 0:
    if m2.Z == 0:
        # print('Unable to build 2-refracted ray m2')
        return(res)
    # print('Refracted ray: (',m2.X,m2.Y,m2.Z,')')
    
    # incidence point coordinate search
    t = (Zd-P2.Z)/m2.Z
    P3x = m2.X*t+P2.X
    P3y = m2.Y*t+P2.Y
    # print('Incident angle cosine = ',m2.Z)
    # print('Coordinates of the refracted ray incidence:')
    # print('P3x = ',P3x,'; P3y = ',P3y)
    
    res[0]=P3x
    res[1]=P3y
    res[2]=m2.Z
    return(res)

# #initial values:

# # point source coordinates (z0=0)
# P0 = prim.MyPoint(1.2,1.5,0)
# m0 = prim.Vector(P0,0,0,1)
# # sensor coordinate
# Zd = 6
# # sphere coordinates
# xc = 1.5
# yc = 2
# zc = 1.5
# R = 1.5
# Pc = prim.MyPoint(xc,yc,zc)
# Sp = prim.Sphere(Pc,R)

# dF = GetAngle(m0,Sp,Zd)
# print('Coordinates and angle cosine =',dF)