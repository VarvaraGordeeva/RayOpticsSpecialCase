# Description of primitive objects like point, vector, sphere, for the problem of refracted ray search

from cmath import pi
from re import X
from tkinter import CENTER
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

import numpy as np

class MyPoint():
    def __init__(self,X,Y,Z):
        self.X = X
        self.Y = Y
        self.Z = Z

class Vector():
    def __init__(self,start,X,Y,Z):
        if(isinstance(start,MyPoint)):
            self.start = start
        else:
            print(f"There should be MyPoint type")
        self.X = X
        self.Y = Y
        self.Z = Z
    
    def length(self):
        L = math.sqrt(self.X**2+self.Y**2+self.Z**2)
        return L

class Point2D():
    def __init__(self,X,Y):
        self.X = X
        self.Y = Y

class Vector2D():
    def __init__(self,start,X,Y):
        if(isinstance(start,Point2D)):
            self.start = start
        else:
            print(f"There should be Point2D type")
        self.X = X
        self.Y = Y
    
    def length(self):
        L = math.sqrt(self.X**2+self.Y**2)
        return L

class Sphere():
    def __init__(self,center,R):
        if (isinstance(center,MyPoint)):
            self.center = center
        else:
            print(f"There should be MyPoint type") 
        self.radius = R
    
    def intersect(self,obj):
        P = 0
        eps = 1e-8
        if (isinstance(obj,MyPoint)):
            delta = (obj.X-self.center.X)**2+(obj.Y-self.center.Y)**2+(obj.Z-self.center.Z)**2-self.radius**2
            if (delta<=eps):
                print(f"The point is on the sphere")
                P = obj
        if(isinstance(obj,Vector)):
            A = (obj.X-obj.start.X)**2+(obj.Y-obj.start.Y)**2+(obj.Z-obj.start.Z)**2
            B = 2*((obj.X-obj.start.X)*(obj.start.X-self.center.X)
                  +(obj.Y-obj.start.Y)*(obj.start.Y-self.center.Y)
                  +(obj.Z-obj.start.Z)*(obj.start.Z-self.center.Z))
            C = ((obj.start.X-self.center.X)**2+(obj.start.Y-self.center.Y)**2
                 +(obj.start.Z-self.center.Z)**2-self.radius**2)
            D = B**2-4*A*C
            if D>0:
                t1 = (-B+D**0.5)/2/A
                t2 = (-B-D**0.5)/2/A
                p1 = MyPoint((obj.X-obj.start.X)*t1+obj.start.X,
                    (obj.Y-obj.start.Y)*t1+obj.start.Y,
                    (obj.Z-obj.start.Z)*t1+obj.start.Z)
                p2 = MyPoint((obj.X-obj.start.X)*t2+obj.start.X,
                    (obj.Y-obj.start.Y)*t2+obj.start.Y,
                    (obj.Z-obj.start.Z)*t2+obj.start.Z)
                l1 = math.sqrt((p1.X-obj.start.X)**2+(p1.Y-obj.start.Y)**2+(p1.Z-obj.start.Z)**2)
                l2 = math.sqrt((p2.X-obj.start.X)**2+(p2.Y-obj.start.Y)**2+(p2.Z-obj.start.Z)**2)
                if l1<l2:
                    P = p1
                else:
                    P = p2
        return P

    def intersect1(self,obj): # returns a section vector for an initial vector and a sphere
        V = 0
        P = 0
        
        if(isinstance(obj,Vector)):
            A = (obj.X-obj.start.X)**2+(obj.Y-obj.start.Y)**2+(obj.Z-obj.start.Z)**2
            B = 2*((obj.X-obj.start.X)*(obj.start.X-self.center.X)
                  +(obj.Y-obj.start.Y)*(obj.start.Y-self.center.Y)
                  +(obj.Z-obj.start.Z)*(obj.start.Z-self.center.Z))
            C = ((obj.start.X-self.center.X)**2+(obj.start.Y-self.center.Y)**2
                 +(obj.start.Z-self.center.Z)**2-self.radius**2)
            D = B**2-4*A*C
            if D>0:
                t1 = (-B+D**0.5)/2/A
                t2 = (-B-D**0.5)/2/A
                p1 = MyPoint((obj.X-obj.start.X)*t1+obj.start.X,
                    (obj.Y-obj.start.Y)*t1+obj.start.Y,
                    (obj.Z-obj.start.Z)*t1+obj.start.Z)
                p2 = MyPoint((obj.X-obj.start.X)*t2+obj.start.X,
                    (obj.Y-obj.start.Y)*t2+obj.start.Y,
                    (obj.Z-obj.start.Z)*t2+obj.start.Z)
                l1 = math.sqrt((p1.X-obj.start.X)**2+(p1.Y-obj.start.Y)**2+(p1.Z-obj.start.Z)**2)
                l2 = math.sqrt((p2.X-obj.start.X)**2+(p2.Y-obj.start.Y)**2+(p2.Z-obj.start.Z)**2)
                if l1<l2:
                    P = p1
                    E = p2
                else:
                    P = p2
                    E = p1
                V = Vector(P,E.X,E.Y,E.Z)    
        else:
            print(f"There should be Vector type")   
        return V

#

def present(obiekt,c):
    
    fig = plt.figure(figsize = (8,8))
    ax = plt.axes(projection='3d')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    ax.set_zlim([-2,2])
    
    if (isinstance(obiekt,MyPoint)):
        ax.scatter(obiekt.X,obiekt.Y,obiekt.Z,color=c,s=300)
    
    if (isinstance(obiekt,Sphere)):
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
        x = obiekt.radius*np.cos(u)*np.sin(v)+obiekt.center.X
        y = obiekt.radius*np.sin(u)*np.sin(v)+obiekt.center.Y
        z = obiekt.radius*np.cos(v)+obiekt.center.Z
        #ax.plot_wireframe(x, y, z, color=c)
        ax.plot_surface(x, y, z, color=c)

    if(isinstance(obiekt,Vector)):
        #x = [obiekt.start.X,obiekt.X+obiekt.start.X]
        #y = [obiekt.start.Y,obiekt.Y+obiekt.start.Y]
        #z = [obiekt.start.Z,obiekt.Z+obiekt.start.Z]
        x = [obiekt.start.X,obiekt.X]
        y = [obiekt.start.Y,obiekt.Y]
        z = [obiekt.start.Z,obiekt.Z]
        ax.plot(x,y,z,color=c)

def mulcr(obj1,obj2): # cross multiplication of vectors

    if (isinstance(obj1,Vector))and(isinstance(obj2,Vector)):
        a = [obj1.X,obj1.Y,obj1.Z]
        b = [obj2.X,obj2.Y,obj2.Z]
        c = np.cross(a,b)
        P = Vector(obj1.start,c[0],c[1],c[2])
    else:
            print(f"There should be Vector type")
    return P

# P = MyPoint(0,0,0)
# Sp = Sphere(P,1)
# present(P,"r")
# plt.show()