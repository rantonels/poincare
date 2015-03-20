#!/usr/bin/python

import numpy as np


P = 3
Q = 7

iterations = 9

def apply_mobius(M,z):
    return (M[0,0]*z + M[0,1])/(M[1,0]*z + M[1,1])

class Point:
	def __init__(self,ID,M):
		self.ID = ID
		self.M = M
		self.done = False

        def position(self):
            return apply_mobius(self.M,0.+0.j)

def distance(u,v):
    up = u.position()
    vp = v.position()

    delta = 2 * abs(up - vp)**2 / (1-abs(up)**2)*(1-abs(vp)**2)

    return np.arccosh(1+delta)


points = {}

points[0] = Point(0,np.eye(2))

links = set([])


#c = np.arccosh(1/ ( np.tan(np.pi/P) * np.tan(np.pi/Q)))
#
#
#a = np.arcsinh(np.sinh(c)*np.sin(np.pi/P))
#
#
#a2 = 2*a
#
#
#y = np.cosh(a2)
#
#aa = np.sqrt(y-1)/np.sqrt(y+1)

s = np.sin(np.pi/P)
c = np.cos(np.pi/Q)

r = 1/np.sqrt(((c*c)/(s*s))-1)

d = 1/np.sqrt(1-((s*s)/(c*c)))

l = r*r/d

aa = d-l

aa = (aa + 0j)

T = np.matrix( [[ 1 ,  aa],  [aa.conjugate(),1]])

B = np.matrix( [[ -1, 0 ],[0,1]] )

for it in range(iterations):
    print it
    undone = [pid for pid in points if not points[pid].done]
    for pid in undone:
        p = points[pid]
        for i in range(P):
            R = np.matrix( [[ np.exp(i*2*np.pi/P * 1j)  , 0 ], [0, 1]])
            
            # one applies, in order: 180 rotation, translation of unit, rotation about parent point, transformation of parent

            F = p.M * R * T * B

            nupoint = Point(len(points),F)

            #search for point in db
            found = False
            for i in points:
                if distance(nupoint,points[i]) < 0.01:
                    found = True
                    companion_index = i
                    break

            if not found:
                points[nupoint.ID] = nupoint
                links.add( (p,nupoint) ) 
            else:
                links.add( (p,points[companion_index]) ) 



        points[pid].done = True

f = open('points','w')
for pid in points:
    z = points[pid].position()
    f.write("%f\t%f\t%f\n"%(pid,z.real,z.imag))
f.close()

f = open('links','w')
for l in links:
    p0 = l[0].position()
    p1 = l[1].position()
    diff = p1-p0
    f.write("%f\t%f\t%f\t%f\n"% (   p0.real, p0.imag, p1.real, p1.imag) )
f.close()
