from geo import *
from math import *

def checkf(value,correct_value,threshold = 1e-14):
    if abs(value - correct_value) > threshold:
        raise RuntimeError("Expected %g, got %g" % (correct_value,value))
    else:
        print "ok",

print "Testing:",

# Test construction and measurements

origin = Point(0,0,0)
xaxis = Line(origin,Point(1,0,0))
yaxis = Line(origin,Point(0,1,0))
zaxis = Line(origin,Point(0,0,1))
xy_plane = Plane(origin,Point(1,0,0),Point(0,1,0))

p1 = Point(1,1,1)
p2 = Point(1,0,0)
checkf(p1.distance_to(p2),sqrt(2))

L = Line(p1,p2)
checkf(L.distance_to(p1),0)

checkf(p1.distance_to(xy_plane),1)

checkf(L.angle_to(xy_plane),pi/4)

# Make sure that the fitted line coincides with L for three points on L.
Lfit = Line(p1,p2,Point(1,0.5,0.5))
checkf(L.angle_to(Lfit),0)
checkf(L.distance_to(Lfit),0)

# Midpoints
p = Point(-1,-2,-3).midpoint_to(Point(1,2,3))
checkf(p.distance_to(origin),0)

p = xaxis.midpoint_to(yaxis)
checkf(p.distance_to(origin),0)

p = zaxis.midpoint_to(xy_plane)
checkf(p.distance_to(origin),0)




#test movements

# pure translation
mt = Movement(origin,p1)
checkf(p1.moved(mt).distance_to(origin),2*p1.distance_to(origin))

#pure rotation
mr = Movement(xaxis,yaxis) # rotate xaxis onto yaxis
checkf(xaxis.moved(mr).distance_to(yaxis),0)
checkf(origin.moved(mr).distance_to(origin),0)


print "done."


