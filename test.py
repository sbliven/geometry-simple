from geo import *
from math import *

def checkf(value,correct_value,threshold = 1e-14):
    if abs(value - correct_value) > threshold:
        raise RuntimeError("Expected %g, got %g" % (correct_value,value))

origin = Point(0,0,0)
p1 = Point(1,1,1)
p2 = Point(1,0,0)
checkf(p1.distance_to(p2),sqrt(2))

L = Line(p1,p2)
checkf(L.distance_to(p1),0)

xy_plane = Plane(origin,Point(1,0,0),Point(0,1,0))
checkf(p1.distance_to(xy_plane),1)

checkf(L.angle_to(xy_plane),pi/4)


