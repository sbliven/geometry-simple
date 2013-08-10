"""Geometry-simple

http://code.google.com/p/geometry-simple/

Created by uekstrom
"""
#The MIT License (MIT)
#
#Copyright (c) <year> <copyright holders>
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

from numpy import *
import math, numpy.linalg, copy
import fractions

# Expose only some parts to "from geo import *"
__all__ = ["use_degrees","use_radians","Point","Line","Plane","Movement"]

angular_unit = 1.0 #units/radian

def use_degrees():
    """ Use degrees for all geo functions """
    global angular_unit
    angular_unit = 180.0/math.pi

def use_radians():
    """ Use radians for all geo functions (default)"""
    global angular_unit
    angular_unit = 1.0

def dot(x,y):
    """ Dot product of arrays x and y"""
    return inner(x,y)

def abs2(x):
    """ the sum of squares of array x """
    return sum(x**2)

def normalized(x):
    """ normalize array x to length 1"""
    return x/math.sqrt(abs2(x))

def orthogonalized_to(x,d):
    """Return a copy of array x orthogonalized to array d != 0"""
    d = normalized(d)
    return x - dot(x,d)*d

def dual(v):
    """Return two unit vectors orthogonal to array v"""
    if abs2(v) > 1e-20:
        if v[0] < 0.7:
            n1 = normalized(orthogonalized_to(array([1,0,0],'d'),v))
        else:
            n1 = normalized(orthogonalized_to(array([0,1,0],'d'),v))
        n2 = cross(v,n1)
        return [n1,n2]
    else:
        return [array([1,0,0],'d'),array([0,1,0],'d')]

def qmul(q1,q2):
    """Take q1 and q2 to be quaternions, and multiply them accordingly"""
    v1 = q1[1:]
    v2 = q2[1:]
    x = q1[0]*v2 + q2[0]*v1 + cross(v1,v2)
    return array([q1[0]*q2[0]-dot(v1,v2),x[0],x[1],x[2]])

def qconj(q):
    qc = -q
    qc[0] *= -1
    return qc

def qrotate(q,v):
    qv = array([0,v[0],v[1],v[2]])
    return qmul(q,qmul(qv,qconj(q)))[1:]

def qrotor(axis,angle):
    axis = math.sin(angle/2)*normalized(axis)
    return array([math.cos(angle/2), axis[0], axis[1], axis[2]])

class Point(object):
    def __init__(self, *x_or_xyz):
        """Create a point from a list of 3 coordinates or 3 individual
        coordinates."""
        if len(x_or_xyz) == 3:
            self.r = array(x_or_xyz,'d')
        elif len(x_or_xyz) == 1 and len(x_or_xyz[0]) == 3:
            self.r = array(x_or_xyz[0],'d')
        else:
            raise TypeError("Invalid arguments to Point()")
    def moved(self, m):
        return m.on_point(self)
    def distance_to(self, obj):
        if isinstance(obj,Point):
            return math.sqrt(abs2(self.r - obj.r))
        else:
            raise TypeError("Invalid type in Point.distance_to()")
    def __str__(self):
        return "(%g, %g, %g)" % tuple(map(float,self.r))
    def __repr__(self):
        return "Point(%g, %g, %g)" % tuple(map(float,self.r))
    def coordinates(self):
        return map(float,self.r)
    def projected_on(self, obj):
        if isinstance(obj,Line):
            return Point(obj.r + dot(obj.t,self.r - obj.r)*obj.t)
        elif isinstance(obj,Plane):
            dr = self.r - obj.r
            return Point(obj.r + orthogonalized_to(dr,obj.n))
        else:
            raise TypeError("Cannot project point onto object")
    def distance_to(self, obj):
        if isinstance(obj,Point):
            return math.sqrt(abs2(self.r - obj.r))
        else:
            p = self.projected_on(obj)
            return self.distance_to(p)
    def midpoint_to(self,obj):
        """Return a point in the middle of the shortest line connecting this and obj."""
        if isinstance(obj,Point):
            return Point(0.5*(self.r + obj.r))
        else:
            return obj.midpoint_to(self)

def pointset_mass_distribution(points):
    cm = zeros((1,3))
    for p in points:
        cm += p.r
    cm /= len(points)
    A = asmatrix(zeros((3,3)))
    for p in points:
        r = asmatrix(p.r - cm)
        A += r.transpose()*r
    return asarray(cm).reshape(3),A

class Line(object):
    def __init__(self, *points):
        """Create an infinite line from at least two points.
        Accepts either two points or a list of points. If more than
        two points are given the line will be a least square fit of
        the point set, in which case the (sign of the) direction is
        undefined."""
        if len(points) == 1:
            points = points[0]
        if len(points) == 2:
            self.r = array(points[0].r)
            self.r2 = array(points[1].r)
            self.t = normalized(self.r2 - self.r)
        elif len(points) > 2:
            r_cm, A = pointset_mass_distribution(points)
            # assume eigh() returns sorted vectors, take the one with
            # largest eigenvalue
            val, vec = numpy.linalg.eigh(A)
            self.t = asarray(vec[:,2]).reshape(3)
            l = math.sqrt(val[2]/2)
            self.r = r_cm
            self.r2 = r_cm + l*self.t
        else:
            raise ValueError("Too few arguments to Line()")

    def points(self):
        """Return two points defining the line"""
        return [Point(self.r),Point(self.r2)]

    def moved(self, m):
        p = self.points()
        return Line(m.on_point(p[0]),m.on_point(p[1]))

    def projected_on(self, plane):
        p = self.points()
        return Line(p[0].projected_on(plane),
                    p[1].projected_on(plane))

    def distance_to(self, obj):
        if isinstance(obj,Point):
            return obj.distance_to(self)
        elif isinstance(obj,Line):
            d = obj.r - self.r
            n = cross(self.t,obj.t)
            if abs2(n) < 1e-16: # parallel lines
                return math.sqrt(abs2(cross(d,self.t)))
            else:
                return abs(dot(d,n))/math.sqrt(abs2(n))
        else:
            # Line-plane distance is only non-zero for exactly parallel objects
            # Because of numerical errors this is unlikely to happen, so always fail.
            raise ValueError("Will not calculate line-plane distance")

    def angle_to(self, obj):
        if isinstance(obj,Line):
            return angular_unit*math.acos(min(1,abs(dot(self.t,obj.t))))
        elif isinstance(obj,Plane):
            return angular_unit*(math.pi/2 - math.acos(min(1,abs(dot(self.t,obj.n)))))
        else:
            raise ValueError("Cannot calculate angle to object of this type")

    def midpoint_to(self,obj):
        """Return a point in the middle of the shortest line connecting this and obj."""
        if isinstance(obj,Point):
            return obj.midpoint_to(obj.projected_on(self))
        elif isinstance(obj,Line):
            d = obj.r - self.r
            t1t2 = dot(self.t,obj.t)
            if abs(abs(t1t2)-1) < 1e-12: #parallel case                
                d = orthogonalized_to(d,self.t)
                return Point(self.r + 0.5*d)
            else:
                t1d = dot(d,self.t)
                t2d = dot(d,obj.t)
                s = (t1t2*t2d - t1d)/(t1t2**2-1)
                u = (t1t2*t1d - t2d)/(t1t2**2-1)
                return Point(0.5*(obj.r + u*obj.t  + self.r + s*self.t))                
        else:
            return obj.midpoint_to(self)

    def __repr__(self):
        p = self.points()
        return "Line(%s, %s)" % (repr(p[0]),repr(p[1]))

    def __str__(self):
        return repr(self)

    def dual(self):
        """Return a plane such that plane.normal() == self"""
        d = dual(self.t)
        return Plane(Point(self.r),Point(self.r + d[0]), Point(self.r + d[1]))

class Plane(object):
    """ Represents a 3D plane

    Stored internally as a point (self.r) and a normal vector (self.n)
    """

    def __init__(self, *points):
        """Create a plane.
        
        Accepts multiple types of arguments, which can either be passed as
        multiple arguments or as a list:
         * A Line and a Point. Creates a Plane purpendicular to the line and
           containing the point.
         * 3 Points. Creates a Plane containing all three points
         * >3 Points. Fit a plane to the points using linear regression. The
           sign of the normal is undefined.
        """
        # accept arguments in a list
        if len(points) == 1:
            points = points[0]

        if len(points) == 2:
            # point and line
            if isinstance(points[0],Line) and isinstance(points[1],Point):
                p = points[1]
                l = points[0]
            elif isinstance(points[1],Line) and isinstance(points[0],Point):
                p = points[0]
                l = points[1]
            else:
                raise ValueError("Invalid arguments to Plane(%s)"%",".join([p.__class__.__name__ for p in points]))
            self.r = p.r
            self.n = l.t
        elif len(points) == 3:
            self.r = points[0].r
            n = cross(points[1].r - points[0].r,
                            points[2].r - points[0].r)
            if abs2(n) < 1e-28:
                raise ValueError("Degenerate points in Plane()")
            self.n = normalized(n)
        elif len(points) > 3:
            self.r, A = pointset_mass_distribution(points)
            # again assume that the vectors are sorted by value
            val, vec = numpy.linalg.eigh(A)
            self.n = asarray(vec[:,0]).reshape(3)
        else:
            raise ValueError("Too few arguments to Plane()")

    def points(self):
        """ Return three points on the plane.

        The first point defines an origin of the plane, while the second two
        points form an orthogonal basis relative to the first point.
        """
        d = dual(self.n)
        return [Point(self.r),Point(self.r + d[0]),Point(self.r + d[1])]

    def moved(self, m):
        p = self.points()
        return Plane(m.on_point(p[0]),m.on_point(p[1]),m.on_point(p[2]))

    def distance_to(self, obj):
        """ Calculates the distance to a point """
        if isinstance(obj,Point):
            return obj.distance_to(self)
        else:
            raise ValueError("Will not calculate line-plane or plane-plane distance")
    def angle_to(self, obj):
        """ Calculates the angle formed between this plane and another Plane or Line (0 to pi/4)"""
        if isinstance(obj,Line):
            return obj.angle_to(self)
        elif isinstance(obj,Plane):
            return angular_unit*math.acos(min(1,abs(dot(self.n,obj.n))))
        else:
            raise ValueError("Cannot calculate angle to object of this type")

    def midpoint_to(self,obj):
        """Return a point in the middle of the shortest line connecting this and obj."""
        if isinstance(obj,Point):
            return obj.midpoint_to(obj.projected_on(self))
        elif isinstance(obj,Line):
            if abs(abs(dot(self.n,obj.t)) - 1) < 1e-12: #line normal to the plane
                d = orthogonalized_to(obj.r - self.r,self.n)
                return Point(self.r+d)
            else:
                return obj.midpoint_to(obj.projected_on(self))
        else:
            raise NotImplemented("Plane-Plane midpoint not implemented")

    def normal(self):
        """Return a line normal to the plane"""
        return Line(Point(self.r),Point(self.r+self.n))

    def __repr__(self):
        """ Three-point representation of a Plane """
        p = self.points()
        return "Plane(%s, %s, %s)" % (repr(p[0]),repr(p[1]),repr(p[2]))

    def __str__(self):
        """ Equation representation of the Plane. Tries to simplify
        coefficients where possible. """
        #degenerate cases: normal along an axis
        absn = abs(self.n)
        if absn[0] < 1e-14:
            if absn[1] < 1e-14:
                # normal to z
                return "z = {}".format(self.r[2])
            elif absn[2] < 1e-14:
                # normal to y
                return "y = {}".format(self.r[1])
        elif absn[1] < 1e-14 and self.n[2] < 1e-14:
            # normal to x
            return "x = {}".format(self.r[0])

        coefs = self.coef(100)
        str = ""
        for coef,var in zip(coefs,"xyz"):
            if abs(coef)>0:
                if str:
                    if abs(coef) != 1:
                        str += " {} {}*{}".format("+" if coef>=0 else "-", abs(coef), var)
                    else:
                        str += " {} {}".format("+" if coef>=0 else "-", var)
                else:

                    if abs(coef) != 1:
                        str += "{}*{}".format(coef,var)
                    else:
                        str += var
        str += " = {}".format(coefs[3])
        return str

    def separates(self,p1,p2):
        """Test if the plane separates (is inbetween) two points"""
        return dot(p1.r - self.r,self.n)*dot(p2.r - self.r,self.n) < 0

    def coef(self,maxcoef=None):
        """coef(array,int) -> array
        Returns a tuple with coefficients (a,b,c,d) which satisfy the equation

            ax+by+cz=d.

        Since the coordinates can be scaled arbitrarily, by default they are
        normalized such that sqrt(a**2 + b**2 + c**2) = 1, where the sign is
        chosen based on the orientation of the plane.

        As an alternative, when maxcoef is specified the method attempts to
        find integer coefficients which approximate the plane. If the real
        coefficients of the plane are irrational or greater than maxcoef, the
        method will give up and return real-valued coefficients.
        """
        tol = 1e-14 #machine precision

        n = append(self.n,dot(self.n,self.r)) #real-valued abcd
        if maxcoef > 0:
            nr = n*max(abs(n)) # should be rational: [a,b,c]/|v|^2*max(a,b,c)
            nf = [fractions.Fraction(c).limit_denominator(maxcoef**2) for c in nr]
            if all([abs(nr[i]-nf[i])<tol for i in xrange(len(nr))]):
                #Good rational approximation for all coef
                factor = reduce(fractions.gcd, nf )
                vi = [f/factor for f in nf]
                assert( all([f.denominator == 1 for f in vi]) )
                n = array([f.numerator for f in vi])
            #No rational approximation, so use real

        # make first coef positive by convention
        if n[0] < 0:
            n = -n

        return n
    def orientation(self,point):
        """ orientation(Point) -> int

        Returns 1,0, or -1 based on whether the input point is above, on, or
        below the given plane. The orientation of the plane depends on its
        conditions at creation, but should remain consistent following
        movements.
        """
        p = point.r - self.r #vector from our origin to the point
        d = dot(self.n,p)
        return int(sign(d))


class Movement(object):
    """ Represents an affine transform between two objects

    Stored internally as a translation (self.dr) and a rotation quaternion (self.q)
    """
    def __init__(self, from_obj, to_obj):
        """Create an affine transformation (rotation+translation) that
        takes from_obj to to_obj in some sense.

        From    To      Movement
        =====   =====   ========
        Point   Point   Translate point onto point
        Point   Line    Translate point onto closest point on line
        Point   Plane   Translate point onto closest point on plane
        Line    Line    Rotate line1 to be parallel to line2, then translate to superimpose
        Line    Plane   Align line with its projection onto plane
        Plane   Plane   Rotate plane1 to be parallel to plane2, then translate to superimpose
        """
        # self.dr is always defined
        self.q = None # rotation quaternion, if applicable
        if isinstance(from_obj,Point) and isinstance(to_obj,Point):
            self.dr = to_obj.r - from_obj.r
        elif isinstance(from_obj,Point) and isinstance(to_obj,Line):
            # move point to closest point on line
            p = from_obj.projected_on(to_obj)
            self.dr = p.r - from_obj.r
        elif isinstance(from_obj,Point) and isinstance(to_obj,Plane):
            # move point to nearest point in plane
            p = from_obj.projected_on(to_obj)
            self.dr = p.r - from_obj.r
        elif isinstance(from_obj,Line) and isinstance(to_obj, Point):
            Movement.__init__(self,to_obj,from_obj)
            self.dr = -self.dr
        elif isinstance(from_obj,Line) and isinstance(to_obj,Line):
            # superimpose lines
            if dot(from_obj.t,to_obj.t) < 1 - 1e-14:
                self.q = qrotor(cross(from_obj.t,to_obj.t),
                                math.acos(dot(from_obj.t,to_obj.t)))
                self.dr = orthogonalized_to(to_obj.r - qrotate(self.q,from_obj.r),to_obj.t)
            else:
                self.dr = orthogonalized_to(to_obj.r - from_obj.r,to_obj.t)
        elif isinstance(from_obj,Line) and isinstance(to_obj,Plane):
            # move line onto its projection in the plane
            lp = from_obj.projected_on(to_obj)
            return Movement.__init__(self,from_obj,lp)
        elif isinstance(from_obj,Plane) and isinstance(to_obj,Point):
            # inverse of Point to plane motion
            Movement.__init__(self,to_obj,from_obj)
            self.q = qconf(self.q)
            self.dr = -qrotate(self.q,self.dr)
        elif isinstance(from_obj,Plane) and isinstance(to_obj,Line):
            # inverse of Point to line motion
            Movement.__init__(self,to_obj,from_obj)
            self.q = qconf(self.q)
            self.dr = -qrotate(self.q,self.dr)
        elif isinstance(from_obj,Plane) and isinstance(to_obj,Plane):
            # superimpose planes
            if dot(from_obj.n,to_obj.n) < 1 - 1e-14:
                self.q = qrotor(cross(from_obj.n,to_obj.n),
                                math.acos(dot(from_obj.n,to_obj.n)))
                self.dr = to_obj.r - qrotate(self.q,from_obj.r)
            else:
                self.dr = orthogonalized_to(to_obj.r - from_obj.r,to_obj.n)
        else:
            raise TypeError("Invalid arguments to Movement(%s,%s)"%(from_obj.__class__.__name__,to_obj.__class__.__name__))
    def on_point(self, p):
        """Private function to move points"""
        r = p.r
        if self.q != None:
            r = qrotate(self.q,r)
        if self.dr != None:
            r = r + self.dr
        return Point(r)
    def inverse(self):
        minc = copy.deepcopy(self)
        if self.q != None:
            minc.q = qconj(self.q)
            minc.dr = -qrotate(minc.q,self.dr)
        else:
            minc.dr = -self.dr
        return minc
    def composed(self, first_movement):
        mnew = copy.deepcopy(self)
        if self.q != None:
            if first_movement.q != None:
                mnew.q = qmul(self.q,first_movement.q)
            mnew.dr = self.dr + qrotate(self.q,first_movement.dr)
        elif first_movement.q != None:
            mnew.q = first_movement.q
            mnew.dr = self.dr + first_movement.dr
        return mnew
    def is_pure_translation(self):
        return self.q == None or abs(self.q[0] - 1) < 1e-12



