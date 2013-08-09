===============
geometry-simple
===============

3D geometry library for python. The goal is to make geometric operations simple
without a deep knowledge of linear algebra.


Installation
------------

Download the latest version of the package from
https://raw.github.com/sbliven/geometry-simple/master/dist/geometry-simple-0.2.tar.gz

In the resulting package, run

    sudo python setup.py install

This will install the package for all users.


Dependencies
------------

The library requires numpy (http://www.numpy.org) for the underlying linear
algebra.


Examples
--------

Import the library::

    from geo import *

Create some geometric objects::

    origin = Point(0,0,0)
    xline = Line( origin, Point(10,0,0) )
    yzplane = Plane( origin, Point(0,1,0), Point(0,0,1) )

Calculate distances and angles::

    dist = Point(5,5,5).distance_to(yzplane)
    use_degrees()
    ang = xline.angle_to( yzplane )

Use linear regression to fit lines and planes::

    from random import gauss
    line_points = [ Point( gauss(x,.1), gauss(2*x,.1), 0) for x in xrange(10) ] #points around y=2x
    line = Line(line_points)
    
    plane_points = [ Point( gauss(0,10), gauss(0,10), gauss(0,.1)) for x in xrange(20) ] #points in the xy plane
    plane = Plane(plane_points)
    
Easily apply affine transforms::

    translation = Movement(origin, Point(5,-5,1))
    Point(1,1,1).moved(translation)
    
    rot90 = Movement(xline, Line(origin, Point(0,1,0)) )
    Point(1,5,0).moved(rot90)


History
-------

This library was created by uekstrom in 2008. The original code can be found at
http://code.google.com/p/geometry-simple/

The code was expanded and moved to github by Spencer Bliven in 2013.

