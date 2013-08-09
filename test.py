#!/usr/bin/env python
"""Geometry-simple Tests

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

from geo import *
from math import *
import unittest


class GeoTest(unittest.TestCase):

    def assertListAlmostEqual(self,first, second, msg=None, places=7, delta=None):
        self.assertEqual(len(first),len(second),msg)
        for i in xrange(len(first)):
            self.assertAlmostEqual(first[i],second[i],places,msg,delta)

    def setUp(self):
        self.longMessage = True

        self.origin = Point(0,0,0)
        self.xaxis = Line(self.origin,Point(1,0,0))
        self.yaxis = Line(self.origin,Point(0,1,0))
        self.zaxis = Line(self.origin,Point(0,0,1))
        self.p1 = Point(1,1,1)
        self.p2 = Point(1,0,0)
        self.l = Line(self.p1,self.p2)
        self.xy_plane = Plane(self.origin,Point(1,0,0),Point(0,1,0))

    def test_Plane_construction(self):
        # y=1 plane
        p = Plane(Point(1,1,1), Point(2,1,0), Point(3,1,3))
        self.assertListEqual(list(p.r),[1,1,1],"Wrong origin")
        self.assertListEqual(list(p.n),[0,-1,0],"Wrong normal")
        # Opposite orientation
        p = Plane(Point(1,1,1), Point(3,1,3), Point(2,1,0))
        self.assertListEqual(list(p.r),[1,1,1],"Wrong origin")
        self.assertListEqual(list(p.n),[0,1,0],"Wrong normal")

        # Point and line
        p = Plane(Point(1,1,1), Line(Point(5,0,5),Point(5,2,5)))
        self.assertListAlmostEqual(list(p.r),[1,1,1],"Wrong origin")
        self.assertListAlmostEqual(list(p.n),[0,1,0],"Wrong normal")
        p = Plane(Line(Point(7,-2,5),Point(7,0,4)), Point(2,1,1) )
        self.assertListAlmostEqual(list(p.r),[2,1,1],"Wrong origin")
        self.assertListAlmostEqual(list(p.n),[0,2/sqrt(5),-1/sqrt(5)],"Wrong normal")

    def test_measurements(self):
        """ Test construction and measurements """

        self.assertAlmostEqual(self.p1.distance_to(self.p2),sqrt(2))

        self.assertAlmostEqual(self.l.distance_to(self.p1),0)

        self.assertAlmostEqual(self.p1.distance_to(self.xy_plane),1)

        self.assertAlmostEqual(self.l.angle_to(self.xy_plane),pi/4)

    def test_lfit(self):
        """ Make sure that the fitted line coincides with l for three points on l. """
        Lfit = Line(self.p1,self.p2,Point(1,0.5,0.5))
        self.assertAlmostEqual(self.l.angle_to(Lfit),0)
        self.assertAlmostEqual(self.l.distance_to(Lfit),0)

        """ Midpoints """
        p = Point(-1,-2,-3).midpoint_to(Point(1,2,3))
        self.assertAlmostEqual(p.distance_to(self.origin),0)

        p = self.xaxis.midpoint_to(self.yaxis)
        self.assertAlmostEqual(p.distance_to(self.origin),0)

        p = self.zaxis.midpoint_to(self.xy_plane)
        self.assertAlmostEqual(p.distance_to(self.origin),0)



    def test_movements(self):
        """ test movements """

        # pure translation
        mt = Movement(self.origin,self.p1)
        self.assertAlmostEqual(self.p1.moved(mt).distance_to(self.origin),2*self.p1.distance_to(self.origin))

        #pure rotation
        mr = Movement(self.xaxis,self.yaxis) # rotate self.xaxis onto self.yaxis
        self.assertAlmostEqual(self.xaxis.moved(mr).distance_to(self.yaxis),0)
        self.assertAlmostEqual(self.origin.moved(mr).distance_to(self.origin),0)

    def test_coef(self):
        """ Test the coef method """
        p = Plane(self.origin,Point(1,0,0),Point(0,1,0))
        r = (0,0,1,0)
        c = p.coef(100)
        self.assertListAlmostEqual(c,r,"Expected coef {0}, got {1} for {2}".format(r,c,p))

        p = Plane(Point(1,1,1),Point(5,1,0),Point(0,1,-8))
        r = (0,1,0,1)
        c = p.coef(100)
        self.assertListAlmostEqual(c,r,"Expected coef {0}, got {1} for {2}".format(r,c,p))

        p = Plane(Point(-5,1,1),Point(-5,1,0),Point(-5.00000000,2,-8))
        r = (1,0,0,-5)
        c = p.coef(100)
        self.assertListAlmostEqual(c,r,"Expected coef {0}, got {1} for {2}".format(r,c,p))

        p = Plane(Point(1,1,1),Point(1,1,0),Point(2,0,1))
        r = (1,1,0,2)
        c = p.coef(100)
        self.assertListAlmostEqual(c,r,"Expected coef {0}, got {1} for {2}".format(r,c,p))
        p = Plane(Point(0,0,0),Point(-1,1,0),Point(0,-1,1))
        r = (1,1,1,0)
        c = p.coef(100)
        self.assertListAlmostEqual(c,r,"Expected coef {0}, got {1} for {2}".format(r,c,p))
        #opposite orientation
        p = Plane(Point(0,0,0),Point(0,-1,1),Point(-1,1,0))
        r = (1,1,1,0)
        c = p.coef(100)
        self.assertListAlmostEqual(c,r,"Expected coef {0}, got {1} for {2}".format(r,c,p))

if __name__ == '__main__':
    unittest.main()
