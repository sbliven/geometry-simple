#!/usr/bin/env python
from __future__ import with_statement
from distutils.core import setup
import sys

if sys.version_info < (2,6):
    print("Error: Python 2.6 is required.")
    sys.exit(1)

with open('README.md') as readme:
    long_description = readme.read()

name = 'geometry-simple'
version = '0.2'
setup(name=name,
    version=version,
    description='3D geometry library for python',
    long_description=long_description,
    author='Spencer Bliven',
    author_email='spencer@bliven.us',
    url='https://github.com/sbliven/geometry-simple',
    download_url="https://raw.github.com/sbliven/geometry-simple/master/dist/%s-%s.tar.gz"%(name,version),
    py_modules=['geo'],
    requires=['numpy'],
    license='MIT',
)

