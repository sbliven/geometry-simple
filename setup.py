#!/usr/bin/env python

from distutils.core import setup

with open('README') as readme:
    long_description = readme.read()

setup(name='geometry-simple',
    version='0.2',
    description='3D geometry library for python',
    long_description=long_description,
    author='Spencer Bliven',
    author_email='spencer@bliven.us',
    url='https://github.com/sbliven/geometry-simple',
    download_url='https://github.com/sbliven/geometry-simple/dist/',
    py_modules=['geo'],
    requires=['numpy'],
    license='MIT',
)

