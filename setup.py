#!/usr/bin/env python
from setuptools import setup, find_packages

# ensure you are using the correct version of setuptools: pip install setuptools==58.2.0

setup(
   description='Drag models',
   package_dir={"": "."},
   packages=find_packages(
       where=".",
   ),
   #scripts=['drag_models/plot_drag_models.py'],
)