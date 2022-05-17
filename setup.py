from setuptools import setup

setup(name="pyrough",
version="1.0",
description="A tool for rough samples constructions",
author="Hugo Iteney",
author_email="hugo.iteney@im2np.fr",
packages=[],
install_requires=[
    'numpy>=1.19.5',
    'pygmsh>=7.1.12',
    'meshio>=4.4.6',
    'wulffpack>=1.1',
    'ase>=3.21.1'
]
      )