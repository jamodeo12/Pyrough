from setuptools import setup

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(name="pyrough",
version="1.0",
description="A tool for rough samples constructions",
author="Hugo Iteney",
author_email="hugo.iteney@im2np.fr",
packages=[],
install_requires=requirements
]
      )

#ATOMSK_PATH=''/usr/local/bin/atomsk''
#GMSH_PATH=/home/hiteney/.local/bin/gmsh