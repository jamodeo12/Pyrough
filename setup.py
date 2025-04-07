from setuptools import find_packages, setup

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="Pyrough",
    version="1.2",
    description="a tool to build 3D samples with rough surfaces for atomistic and finite-element simulations",
    author="Jonathan Amodeo et al.",
    url="https://github.com/jamodeo12/Pyrough",
    author_email="jonathan.amodeo@cnrs.fr",
    packages=find_packages(),
    install_requires=requirements,
)
