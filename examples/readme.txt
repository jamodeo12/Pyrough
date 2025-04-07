#! /bin/bash
# J. Amodeo, Thu Apr  3 14:44:37 CEST 2025
# example file description

box.json            : typical example of Box keyword for the design of a  block with a single rough surface, here applied FCC Au with <001> cubic orientation
box_110.json        : same, but with [110], [-110] and [001] orientation
box_Mg_****.json    : same but applied to HCP Mg w/ various orientation

cube.json           : typical example of Cube keyword for the design of a cube with all surfaces made rough, here applied to FCC Au with <001> cubic orientation 

cwire.json          : typical example of cWire keyword for the design of a rough wire w/ a circular cross section, here applied to FCC Au
fWire.json          : typical example of fWire keyword for the design of a rough wire w/ a facetted cross section, here applied to FCC Au

grain.json          : typical example of Grain keyword to design rough grain boundaries, here applied to FCC/FCC <001>/<001> interface
grain_fcc100_bcc100.json  : same, but applied to FCC/BCC <001>/<001> interface

sphere.json         :  typical example of sphere keyword, to design a rough sphere using spherical harmonics. Here applied to FCC Cu

wulff.json          : typical example of Wulff keyword, to design a Wulff-shaped particle with rough surfaces, here applied to FCC Au

precinmatrix_110bcc_in_100fcc.json      : typical example of precinmatrix keyword to design a nanostructured sample (matrix + precipitate ourgh or not) at the atomic scale. Here is a <110>-oriented BCC Fe cube inserted into a an <100>-oriented FCC Ag matrix
precinmatrix_bcc_sphere_in_hcp.json     : same, but inserting a (rough or not)  BCC Fe sphere into an HCP Mg matrix

surface/            : this folder contains an example of the -surface option application
