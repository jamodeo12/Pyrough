#! /bin/bash
# J. Amodeo, H. Iteney, Fri Dec  1 11:08:39 CET 2023
# Here is an example of usage of the pyrough -surface option applied to a polystyrene particle characterized by AFM.
# This is the example discussed in the article Iteney et al., Comput. Phys. Commun. (2023) 108958

The experimental surface.png file is a screenshot of Figure in Yamamoto et al. J. Colloid Interface Sci. 292 (2005)
Details can be found in the original paper, size is 200 nm x 200 nm, min and max heights are -7.5 nm and +7.5 nm. 

The -surface of Pyrough will allow to derive the roughness parameters of the surface, including the RMS, H, A and B

python ../../Pyrough.py -surface surface.png 200 -7.5 7.5 

As outpus, you should obtain 3 plots of the original image, PSD plot and a digital twin (saved in Rough_data.csv), as well as a parameters you can use as inputs in Pyrough to design new digital twins.

#====== > Extraction of new inputs for Pyrough ...
#Here are the parameters you can use as inputs in Pyrough to generate other digital twins of your experimental image:
#RMS :  4.898416269531343
#H :  0.62
#A :  31.0
#B :  31.0

'N.B.:' Please note that -surface option better works with squared images. 
Rectangular a x b images will be made square by by Pyrough :) and only a min(a,b) x min(a,b) subimage will be processed.
