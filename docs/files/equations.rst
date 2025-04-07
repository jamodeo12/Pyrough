Rough surface generation
========================

The characterization of a rough surface can be done with respect to its spatial frequency content.
This can be transformed into a constructive method by using a sum of trigonometric functions similar to a Fourier series.
Each term in such a sum represents a certain frequency of spatial oscillation.


In *Pyrough*, a discrete set of spatial frequencies :math:`\nu_x = a` and :math:`\nu_y =  b` (where :math:`a` and :math:`b` are integers) is used to rationalize the range of investigated frequencies.
:math:`A` and :math:`B` are defined as the respective high-frequency cutoffs for :math:`a` and :math:`b` so that :math:`a \in [-A;A]` and :math:`b \in [-B;B]`.
Thus, the shortest wavelengths are :math:`\lambda_{x,min} = \frac{1}{A}` and :math:`\lambda_{y,min} = \frac{1}{B}` along :math:`\vec{u_x}` and :math:`\vec{u_y}` directions, respectively.
:math:`a` and :math:`b` can be positive or negative to ensure oscillations in both directions.
A rough surface :math:`h(x,y)` can be described by a sum of elementary waves as,

.. math::

	h(x,y) = \sum_{a=-A}^{A}  \sum_{b=-B}^{B} \alpha_{a,b}cos[2\pi(a x + b y) + \phi]

where :math:`\alpha_{a,b}` is the associated amplitude to each elementary wave.

Two more contributions are made in order to allow *Pyrough* to generate self-affine rough surfaces that are randomly perturbed.
First, the phase is randomly perturbed using :math:`\phi = U(a,b)` where :math:`U` states for a uniform distribution on an interval of length :math:`\pi`.
Also, random perturbations and self-affine aspects are implemented within :math:`\alpha_{a,b}`.
One can for example choose :math:`\alpha_{a,b}` as a zero-centered Gaussian distribution to get a smooth but random variation in amplitudes without constraining the magnitude i.e., :math:`\alpha_{a,b} = G(a,b){(a^2+b^2)}^{-(1+\eta)}` where :math:`G` states for a reduced centered normal law and :math:`{(a^2+b^2)}^{-(1+\eta)}` traduces the self-affine aspect of the surface.
Finally, the construction of a randomly perturbed self-affine surface can be modeled as,

.. math::

	h(x,y) = C_1\sum_{a=-A}^{A} \sum_{b=-B}^{B} G(a,b) (a^2+b^2)^{-(1+\eta)} cos[2\pi(a x + b y) + U(a,b) ]

in which :math:`\eta` allows to monitor the roughness degree and :math:`C_1` is a normalization factor introduced to fit the surface heights to the sample dimensions.

.. figure:: ./pictures/eta_influence.png
    :scale: 50%
    :align: center

    Random numerical self-affine rough surfaces for various Hurst exponents. a) :math:`\eta = 0.`, b) :math:`\eta = 0.25`, c) :math:`\eta = 0.50`, d) :math:`\eta = 0.75`, e) PSD of the various *h* profiles. For all calculations, A=50, B=50 and h(x,y) values are normalized between 1 and -1.
