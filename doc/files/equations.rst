Rough surface generation
========================

The characterization of a rough surface can be done with respect to its spatial frequency content. This can be transformed into a constructive method by using a sum of trigonometric functions similar to a Fourier series. Each term in such a sum represents a certain frequency of spatial oscillation.


Spatial frequencies
-------------------

The mathematical expression of oscillations over time is given by :math:`cos(2\pi f t)` with :math:`f` the time frenquency in Hz.
In the case of spatial oscilations, we simply have to replace time variables by spatial ones: 

.. math::

	cos(2 \pi \nu x)

with :math:`\nu` the spatial frequency in :math:`m^{-1}`.

In the case of 2D and Cartesian coordinates, we have:

.. math::

	cos(2 \pi ( \nu_x x + \nu_y y)) = cos(\vec{k}.\vec{x})

where :math:`\vec{k} = (k_x,k_y) = (2\pi \nu_x,2\pi \nu_y)` is the wave vector representing the direction of the wave and :math:`\vec{x} = (x,y)` the position vector.

Elementary waves
----------------

A rough surface :math:`z(x,y)` is composed of several elementary waves: 

.. math::

	cos(\vec{k}.\vec{x} + \phi) 

where :math:`\phi` is a phase angle. 

Due to the relationship :math:`sin(\theta) = cos(\frac{\pi}{2} - \theta)`, there is no influence of using either sines or cosines.

For a random surface construction, the phase angle :math:`\phi` can be any value in the interval :math:`[0;\pi]` or :math:`[-\frac{\pi}{2};\frac{\pi}{2}]`. 
It follows then a uniform random distribution in such an interval of length :math:`\pi`.
It is important to notice that there may be end-point or wrap-around effects if we choose an interval with a size bigger than :math:`\pi`.
This is due to the cosine function being its own mirror image in steps of :math:`\pi`, according to :math:`cos(\pi - \theta) = -cos(\theta)`.

In order to include this formulation in computational simulations, a discrete set of spatial frequencies is needed : :math:`\nu_x = m` and :math:`\nu_y = n`, where :math:`m` and :math:`n` are integers.

We define :math:`M` and :math:`N` as respective high-frequency cutoffs for :math:`m` and :math:`n` so that :math:`m \in ⟦-M;M⟧` and :math:`n \in ⟦-N;N⟧`. This means that the shortest wavelengths being represented are respectively :math:`\lambda_{x,min} = \frac{1}{M}` and :math:`\lambda_{y,min} = \frac{1}{N}` for :math:`x` and :math:`y` directions.

A rough surface is then constituted of elementary waves of the form:

.. math::

	cos( \vec{k_{m,n}}.\vec{x} + \phi) = cos(2\pi (mx + ny) + \phi)

To get a method of synthesizing a surface with no preferred direction of oscillation, we let :math:`m` and :math:`n` take positive and negative values with equal probabilities

Associated Amplitudes
---------------------

Each elementary wave has an associated amplitude so that each wave component has the form :math:`A_{m,n} cos( \vec{k_{m,n}}.\vec{x} + \phi)`.

The final surface is a sum over the wave components:

.. math ::
	
	z(\vec{x}) = \sum_{m,n} A_{m,n} cos(\vec{k_{m,n}}.\vec{x} + \phi)

In order to obtain a more realistic rough surface, two more contributions should be made :

	* Since :math:`\phi` can take any value in an interval of length :math:`\pi` with the same probability, :math:`\phi = U(m,n)` where :math:`U` states for a random function with a uniform distribution. 
	* Autocorrelation is implemented through :math:`A_{m,n}` which follows a Gaussian law. A Gaussian, or normal, distribution is chosen to get a smooth but random variation in amplitudes with no limit on the magnitude. Autocorrelation is important to modulate the Gaussian distribution and get credible roughness. It avoids or limits the presence of too rough peaks and valleys. Therefore, :math:`A_{m,n} = G(m,n).(m^2 + n^2)^{-\frac{b}{2}}` where :math:`G` states for a random function with a Gaussian distribution and :math:`(m^2 + n^2)^{-\frac{b}{2}}` is the autocorrelation function.

To sum up, we can represent rough surfaces by the final equation: 

.. math::

	z(x,y) = C_1 \sum^{M}_{m = -M} \sum^{N}_{n = -N} G(m,n).(m^2 + n^2)^{-\frac{b}{2}} cos[2 \pi (mx + ny) + U(m,n)]

Where :math:`C_1` is a normalization factor and :math:`b` is the spectral exponent that parameterizes how fast the high frequencies are attenuated.

.. figure:: ./pictures/Influence_of_b.png
    :align: center

    Surfaces generated on the square :math:`[0,1]×[0,1]` by superposing 50 frequencies with spectral exponents :math:`b = 1.5`, :math:`b = 1.8`, :math:`b = 2` and :math:`b = 2.2`, clockwise from the top-left image.



Periodicity
-----------

Due to its definition, the function :math:`z(x,y)` is periodic:

.. math::

	z(x+1,y) = C_1 \sum^{M}_{m = -M} \sum^{N}_{n = -N} G(m,n).(m^2 + n^2)^{-\frac{b}{2}} cos[2 \pi (mx + ny) + U(m,n) + 2m\pi] = z(x,y)

The same demonstration for :math:`z(x,y+1)` justifies the periodicity of the function on the interval :math:`[0;1]×[0;1]`.
The surface is generated over a rectangle :math:`[a; a + 1]×[b; b + 1]` or smaller in order to avoid spatial periodicity and pattern repetitions.
