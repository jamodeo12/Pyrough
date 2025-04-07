Generation of rough shapes
=============================

In this chapter, the strategy for applying roughness to simple shapes is explained.

Box
----
A box of dimensions *width*, *length* and *height* is first constructed.
Surface nodes are the ones at the maximum height of the pristine sample mesh.
The roughening procedure consists in displacing each surface node along the *z*-axis according to *h(x,y)*.

.. literalinclude:: ../../examples/box.json
   :language: python

Sphere
------
Surface roughness is applied to each surface node of the sphere by summing spherical harmonics.

.. literalinclude:: ../../examples/sphere.json
   :language: python

Wire
----
A cylinder of radius *radius* and length *length* is constructed.
Surface nodes are those positioned on the external radius *r* of the cylinder.
Thus, each surface node from the wire mesh is moved along *r*-axis by its corresponding displacement vector from the generated rough surface.

.. literalinclude:: ../../examples/cwire.json
   :language: python

Faceted wire
------------
The initial mesh is constructed from centered regular polygon base with *nfaces* sides and of length *length*.
Surface nodes are here translated along the facet normal :math:`\vec{n}` they belong to.

.. literalinclude:: ../../examples/fwire.json
   :language: python

Wulff-shaped faceted particle
---------------------------------
One rough surface per facet is generated in this case.
Each node is translated along its facet normal :math:`\vec{n}` using *h(x,y)* to build the rough surface.
The process is repeated for each facet of the sample.

.. literalinclude:: ../../examples/wulff.json
   :language: python

Cube
----
A pristine cube of edge length *length* is generated.
Surface roughness is applied in the same way as in the case of wulff shapes.

.. literalinclude:: ../../examples/cube.json
   :language: python

Grain boundary
--------------
A 3D block made composed of two grains separated by a rough and wavy grain boundary, where each grain can have its own crystal structure and orientation.

.. literalinclude:: ../../examples/grain.json
   :language: python

Precipitate embedded in a matrix
--------------------------------
A precipitate (rough or not) embedded into a 3D-periodic matrix.
The precipitate can take any of the shapes already managed by Pyrough (sphere, cube, cylinder, wulff-shape, *etc.*).

.. literalinclude:: ../../examples/precinmatrix_110bcc_in_100fcc.json
   :language: python


