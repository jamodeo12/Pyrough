Installation
============

Required Python packages 
------------------------

* `NumPy <https://numpy.org/doc/stable/index.html>`_
* `pygmsh <https://pypi.org/project/pygmsh/>`_
* `meshio <https://pypi.org/project/meshio/>`_
* `ase <https://wiki.fysik.dtu.dk/ase/index.html>`_
* `WulffPack <https://wulffpack.materialsmodeling.org/>`_

Additionnal requirement 
------------------------

* `Atomsk <https://atomsk.univ-lille.fr/fr/>`_

Download the code 
-----------------

.. code-block:: python

	git clone git@github.com:jamodeo12/Pyrough.git

Run Pyrough
-----------

This command line will enable to run Pyrough for the generation of 3D rough shape:

.. code-block:: python

	python Pyrough.py path/to/the/file/input_file.json

Tu run the surface analysis process, simply add the *-surface* option, then the surface to analyse and the heights interval:

.. code-block:: python

	python Pyrough.py -surface Surf.png zmin zmax



