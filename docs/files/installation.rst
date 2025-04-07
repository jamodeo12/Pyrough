Installation
============

Required Python packages 
------------------------

* `NumPy <https://numpy.org/doc/stable/index.html>`_
* `Gmsh <https://gmsh.info/>`_
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

Run Pyrough with options
------------------------

Tu run image surface analysis and generate digital twin, add the *-surface* keyword, the image file to  analyze and the heights interval, when running Pyrough:

.. code-block:: python

    python Pyrough.py -surface Surf.png size zmin zmax

To test the installation of Pyrough, you can use the *-test_pyrough_execution* option which checks the correct computation for each .json file in the examples/ folder:

.. code-block:: python

    python Pyrough.py -test_pyrough_execution
