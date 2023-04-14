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

	git clone git@github.com:HugoIte/Pyrough.git

Run Pyrough from anywhere 
-------------------------

To run Pyrough in any folder, user can make it executable : 

.. code-block:: python

	chmod -x $PATH/Pyrough.py

And then create a symbolic link from *~/usr/bin* : 

.. code-block:: python

	ln -s $PATH/Pyrough.py pyrough

This command line will then run Pyrough from any folder :

.. code-block:: python

	pyrough input_file.json



