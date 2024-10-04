# Surface growth in tumour spheroids

This repository contains the source code used in the simulations described in [[1]](#1), along with experimental data extracted from the work in [[2]](#2).

## Dependencies

The code is written in Python 3 and tested with version 3.12. The following additional libraries are required:
* FEniCS (https://fenicsproject.org/)
* Numpy (https://numpy.org/)
* Scipy (https://scipy.org/)
* gmsh (https://gmsh.info/)
* BiFEniCS (https://github.com/riccobelli/bifenics)

## Repository structure

The repository is organized as follows:
* `2D_disk` contains the code for simulating the opening of a spheroid in 2D.
    * `disk.py` implements the class for simulating the opening of a cut disk.
    * `mesh.py` uses GMSH to generate the mesh and convert it into an xml format.
* `3D_sphere` contains the code for simulating the opening of a spheroid in 3D.
    * `sphere.py` implements the class for simulating the opening of a cut sphere.
    * `mesh.py` uses GMSH to generate the mesh and convert it into an xml format.
* `data` includes the data extracted from [[2]](#2) using a plot digitizer and reported in the paper [[1]](#1).
    * `data_plot_2D.csv` contains the data extracted from the bottom panel of Fig. 2a in [[2]](#2).
    * `data_plot_3D.csv` contains the data extracted from the top panel of Fig. 2a in [[2]](#2).

## Citing

If you find this code useful for your work, please cite [[1]](#1).

## License

The source code contained in this repository is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 2.1 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [GNU Lesser General Public License](https://www.gnu.org/licenses/lgpl-3.0.html) for more details.

## References
<a id="1">[1]</a>
D. Riccobelli, "Elastocapillarity-driven surface growth in tumour spheroids", *under review*.

<a id="2">[2]</a>
L. Guillaume et al. "Characterization of the physical properties of tumor-derived spheroids reveals critical insights for pre-clinical studies." *Scientific reports* 9.1 (2019): 6597.

## Author
Davide Riccobelli
