com1DFA DFA-Kernel numerics
============================

Mesh and interpolation
-----------------------
The numerical method used in com1DFA mixes particle methods and
grid methods. Mass and momentum are tracked using particles but flow
depth is tracked using the grid. The grid is also to access topographic information
(surface elevation, normal) as well as for displaying results. Therefore it is
important to define properly the mesh and the interpolation method to

Mesh
~~~~~~

For practical reasons, a 2D rectilinear grid is used. Indeed the topographic
input information is read from 2D raster files which corresponds exactly to a
2D rectilinear grid. Moreover, as we will see in the following sections,
2D rectilinear grids are very convenient for interpolations as well as for
particle tracking. The 2D rectilinear grid is composed of :math:`N_{rows}` and
:math:`N_{cols}` rows and columns of square cells (of side length :math:`csz`)
as described in :numref:`rasterGrid`. Each cell has a center (also referred to
as node or grid point) and four vertices. Grid values are defined at cell centers.

.. _rasterGrid:

.. figure:: _static/rasterGrid.png
        :width: 90%

        Rectangular grid

Cell normals
""""""""""""""
There are many different methods available for computing normal vectors
on a 2D rectilinear grid. Several options are available in com1DFA.

The first one consists in computing the cross product of the diagonal vectors
between four cell centers. This defines the normal vector at the vertices. It is
then possible to interpolate the normal vector the cell centers from the ones
at the vertices.

The other methods use the plane defined by the different adjacent triangles to
a cell center. Each triangle has a normal and the cell normal is the average
of the triangles normals.

.. _meshNormal:

.. figure:: _static/meshNormal.png
        :width: 90%

        Rectangular grid

Cell area
"""""""""""



Interpolation
~~~~~~~~~~~~~~
In the DFA kernel, mass, flow depths, velocity fields can be defined at particle
location or on the grid. We need a method to be able to go from particle property
to grid field values and from grid values to particle property.

Grid to particle
"""""""""""""""""""

On a 2D rectilinear grid, scalar and vector fields defined on grid points
can be evaluated anywhere within the mesh using a bilinear interpolation
between grid points. Evaluating a vector field simply consists in evaluating
the three components as scalar fields.

The bilinear interpolation consists in successive linear interpolations
in both :math:`x` and :math:`y` using the four nearest grid points,
two linear interpolations in the first direction (in our case in the
:math:`y` direction in order to evaluated :math:`f_{0v}` and :math:`f_{1v}`)
followed by a second linear interpolation in the second direction
(:math:`x` in our case to finally evaluate :math:`f_{uv}`) as shown on :numref:`BilinearInterp`:

.. math::
    \begin{aligned}
    f_{0v} = & (1-v)f_{00} + vf_{01}\\
    f_{1v} = & (1-v)f_{10} + vf_{11}
    \end{aligned}

and

.. math::
    \begin{aligned}
    f_{uv} = & (1-u)f_{0v} + uf_{1v}\\
           = & (1-u)(1-v)f_{00} + (1-u)vf_{01} + u(1-v)f_{10} + uvf_{11}\\
                  = & w_{00}f_{00} + w_{01}f_{01} + w_{10}f_{10} + w_{11}f_{11}
    \end{aligned}

the :math:`w` are the bilinear weights.


.. _BilinearInterp:

.. figure:: _static/BilinearInterp.png
        :width: 90%

        Bilinear interpolation on in a unit cell.


Particles to grid
"""""""""""""""""""
Going from particle property to grid value is also based on bilinear interpolation and
weights but require a bit more care in order to conserve mass and momentum balance.
Flow depth and velocity fields are determined on the grid using, as intermediate step
mass and momentum fields. First, mass and momentum grid fields can be evaluated by
summing particles mass and momentum. This can be donne using the bilinear
weights :math:`w` defined in the previous paragraph (here :math:`f` represents
the mass or momentum):

.. math::
    \begin{aligned}
    f_{00} = & w_{00}f_{uv}\\
    f_{01} = & w_{01}f_{uv}\\
    f_{10} = & w_{10}f_{uv}\\
    f_{11} = & w_{11}f_{uv}
    \end{aligned}

This method ensures that the total mass and momentum of the particles is
preserved (the mass and momentum on the grid will sum up to the same total)
Flow depth and velocity grid fields can then be deduced from the mass and momentum
fields and the cell area (real area of each grid cell).

Neighbor search
------------------


SPH gadient
--------------


Artificial viscosity
--------------------
