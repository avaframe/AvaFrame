com1DFA DFA-Kernel numerics
============================

In :ref:`theoryCom1DFA`, the mass and momentum equation were derived using a Lagrangian approach.
In order to solve this set of equations numerically, we employ a mixed of particle and grid approach.
We discretize the material into particles and solve the momentum equation for these particles
(:cite:`Sa2007,SaGr2009`).
The pressure gradients are computed using a SPH method
(**S**\ moothed **P**\ article **H**\ ydrodynamics :cite:`Mo1992`). Some of the
advantages of the SPH method are that free surface flows, material boundaries and
moving boundary conditions are considered implicitly. In addition, large
deformations can be modeled due to the fact that the method is not based
on a mesh.
However, we use a grid to compute several parameters that are required for the computations as
for example surface normal vectors and flow thickness.


Numerical discretization: a particle grid approach
-----------------------------------------------------

Space discretization: The particular momentum equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Discretizing the material into particles (particle quantities are denoted by the subscript :math:`k`, e.g.
:math:`m_k = \rho_0 V_k` is the mass of particle :math:`k`) leads to the following continuity equation:

.. math::
	\frac{\mathrm{d}}{\mathrm{d}t} m_k = A_k^\text{ent} q^{\text{ent}}
  :label: eq-particle-continuity-balance

By assuming that the Lagrangian control volume :math:`V` can be represented by a particle,
we can derive from :eq:`momentum-balance6` to the particular momentum equation in the normal direction and in the tangent plane
(the entrainment force is detailed in Appendix~\ref{ap-entrainment-force}):

.. math::
  \left\{
  \begin{aligned}
    & p_k^b = \rho_0 \, h_k \, g_k^\text{eff}\\
    &m_k \frac{\mathrm{d}\overline{\mathbf{u}}_k}{\mathrm{d}t} = A_k^b\boldsymbol{\tau^b}
  	- m_k \, g^\text{eff}_k \, \boldsymbol{\nabla}_{s} h	+ m_k \mathbf{g}_s	+ \mathbf{F}_k^{\text{ext}}
    - \overline{\mathbf{u}}_k A_k^\text{ent} q_k^{\text{ent}}
    - m_k \left( \overline{\mathbf{u}}_k \cdot \frac{\mathrm{d}\mathbf{v}_{3,k}}{\mathrm{d}t} \right)\mathbf{v}_{3,k}
  \end{aligned}
  \right.
  :label: eq-momentum-particle

In this equation (:eq:`eq-momentum-particle`), flow thickness gradient, basal friction and
curvature acceleration terms need to be further developed and discretized.

Time discretization
~~~~~~~~~~~~~~~~~~~~~

The momentum equation is solved numerically in time using an Euler time scheme.
The time derivative of any quantity :math:`f` is approximated by:

.. math::
  \frac{\mathrm{d}f_k}{\mathrm{d}t} \approx
  \frac{f_k^{n+1} - f_k^n}{\Delta t}

where :math:`\Delta t` represents the time step and :math:`f^n = f(t^n)`, :math:`t^n = n \Delta t`.
For the velocity this reads:

.. math::
  \frac{\mathrm{d}\overline{\mathbf{u}}_k}{\mathrm{d}t} \approx
  \frac{\overline{\mathbf{u}}_k^{n+1} - \overline{\mathbf{u}}_k^n}{\Delta t}

The position of the particles is then updated using a centered Euler scheme:

.. math::
  \mathbf{x}_{k}^{n+1} = \mathbf{x}_{k}^{n} + \frac{\Delta t}{2m}\left(\overline{\mathbf{u}}^{n+1}_{k} + \overline{\mathbf{u}}^{n}_{k}\right)


Taking the forces into account is done in two subsequent steps as forces acting on the particles can be
sorted into driving forces and friction forces.
Friction forces act against the particle motion only affecting the magnitude of the velocity.
They can in no case become driving forces..
This is why in a first step the velocity is updated with the driving forces before updating in a
second step the velocity magnitude applying the friction force.

Grid and interpolation
-----------------------

Topography information is usually provided in a raster format which corresponds to a regular rectilinear
mesh, from hereon referred to as grid.
In order to get information on the surface elevation and normal vectors, the topography information
needs to be interpolated at the particle locations, and this needs to be repeated at each time step
since the particles are moving.
Similarly, the particle properties such as mass or momentum, which translate into flow thickness and
velocity, also need to be interpolated onto the grid.
Grid velocity is needed to compute the artificial viscosity term, ensuring numerical stability,
see Sect.~\ref{sec-numerical-stability}.
Grid flow thickness is used to compute the particle flow thickness which is required for computing
the friction force.
Due to the utilized regular rectilinear mesh a bilinear interpolation method is applied for its
simplicity and suitability.
It also ensures the conservation of mass or momentum when interpolating from particles to grid and back.

Here is a description of the grid and the interpolation method that is used to
switch from particle to grid values and the other way around.

Grid
~~~~~~

For practical reasons, a 2D rectilinear mesh (grid) is used. Indeed the topographic
input information is read from 2D raster files (with :math:`N_{y}` and :math:`N_{x}`
rows and columns) which correspond exactly to a
2D rectilinear mesh. Moreover, as we will see in the following sections,
2D rectilinear meshes are very convenient for interpolations as well as for
particle tracking. The grid is composed of :math:`N_{y}` and
:math:`N_{x}` rows and columns of square cells (of side length :math:`csz`)
and :math:`N_{y}+1` and :math:`N_{x}+1` rows and columns of vertices
as described in :numref:`rasterGrid`. Each cell has a center and four vertices.
The data read from the raster file is assigned to the cell centers. Note that
although this is a 2D grid, as we use a terrain-following coordinate system to perform
our computations, this 2D grid is oriented in 3D space and hence the projected side length
corresponds to :math:`csz`, whereas the actual side length and hence also the
:ref:`DFAnumerics:cell area`, depend on the local slope,
expressed by the :ref:`DFAnumerics:Cell normals`.

.. _rasterGrid:

.. figure:: _static/rasterGrid0.png
        :width: 90%

        Rectangular grid

Cell normals
""""""""""""""
There are many different methods available for computing normal vectors
on a grid. Several options are available in com1DFA.

The first one consists in computing the cross product of the diagonal vectors
between four cell centers. This defines the normal vector at the vertices. It is
then possible to interpolate the normal vector at the cell centers from the ones
at the vertices.

The other methods use the plane defined by different adjacent triangles to
a cell center. Each triangle has a normal and the cell center normal is the average
of the triangles normal vectors.

.. _meshNormal:

.. figure:: _static/meshNormal0.png
        :width: 90%

        Grid normal computation

Cell area
"""""""""""
The cell area can be deduced from the grid cellsize and the cell normal.
A cell is a plane (:math:`z = ax+by+c`) of same normal as the cell center:

.. math::
   \mathbf{n} = \frac{1}{\sqrt{1+a^2+b^2}}
   \left|\begin{aligned}
   &-a\\
   &-b\\
   &1
   \end{aligned}
   \right.

Surface integration over the cell extent leads to the area of the cell:

.. math::
   A_{cell} = \iint_{S} \mathrm{d}{S} = \int\limits_{0}^{csz}\int\limits_{0}^{csz}
   \sqrt{1+\frac{\partial z}{\partial x}^2+\frac{\partial z}{\partial y}^2}
   \mathrm{d}{x}\,\mathrm{d}{y} =
   csz^2 \sqrt{1+\frac{\partial z}{\partial x}^2+\frac{\partial z}{\partial y}^2} = \frac{csz^2}{n_z}


Interpolation
~~~~~~~~~~~~~~
In the DFA kernel, mass, flow thickness and flow velocity can be defined at particle
location or on the grid. We need a method to be able to go from particle properties
to grid (field) values and from grid values to particle properties.

Grid to particle
""""""""""""""""""

On a grid, scalar and vector fields defined at cell centers
can be evaluated anywhere within the grid using a bilinear interpolation
between grid cell centers. Evaluating a vector field simply consists in evaluating
the three components as scalar fields.

The bilinear interpolation consists in successive linear interpolations
in both :math:`x` and :math:`y` direction using the four nearest cell centers,
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

the :math:`w` are the bilinear weights. The example given here is for a unit cell.
For no unit cells, the :math:`u` and :math:`v` simply have to be normalized by the
cell size.


.. _BilinearInterp:

.. figure:: _static/BilinearInterp.png
        :width: 90%

        Bilinear interpolation in a unit grid (cell size is 1).


Particles to grid
"""""""""""""""""""
Going from particle property to grid value is also based on bilinear interpolation and
weights but requires a bit more care in order to conserve mass and momentum balance.
Flow thickness and velocity fields are determined on the grid using, as intermediate step,
mass and momentum fields. First, mass and momentum grid fields can be evaluated by
summing particles mass and momentum. This can be donne using the bilinear
weights :math:`w` defined in the previous paragraph (here :math:`f` represents
the mass or momentum and :math:`f_{uv}` is the particle value. :math:`f_{nm}`
, :math:`{n, m} \in \{0, 1\} \times \{0, 1\}`, are the cell center values):

.. math::
    \begin{aligned}
    f_{00} = & w_{00}f_{uv}\\
    f_{01} = & w_{01}f_{uv}\\
    f_{10} = & w_{10}f_{uv}\\
    f_{11} = & w_{11}f_{uv}
    \end{aligned}

The contribution of each particle to the different grid points is summed up to
finally give the grid value. This method ensures that the total mass and
momentum of the particles is preserved (the mass and momentum on the grid will
sum up to the same total). Flow thickness and velocity grid fields can then be deduced
from the mass and momentum fields and the cell area (actual area of each grid cell,
not the projected area).


Flow thickness and its gradient
----------------------------------

SPH method can be used to solve thickness integrated equations where a 2D
(respectively 3D) equation is reduced to a 1D (respectively 2D) one.
This is used in ocean engineering to solve shallow water equations (SWE)
in open or closed channels for example. In all these applications,
whether it is 1D or 2D SPH, the fluid is most of the time,
assumed to move on a horizontal plane (bed elevation is set to a constant).
In the case of avalanche flow, the "bed" is sloped and irregular.
The aim is to adapt the SPH method to apply it to thickness integrated equations
on a 2D surface living in a 3D world.

Flow thickness gradient computation using SPH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to assess the flow thickness gradient, we employ a SPH method (Smoothed Particles Hydrodynamics Method
:cite:`LiLi2010`), where the gradient is directly derived from the particles and does not require any mesh.
In contrast, a mesh method or a MPM (Material Point Method) would directly use a mesh formulation to
approximate the gradient or interpolate the particles properties on an underlying mesh and
then approximate the gradient of the flow thickness using a mesh formulation.

In theory, a SPH method does not require any mesh to compute the gradient.
However, applying this method requires finding neighbor particles.
This process can be sped up with the help of an underlying grid,  different neighbor search methods
are presented in :cite:`IhOrSoKoTe2014`, a "uniform grid method" is used in this paper.

The SPH method is used to express a quantity (the flow thickness in our case) and its gradient at
a certain particle location as a weighted sum of its neighbors properties.
The principle of the method is well described in :cite:`LiLi2010` and the basic formula reads:

.. math::
  \begin{aligned}
  f_{k} \simeq \langle f_{k}\rangle &= \sum\limits_{l}f_{l}A_{l}\,W_{kl}\\
  \boldsymbol{\nabla} f_{k} \simeq \langle \boldsymbol{\nabla} f_{k}\rangle &= -\sum\limits_{l}f_{l}A_{l}\,\boldsymbol{\nabla} W_{kl}
	\end{aligned}
  :label: eq-sph-formulation

Where :math:`W` represents the SPH-Kernel function (we employ the spiky kernel, see
:eq:`eq-kernel-function`) an the subscript :math:`l` denotes the neighbor particles to
particle :math:`k`.
This kernel function is designed to satisfy the unity condition, be an
approximation of the Dirac function and have a compact support domain
(:cite:`LiLi2010`).

:eq:`eq-sph-formulation` gives for the flow thickness:

.. math::
  h_{k}  \simeq \langle h_{k}\rangle &= \frac{1}{\rho_0}\,\sum\limits_{l}{m_{l}}\,W_{kl}\\
  \boldsymbol{\nabla}h_{k} \simeq \langle \boldsymbol{\nabla} h_{k}\rangle &= -\frac{1}{\rho_0}\,\sum\limits_{l}{m_{l}}\,\boldsymbol{\nabla}W_{kl}
  :label: sph formulation for fd

This method is usually either used in a 3D space where particles move
freely in this space and where the weighting factor for the summation is
the volume of the particle or on a 2D horizontal plane where the weighting
factor for the summation is the area of the particle and the gradient is
2D.
Here we want to compute the gradient of the flow thickness on a 2D surface
(the topography) that lives in 3D. The method used is analog to the SPH
gradient computation on the 2D horizontal plane but the gradient is 3D
and tangent to the surface (colinear to the local tangent plane).
The theoretical derivation in the following section shows that the SPH
computation is equivalent in applying the 2D SPH method in the local
tangent plane instead of in the horizontal plane.

.. _standard-method:

Standard method
""""""""""""""""

Let us start with the computation of the gradient of a scalar function
:math:`f \colon \mathbb{R}^2 \to \mathbb{R}` on a horizontal plane.
Let :math:`P_k=\mathbf{x}_k=(x_{k,1},x_{k,2})` and :math:`Q_l=\mathbf{x}_l=(x_{l,1},x_{l,2})` be
two points in :math:`\mathbb{R}^2` defined by their coordinates in the Cartesian coordinate system
:math:`(P_k,\mathbf{e_1},\mathbf{e_2})`. :math:`\mathbf{r}_{kl}=\mathbf{x}_k-\mathbf{x}_l` is the
vector going from :math:`Q_l` to :math:`P_k` and :math:`r_{kl} = \left\Vert \mathbf{r}_{kl}\right\Vert`
the length of this vector. Now consider the kernel function :math:`W`:

.. math::
  \left.
  \begin{aligned}
  W \colon \mathbb{R}^2 \times \mathbb{R}^2 \times \mathbb{R} &\to \mathbb{R}\\
  (P_k, Q_l, r_0) &\mapsto W(P_k, Q_l, r_0)
  \end{aligned}
  \right.\quad, r_0\in\mathbb{R} \mbox{ is the smoothing kernel length}

In the case of the spiky kernel, :math:`W` reads (2D case):

.. math::
   \begin{aligned}
   W_{kl} = &W(\mathbf{x_k},\mathbf{x_l},r_0) = W(\mathbf{x_k}-\mathbf{x_l},r_0) = W(\mathbf{r_{kl}},r_0)\\
   =&\frac{10}{\pi r_0^5}\left\{
   \begin{aligned}
   & (r_0 - \left\Vert \mathbf{r_{kl}}\right\Vert)^3, \quad &0\leq \left\Vert \mathbf{r_{kl}}\right\Vert \leq  r_0\\
   & 0 , & r_0 <\left\Vert \mathbf{r_{kl}}\right\Vert
   \end{aligned}
   \right.
   \end{aligned}
   :label: kernel function


:math:`\left\Vert \mathbf{r_{kl}}\right\Vert= \left\Vert \mathbf{x_{k}}-\mathbf{x_{l}}\right\Vert`
represents the distance between particle :math:`k` and :math:`l` and
:math:`r_0` the smoothing length.

Using the chain rule to express the gradient of :math:`W` in the Cartesian
coordinate system :math:`(x_1,x_2)` leads to:


.. math::
   \boldsymbol{\nabla}W_{kl} = \frac{\partial W}{\partial r}.\boldsymbol{\nabla}r,
   \quad r = \left\Vert \mathbf{r} \right\Vert = \sqrt{(x_{k,1}-x_{l,1})^2 + (x_{k,2}-x_{l,2})^2}
   :label: kernel function gradient 1

with,

.. math::
  \frac{\partial W}{\partial r} = -3\frac{10}{\pi r_0^5}\left\{
  \begin{aligned}
  & (r_0 - \left\Vert \mathbf{r_{kl}}\right\Vert)^2, \quad &0\leq \left\Vert \mathbf{r_{kl}}\right\Vert \leq  r_0\\
  & 0 , & r_0 <\left\Vert \mathbf{r_{kl}}\right\Vert
  \end{aligned}
  \right.

and

.. math::
  \frac{\partial r}{\partial x_{k,i}} = \frac{(x_{k,i}-x_{l,i})}{\sqrt{(x_{k,1}-x_{l,1})^2 + (x_{k,2}-x_{l,2})^2}},
  \quad i=\{1,2\}
which leads to the following expression for the gradient:

.. math::
   \boldsymbol{\nabla}W_{kl} = -3\frac{10}{\pi r_0^5}\left\{
   \begin{aligned}
   & (r_0 - \left\Vert \mathbf{r_{kl}}\right\Vert)^2\frac{\mathbf{r_{kl}}}{r_{kl}}, \quad &0\leq \left\Vert \mathbf{r_{kl}}\right\Vert \leq  r_0\\
   & 0 , & r_0 <\left\Vert \mathbf{r_{kl}}\right\Vert
   \end{aligned}
   \right.
   :label: kernel function gradient

The gradient of :math:`f` is then simply:

.. math::
    \boldsymbol{\nabla}f_{k} = -\sum\limits_{l}f_{l}A_{l}\,\boldsymbol{\nabla}W_{kl}
    :label: sph gradient

2.5D SPH method
""""""""""""""""
We now want to express a function :math:`f` and its gradient on a potentially
curved surface and express this gradient in the 3 dimensional Cartesian
coordinate system :math:`(P_k,\mathbf{e_1},\mathbf{e_2},\mathbf{e_3})`.

Let us consider a smooth surface :math:`\mathcal{S}` and two points
:math:`P_k=\mathbf{x}_k=(x_{k,1},x_{k,2},x_{k,3})` and :math:`Q_l=\mathbf{x}_l=(x_{l,1},x_{l,2},x_{l,3})`
on :math:`\mathcal{S}`. We can define :math:`\mathcal{TP}` the tangent plane
to :math:`\mathcal{S}` in :math:`P_k`. If :math:`\mathbf{u}_k` is the (none zero)
velocity of the particle at :math:`P_k`, it is possible to define the local
orthonormal coordinate system :math:`(P_k,\mathbf{V_1},\mathbf{V_2},\mathbf{V_3}=\mathbf{n})`
with :math:`\mathbf{V_1}=\frac{\mathbf{u}_k}{\left\Vert \mathbf{u}_k\right\Vert}`
and :math:`\mathbf{n}` the normal to :math:`\mathcal{S}` at :math:`P_k`.
Locally, :math:`\mathcal{S}` can be assimilated to :math:`\mathcal{TP}` and
:math:`Q_l` to its projection :math:`Q'_l` on :math:`\mathcal{TP}`.
The vector :math:`\mathbf{r'}_{kl}=\mathbf{x}_k-\mathbf{x'}_l` going from
:math:`Q'_l` to :math:`P_k` lies in :math:`\mathcal{TP}` and can be express
in the plane local basis:

.. math::
  \mathbf{r'}_{kl}=\mathbf{x}_k-\mathbf{x'}_l = v_{kl,1}\mathbf{V_1} + v_{kl,2}\mathbf{V_2}

It is important to define :math:`f` properly and the gradient that will be calculated:

.. math::
  \left.
  \begin{aligned}
  f \colon \mathcal{TP}\subset\mathbb{R}^3 &\to \mathbb{R}\\
  (x_1,x_2,x_3) &\mapsto f(x_1,x_2,x_3) = f(x_1(v_1,v_2),x_2(v_1,v_2)) = \tilde{f}(v_1,v_2)
  \end{aligned}
  \right.
Indeed, since :math:`(x_1,x_2,x_3)` lies in :math:`\mathcal{TP}`, :math:`x_3`
is not independent of :math:`(x_1,x_2)`:

..  .. math::
..   x_3 = \frac{-x_1(\mathbf{e_1}.\mathbf{V_3})-x_2(\mathbf{e_2}.\mathbf{V_3})}{\mathbf{e_3}.\mathbf{V_3}} */

.. math::
  \left.
  \begin{aligned}
  \tilde{f} \colon \mathcal{TP}\subset\mathbb{R}^2 &\to \mathbb{R}\\
  (v_1,v_2) &\mapsto \tilde{f}(v_1,v_2) = \tilde{f}(v_1(x_1,x_2),v_2(x_1,x_2)) = f(x_1,x_2,x_3)
  \end{aligned}
  \right.

The target is the gradient of :math:`\tilde{f}` in terms of the :math:`\mathcal{TP}` variables
:math:`(v_1,v_2)`. Let us call this gradient :math:`\boldsymbol{\nabla}_\mathcal{TP}`.
It is then possible to apply the :ref:`standard-method` to compute this gradient:


.. math::
   \boldsymbol{\nabla}_\mathcal{TP}W_{kl} = \frac{\partial W}{\partial r}.\boldsymbol{\nabla}_\mathcal{TP}r,
   \quad r = \left\Vert \mathbf{r} \right\Vert = \sqrt{v_{kl,1}^2 + v_{kl,2}^2}
   :label: kernel function gradient TP 1

Which leads to:

.. math::
  \boldsymbol{\nabla}_\mathcal{TP}W_{kl} = -3\frac{10}{\pi r_0^5}\frac{(r_0 - \left\Vert \mathbf{r_{kl}'}\right\Vert)^2}{r_{kl}'}\left\{
  \begin{aligned}
  & v_{kl,1}\mathbf{V_1} + v_{kl,2}\mathbf{V_2}, \quad &0\leq \left\Vert \mathbf{r_{kl}'}\right\Vert \leq  r_0\\
  & 0 , & r_0 <\left\Vert \mathbf{r_{kl}'}\right\Vert
  \end{aligned}
  \right.
  :label: kernel function gradient TP 2

.. math::
  \boldsymbol{\nabla}_\mathcal{TP}\tilde{f_{k}} = -\sum\limits_{l}\tilde{f_{l}}A_{l}\,\boldsymbol{\nabla}W_{kl}
  :label: sph gradient

This gradient can now be expressed in the Cartesian coordinate system.
It is clear that the change of coordinate system was not needed:

.. math::
  \boldsymbol{\nabla}_\mathcal{TP}W_{kl} = -3\frac{10}{\pi r_0^5}\frac{(r_0 - \left\Vert \mathbf{r_{kl}'}\right\Vert)^2}{r_{kl}'}\left\{
  \begin{aligned}
  & r_{kl,1}\mathbf{e_1} + r_{kl,2}\mathbf{e_2} + r_{kl,3}\mathbf{e_3}, \quad &0\leq \left\Vert \mathbf{r_{kl}'}\right\Vert \leq  r_0\\
  & 0 , & r_0 <\left\Vert \mathbf{r_{kl}'}\right\Vert
  \end{aligned}
  \right.

The advantage of computing the gradient in the local coordinate system is if
the components (in flow direction or in cross flow direction) need to be treated
differently.


.. _2_5DSPH:

.. figure:: _static/2_5DSPH.png
        :width: 90%

        Tangent plane and local coordinate system used to apply the SPH method


Flow thickness computation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The particles flow thickness is computed with the help of the grid.
The mass of the particles is interpolated onto the grid using a bilinear interpolation method
(described in Sect.~\ref{sec-particle grid-interpolation}).
Then, dividing the mass at the grid cells by the area of the grid cells, while taking the slope of the
cell into account, returns the flow thickness field on the grid.
This property is interpolated back to the particles which leads to the particle flow thickness property.

We do not compute the flow thickness directly from the particles properties (mass and position) using a
SPH method because it induced instabilities.
Indeed, the cases where too few neighbors are found, lead to very small flow thickness which becomes an
issue for flow thickness dependent
friction laws. Note that using such a SPH method would lead to a fully particular method.
But since the flow thickness is only used in some cases for the friction force computation, using a the
previously describe grid method should not affect significantly the computation.

Friction force discretization
---------------------------------
.. \label{sec-discretizing-friction}

Expressing the friction force term in :Eq:`eq-momentum-particle` for a particle reads:

.. math::
	\mathbf{F}_k^\text{fric} = A_k^b \, \boldsymbol{\tau^b} =
	- {\left\Vert\mathbf{F}_k^\text{fric}\right\Vert}_\text{max}
  \, \mathbf{v}_1
	:label: eq-coulomb-friction-particle

This relation stands if the particle is moving. The starting and stopping processes
satisfy a different equation and are handled differently in the numerical
implementation (using the same equation would lead to a non-physical behavior).
This is described in more details in Sect.~\ref{sec-adding-friction}.


Convergence
------------
.. \label{sec-convergence-criterion}

We are looking for a criterion that relates the properties of the spatial and
temporal discretization to ensure convergence of the numerical solution.
Simply decreasing the time step and increasing the spatial resolution,
by decreasing the grid cell size and kernel radius and increasing the number of
particles, does not ensure convergence.
The analysis from :cite:`MoVi2000` carried out on a very similar problem
(hyperbolic non linear transport equation with a particle and SPH method)
shows that the kernel radius size can not be varied
independently from the time step and number of particles.
Indeed, they show that the numerical solution converges towards the solution of
the equation at the following condition:

.. math::
	\left\{
	\begin{aligned}
			r_{\text{part}} &\to 0\\
			r_{\text{kernel}} &\to 0\\
			\frac{r_{\text{part}}^m}{r_{\text{kernel}}^{m+1}} &\to 0\quad m=2
	\end{aligned}
			\right.
			\quad\mbox{and} \quad dt \leq C r_{\text{kernel}}
	:label: eq-ben-moussa

Where :math:`r_{\text{part}}` represents the "size" of a particle
.. (in our case, the "size" is the
.. basal area of a particle :math:`r_{\text{part}} = \sqrt{A^b/\pi}`
.. \AW{so it is not the basal area but rPart that we use so maybe say: 'in our case,
.. the rPart is based on the basal area of the particle following :math:`r_{\text{part}} = \sqrt{A^b/\pi}`' -
.. however as you bring up the equation for it a bit later in more detail maybe just
.. skip the info in the brackets here?})
, :math:`r_{\text{kernel}}` represents the SPH kernel radius, :math:`dt` is the time
step and :math:`C` a constant.
The conditions in Eq.~\ref{eq-ben-moussa} mean that both :math:`r_{\text{part}}`
(particle size) and :math:`r_{\text{kernel}}` (SPH kernel radius) need to go to zero
but also that the particle size needs to go faster to zero than the SPH kernel radius.
Finally, the time step needs to go to zero and this at the same rate as
:math:`r_{\text{kernel}}`.
The particle size can be expressed as a function of the SPH kernel radius:

.. math::
	r_{\text{part}} = \left(\frac{A^b}{\pi}\right)^{1/2} =
	\left(\frac{A_{\text{kernel}}}{\pi n_{\text{ppk}}}\right)^{1/2}
	=  \frac{r_{\text{kernel}}}{n_{\text{ppk}}^{1/2}},
where the particles basal area was assumed to be a circle.

Note that this does not affect the results except adding a different shape factor in
front of this expression.
:math:`n_{\text{ppk}}` is the number of particles per kernel radius and defines the
density of the particles when initializing a simulation.
Let :math:`n_{\text{ppk}}` be defined by a reference number of particles per kernel
radius :math:`n_{\text{ppk}}^0>0`, a reference kernel radius
:math:`r_{\text{kernel}}^0>0` and real exponent :math:`\alpha`:

.. math::
	n_{\text{ppk}} = n_{\text{ppk}}^0\left(\frac{r_{\text{kernel}}}{r_{\text{kernel}}^0}\right)^{\alpha}

This leads to a :math:`r_{\text{part}}`:

.. math::
	r_{\text{part}} = \left(\frac{{r_{\text{kernel}}^0}^\alpha}{n_{\text{ppk}}^0}\right)^{1/2} r_{\text{kernel}}^{1-\alpha/2}

Replacing :math:`r_{\text{part}}` by the previous equation in
:eq:`eq-ben-moussa` leads to the following condition:

.. math::
	\frac{{r_{\text{kernel}}^0}^\alpha}{n_{\text{ppk}}^0} r_{\text{kernel}}^{-1-\alpha} \to 0
	:label: eq-ben-moussa-new

This brings us to the following choice:

.. math::
	\left\{
	\begin{aligned}
			dt &= C_{\text{time}} r_{\text{kernel}}\\
		 n_{\text{ppk}} &= n_{\text{ppk}}^{0} \left(\frac{r_{\text{kernel}}}{r_{\text{kernel}}^0}\right)^{\alpha}
	\end{aligned}
			\right.
  :label: eq-convergence-relation

Which satisfies the convergence criterion if:

.. math::
 	\alpha < -1
	:label: eq-convergence-criterion

Note that this criterion leaves some freedom on the choice of exponent :math:`\alpha` and that there are no constraints on the reference kernel radius :math:`r_{\text{kernel}}^0` and reference number of particles per kernel radius :math:`n_{\text{ppk}}^0`.
Even though it seems logical to require a minimum number of particles per kernel radius so that enough neighbors are available
to get a reasonable estimate of the gradient.
These parameters should be adjusted according to the expected accuracy of the results and/or the computer power available.
Determining the optimal parameter values for :math:`\alpha`, :math:`r_{\text{kernel}}^0` and :math:`n_{\text{ppk}}^0`, for example according to a user's needs in terms of accuracy and computational efficiency, requires a specific and detailed investigation of the considered case.
In the Sect.~\ref{sec-verification}, we will explore model convergence using the condition eq.~\ref{eq-convergence-criterion} with different values of :math:`\alpha`.


Forces discretization
----------------------

Lateral force
~~~~~~~~~~~~~~

The SPH method is introduced when expressing the flow thickness gradient for each
particle as a weighted sum of its neighbors (:cite:`LiLi2010,Sa2007`).
Which leads to, using the relation :eq:`sph formulation`:

.. math::
    \mathbf{F}_{k}^{\text{lat}} =
    - \rho_0 \, g^\text{eff} \, h A^b \boldsymbol{\nabla}_{\!s} h
    = -m_{k}\,g^\text{eff}\,.\,\frac{1}{\rho_0}\,\sum\limits_{l}{m_{l}}\,\left.\boldsymbol{\nabla}_{\!s}W_{kl}\right\rvert_{l}
    :label: lateral force


Bottom friction force
~~~~~~~~~~~~~~~~~~~~~~~
The bottom friction forces on each particle depend on the chose friction model. Using the SamosAT friction model
(using equation :eq:`sigmab` for the expression of :math:`\sigma^{(b)}_{k}`) the formulation of the bottom friction forec
reads:

.. math::
    \mathbf{F}_{k}^{\text{bot}} = A_{k}\,\boldsymbol{\tau}_{k}^b
    = -A_{k}\,\left(\tau_0 + \tan{\delta}\,\left(1+\frac{R_s^0}{R_s^0+R_s}\right)\,\boldsymbol{\sigma}_{k}^b\cdot\mathbf{n^b}
     + \frac{\rho_0\,\mathbf{\overline{u}}_{k}^2}{\left(\frac{1}{\kappa}\,\ln\frac{h}{R} + B\right)^2}\right)\mathbf{v}_1
    :label: bottom force


Added resistance force
~~~~~~~~~~~~~~~~~~~~~~~
The resistance force on each particle reads (where :math:`h^{\text{eff}}_{k}`
is a function of the average flow thickness :math:`h_{k}`):

.. math::
    \mathbf{F}_{k}^{\text{res}}
    = - \rho_0\,A_{k}\,h^{\text{eff}}_{k}\,C_{\text{res}}\,\overline{\mathbf{u}}_{k}^2\,\mathbf{v}_1
    :label: resistance force

Both the bottom friction and resistance force are friction forces. The expression above represent the maximal
friction force that can be added. This maximal force is added if the particles are flowing. If not, the friction force
equals the driving forces. See :cite:`MaVi2003` for more information.

Entrainment force
~~~~~~~~~~~~~~~~~~~~~~~
The term :math:`- \overline{\mathbf{u}}\,\rho_0\,\frac{\mathrm{d}(A\,h)}{\mathrm{d}t}`
related to the entrained mass in :eq:`momentum-balance3` now reads:

.. math::
    - \overline{\mathbf{u}}_k\,\rho_0\,\frac{\mathrm{d}}{\mathrm{d}t}\,\left(A_{k}\,h_{k}\right)
    = - \overline{\mathbf{u}}_k\,A^{\text{ent}}_{k}\,q^{\text{ent}}_{k}


The mass of entrained snow for each particle depends on the type of entrainment involved
(plowing or erosion) and reads:

.. math::
    \rho_0\,\frac{\mathrm{d}}{\mathrm{d}t}\,\left(A_{k}\,h_{k}\right)
    = \frac{\mathrm{d}\,m_{k}}{\mathrm{d}t}
    = A_{k}^\text{ent}\,q_{k}^{\text{ent}}

with

.. math::
    \begin{aligned}
    A_{k}^{\text{plo}} &= w_f\,h_{k}^{\text{ent}}= \sqrt{\frac{m_{k}}{\rho_0\,h_{k}}}\,h_{k}^{\text{ent}}
    \quad &\mbox{and} \quad &q_{k}^{\text{plo}} = \rho_{\text{ent}}\,\left\Vert \overline{\mathbf{u}}_{k}\right\Vert
    \quad &\mbox{for plowing}\\
    A_{k}^{\text{ero}} &= A_{k} = \frac{m_{k}}{\rho_0\,h_{k}}
    \quad &\mbox{and} \quad &q_{k}^{\text{ero}} = \frac{\tau_{k}^{(b)}}{e_b}\,\left\Vert \overline{\mathbf{u}}_{k}\right\Vert
    \quad &\mbox{for erosion}\end{aligned}

Finaly, the entrainment force reads:

.. math::
    \mathbf{F}_k^{\text{ent}} = -w_f\,(e_s+\,q_{k}^{\text{ent}}\,e_d)\mathbf{v}_1

Adding forces
--------------
The different components are added following an operator splitting method.
This means particle velocities are updated successively with the different forces.

Numerical stability
---------------------

Because the lateral shear force term was removed when deriving the model equations
(because of its relative smallness, :cite:`GrEd2014`), :eq:`eq-momentum-balance-approx`
is hyperbolic.
Hyperbolic systems have the characteristic of carrying discontinuities or shocks which
will cause numerical instabilities.
They would fail to converge if for example an Euler forward in time scheme is used
(:cite:`Le1990`).
Several methods exist to stabilize the numerical integration of an hyperbolic system
of differential equations.
All aim at adding some upwinding in the discretization scheme.
Some methods tackle this problem by introducing some upwinding in the discretization
of the derivatives (:cite:`HaLaLe1983, HaHy1983`).
Others introduce some artificial viscosity (as in :cite:`Mo1992`).

Two options are available to add viscosity to stabilize the numerics. The first option
consists in adding artificial viscosity (``viscOption`` = 1). This is the default
method and is used for operational applications. The second option attempts
to adapt the Lax-Friedrich scheme (usually applied to grids) to the particle method
(``viscOption`` = 2). This method Finally, ``viscOption`` = 0 deactivates any viscosity force.

SAMOS Artificial viscosity
~~~~~~~~~~~~~~~~~~~~~~~~~~

The artificial viscosity force acting on particle :math:`k` then reads:

.. math::
  \begin{aligned}
  \mathbf{F_k^{visc}} = &- \frac{1}{2}\rho_0 C_{Lat} A_k^{\text{Lat}}\Vert\mathbf{d\overline{u}}_k\Vert^2
  \frac{\mathbf{d\overline{u}}_k}{\Vert\mathbf{d\overline{u}}_k\Vert}\\
  = & - \frac{1}{2}\rho_0 C_{Lat} A_k^{\text{Lat}}\Vert\mathbf{d\overline{u}}_k\Vert \mathbf{d\overline{u}}_k,
  \end{aligned}

where the velocity difference reads :math:`\mathbf{d\overline{u}}_k = \overline{\mathbf{u}}_k - \widehat{\overline{\mathbf{u}}}_k` (:math:`\widehat{\overline{\mathbf{u}}}_k` represents the averaged velocity of the neighbor particles and is practically the grid velocity interpolated at the particle position).
:math:`C_{Lat}` is a coefficient that controls the viscous force.
It would be the equivalent of :math:`C_{Drag}` in the case of the drag force.
:math:`C_{Lat}` is a numerical parameter that depends on the grid size.

In this expression, let :math:`\mathbf{\overline{u}}_k^{n}` be the velocity at the beginning of the time step and  :math:`{\overline{\mathbf{u}}_k^{n+1}}^\blacktriangle` be the velocity
after adding the numerical viscosity (\fig{\ref{fig-DFA-solver}}).
In the norm term :math:`\Vert\mathbf{d\overline{u}}_k\Vert` the particle and grid velocity at the beginning of the time step are used.
This ensures no implicit relation on the norm term or on the average velocity :math:`\widehat{\overline{\mathbf{u}}}_k`.
On the contrary, an implicit formulation is used in :math:`\mathbf{d\overline{u}}_k` because the new value of the velocity is used there.
The artificial viscosity force now reads:

.. math::
 	\mathbf{F_k^{visc}} =  -\frac{1}{2}\rho_0 C_{Lat} A_k^{\text{Lat}}\Vert\overline{\mathbf{u}}_k^{n} - \widehat{\overline{\mathbf{u}}}_k^{n}\Vert
 	\left(\left.\overline{\mathbf{u}}_k^{n+1}\right.^\blacktriangle - \widehat{\overline{\mathbf{u}}}_k^{n}\right)

Updating the velocity then gives:

.. math::
 	\left.\overline{\mathbf{u}}_k^{n+1}\right.^\blacktriangle = \frac{\overline{\mathbf{u}}_k^{n} - C_{vis}\widehat{\overline{\mathbf{u}}}_k^{n}}{1 + C_{vis}}

with

.. math::
	C_{vis} = \frac{1}{2}\rho_0 C_{Lat}A_k^{\text{Lat}} \Vert\overline{\mathbf{u}}_k^{n} - \widehat{\overline{\mathbf{u}}}_k^{n}\Vert\frac{\Delta t}{m}.

This approach to stabilize the momentum equation (:eq:`eq-momentum-particle`) is not
optimal for different reasons.
Firstly, it introduces a new coefficient :math:`C_{vis}` which is not a physical
quantity and will require to be calibrated.

Secondly, it artificially adds a force that should be described physically.
So it would be more interesting to take the physical force into account in the first
place.

Potential solutions could be  taking the physical shear force into account,
using for example the :math:`\mu`-I rheology (:cite:`GrEd2014, BaBaGr2016`).
Another option would be to replace the artificial viscosity with a purely
numerical artifact aiming to stabilize the equations such as a SPH version of
the Lax-Friedrich scheme as presented in :cite:`AtSo2005`.

Ata Artificial viscosity: an upwind method based on Lax-Friedrichs scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Shallow Water Equations are well known for being hyperbolic transport equations.
They have the particularity of carrying discontinuities or shocks which will cause
numerical instabilities.

A decentering in time allows to better capture the discontinuities.
This can be done in the manner of the Lax-Friedrich scheme as described in :cite:`AtSo2005`,
which is formally the same as adding a viscous force. Implementing it for the SPH method,
this viscous force applied on a given particle :math:`k` can be expressed as follows:

.. math::
  \mathbf{F}_k^\text{viscosity} = \sum_{l} \frac{m_l}{\rho_l} \Pi_{kl} \boldsymbol{\nabla}W_{kl}

with :math:`\Pi_{kl} = \lambda_{kl}(\mathbf{u}_l - \mathbf{u}_k) \cdot
\frac{\mathbf{r}_{kl}}{\vert\vert \mathbf{r}_{kl} \vert\vert}`, and
:math:`\boldsymbol{\nabla}W_{kl}` is the gradient of the kernel function and
is described in :ref:`DFAnumerics:SPH gradient`.

:math:`\mathbf{u}_{kl} = \mathbf{u}_k - \mathbf{u}_l` is the relative velocity
between particle k and l, :math:`\mathbf{r}_{kl} = \mathbf{x}_k - \mathbf{x}_l` is
the vector going from particles :math:`l` to particle :math:`k` and
:math:`\lambda_{kl} = \frac{c_k+c_l}{2}` with :math:`c_k = \sqrt{gh_l}`
the wave speed. The :math:`\lambda_{kl}` is obtained by turning expressions
related to time and spatial discretization parameters into an expression
on maximal speed between both particles in the Lax Friedrich scheme.

Due to the expression of the viscosity force, it makes sense to
compute it at the same place where the SPH pressure force are computed (for this reason, the
``viscOption`` = 2 corresponding to the "Ata" viscosity option is only available
in combination with the ``sphOption`` = 2).


Adding artificial viscosity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If the viscosity option (``viscOption``) is set to 1, artificial viscosity is added first, as described
in :ref:`DFAnumerics:SAMOS Artificial viscosity` (this is the default option). With ``viscOption`` set to 0, no viscosity is added. Finally, if
``viscOption`` is set to 2, artificial viscosity is added during SPH force computation
(so in :ref:`DFAnumerics:Account for driving forces` according to the :math:`\mathbf{F}_k^\text{viscosity}`
computed in :ref:`DFAnumerics:Ata Artificial viscosity: an upwind method based on Lax-Friedrichs scheme`)


Curvature acceleration term
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. \label{sec-curvature-acc-term-estimation}

The last term of the particular momentum equation (:eq:`eq-momentum-particle`)
as well as the effective gravity :math:`g^{\text{eff}}` are the final terms to be
discretized before the numerical integration.
In both of these terms, the remaining unknown is the curvature acceleration term
:math:`\overline{\mathbf{u}}_k \cdot \frac{\mathrm{d}\mathbf{v}_{3,k}}{\mathrm{d}t}`.
Using the forward Euler time discretization for the temporal derivative of the
normal vector :math:`\mathbf{v}_{3,k}` gives:

.. math::
	\left.\frac{\mathrm{d}\mathbf{v}_{3,k}}{\mathrm{d}t}\right|^n \approx
	\frac{\mathbf{v}_{3,k}^{n+1} - \mathbf{v}_{3,k}^n}{\Delta t}

:math:`\mathbf{v}_{3,k}^n` is a known quantity, the normal vector of the bottom surface at
:math:`\mathbf{x}_k^n` wich is interpolated from the grid normal vector values at the
position of the particle :math:`k` at time :math:`t^n`.
:math:`\mathbf{v}_{3,k}^{n+1}` is unknown since :math:`\mathbf{x}_k^{n+1}` is not known yet,
hence we estimate :math:`\mathbf{x}_k^{n+1}` based the position  :math:`\mathbf{x}_k^n` and
the velocity at :math:`t^n`:

.. math::
	\mathbf{x}_k^{n+1} =\mathbf{x}_k^n + \Delta t \left.\overline{\mathbf{u}}_k^{n+1}\right.^\blacktriangle

This position at :math:`t^{n+1}` is projected onto the topography and
:math:`\mathbf{v}_{3,k}^{n+1}` can be interpolated from the grid normal vector values.

Note that the curvature acceleration term is needed to compute the bottom pressure
(:eq:`eq-pressure-distribution`),  which is used for the bottom friction
computation and for the pressure gradient computation.
The curvature acceleration term can lead to a negative value, which means detachment
of the particles from the bottom surface.
In **com1DFA**, surface detachment is not allowed and if pressure becomes
negative, it is set back to zero forcing the material to remain in contact with the
topography.



Account for entrainment
~~~~~~~~~~~~~~~~~~~~~~~~~~

Entrainment is taken into account by:
* First adding the component representing the loss of momentum due to
	acceleration of the entrained mass :math:`- \overline{\mathbf{u}}_{k}\,A^{\text{ent}}_{k}\,q^{\text{ent}}_{k}`.
	The entrained mass by a particle :math:`k` during a time step :math:`\Delta t` reads:

	.. math::
	    \mathrm{d}m_k^{n}  = m_k^{n+1} - m_k^{n} = \Delta t \,A^{\text{ent}}_{k}\,q^{\text{ent}}_{k}

	Which leads by the way to the new mass of particle :math:`m_k^{n+1}`:

	.. math::
	    m_k^{n+1} =  m_k^{n} + \mathrm{d}m_k^{n} = m_k^{n} + \Delta t \,A^{\text{ent}}_{k}\,q^{\text{ent}}_{k}

	Implicitly updating the velocity leads to (if we call :math:`\overline{\mathbf{u}}_k^{n+1}`
	the velocity before adding the momentum loss and
	:math:`\left.\overline{\mathbf{u}}_k^{n+1}\right.^\blacktriangle` the velocity after):

	.. math::
		\left.\overline{\mathbf{u}}_k^{n+1}\right.^\blacktriangle = \overline{\mathbf{u}}_k^{n+1}
		\frac{m_k^{n}}{m_k^{n} + \mathrm{d}m_k^n} = \overline{\mathbf{u}}_k^{n+1}
		\frac{m_k^{n}}{m_k^{n+1}}

* Second by adding the force due to the need to break and compact the
	entrained mass (:math:`\mathbf{F}_k^{\text{ent}}`) as described in :ref:`DFAnumerics:Entrainment force`.

	.. warning::
		ToDO

	.. math::
	    \mathbf{F}_k^{\text{ent}} = -w_f\,(e_s+\,q_{k}^{\text{ent}}\,e_d)\mathbf{v}_1

Account for driving forces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adding the driving forces -gravity force and lateral forces- is done after adding
the artificial viscosity as described on \fig{\ref{fig-DFA-solver}}.
The velocity is updated as follows
(:math:`{\overline{\mathbf{u}}_k^{n+1}}^\bigstar` is the velocity after taking the
driving force into account):

.. math::
	\begin{aligned}
  	{\overline{\mathbf{u}}_k^{n+1}}^\bigstar &= \left.\overline{\mathbf{u}}_k^{n+1}\right.^\blacktriangle
		+ \frac{\Delta t}{m_k}\mathbf{F}_{k}^{\text{drive}}\\
		&= \left.\overline{\mathbf{u}}_k^{n+1}\right.^\blacktriangle
		+ \frac{\Delta t}{m_k} \left(- m_k \, g^\text{eff}_k \, \boldsymbol{\nabla}_{s} h
			+ m_k \mathbf{g}_s  - m_k \left( \left.\overline{\mathbf{u}}_k^{n+1}\right.^\blacktriangle \cdot \left . \frac{\mathrm{d}\mathbf{v}_{3,k}}{\mathrm{d}t}\right|^n \right)\mathbf{v}_{3,k}^n\right)
	\end{aligned}
	:label: eq-adding-driving-force

Account for friction forces
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Both the bottom friction and resistance forces act against the flow. Two methods are available to add these
forces in com1DFA.

An implicit method:
""""""""""""""""""""

If the velocity of the particle :math:`k` reads :math:`{\overline{\mathbf{u}}_k^{n+1}}^\bigstar`
after adding the driving forces, adding the friction force leads to:

.. math::
  \overline{\mathbf{u}}_k^{n+1} = \frac{{\overline{\mathbf{u}}_k^{n+1}}^\bigstar}{1 + \frac{C_{k}^{\text{fric}}\Delta t}{m_k}}

where :math:`\mathbf{F}_k^{\text{fric}} = -C_{k}^{\text{fric}}{\overline{\mathbf{u}}_k^{n+1}}^\bigstar = \mathbf{F}_k^{\text{res}} + \mathbf{F}_k^{\text{bot}}`
(the two forces are described in :ref:`DFAnumerics:Bottom friction force` and :ref:`DFAnumerics:Added resistance force`).

This implicit method has a few draw-backs. First the flow does not start properly if the
friction angle :math:`\delta` is too close to the slope angle. Second, the flow never properly stops, even if the
particles physically should, i.e. particles keep oscillating back and force around their end position.


An explicit method:
""""""""""""""""""""

The method based on :cite:`MaVi2003` addresses these two issues.
The idea is that the friction force acts against motion, hence it only affects the magnitude of the velocity
and can not be a driving force (:cite:`MaVi2003`).
Moreover, the friction force magnitude depends on the particle state, i.e. if it is
flowing or at rest.
If the velocity of the particle :math:`k` reads :math:`{\overline{\mathbf{u}}_k^{n+1}}^\bigstar`
after adding the driving forces, adding the friction force leads, depending on the
sign of :math:`\frac{m_k \left\Vert{\overline{\mathbf{u}}_k^{n+1}}^\bigstar\right\Vert}{\Delta t} - \left\Vert\mathbf{F}_{k}^{\text{fric}}\right\Vert_{max}`
(where :math:`\left\Vert\mathbf{F}_{k}^{\text{fric}}\right\Vert_{max}`
depends on the chosen friction law introduced in Sect.~\ref{sec-discretizing-friction}), to:

* :math:`\left\Vert\mathbf{F}^{\text{fric}}\right\Vert = \left\Vert\mathbf{F}_{k}^{\text{fric}}\right\Vert_\text{max}` and
  :math:`\overline{\mathbf{u}}_k^{n+1} = {\overline{\mathbf{u}}_k^{n+1}}^\bigstar \left(1 - \frac{\Delta t}{m_k} \frac{\left\Vert\mathbf{F}_{k}^{\text{fric}}\right\Vert_\text{max}}{\left\Vert{\overline{\mathbf{u}}_k^{n+1}}^\bigstar\right\Vert}\right)`,
  if :math:`\frac{m_k \left\Vert{\overline{\mathbf{u}}_k^{n+1}}^\bigstar\right\Vert}{\Delta t} >
  \left\Vert\mathbf{F}_{k}^{\text{fric}}\right\Vert_\text{max}`

* :math:`\left\Vert\mathbf{F}_{k}^{\text{fric}}\right\Vert \leq \left\Vert\mathbf{F}_{k}^{\text{fric}}\right\Vert_\text{max}` and the particle stops moving
  :math:`\overline{\mathbf{u}}_k^{n+1} = 0` before the end of the time step, if
	:math:`\frac{m_k \left\Vert{\overline{\mathbf{u}}_k^{n+1}}^\bigstar\right\Vert}{\Delta t} \leq \left\Vert\mathbf{F}_{k}^{\text{fric}}\right\Vert_\text{max}`

This method prevents the friction force to become a driving force and nonphysically
change the direction of the velocity.
This  would lead to oscillations of the particles instead of stopping.
Adding the friction force following this approach (:cite:`MaVi2003`) allows the
particles to start and stop flowing properly.


Reprojection
~~~~~~~~~~~~~

The last term in :eq:`eq-momentum-particle` (accounting for the curvature effects)
adds a non tangential component allowing
the new velocity to lie in a different plane than the one from the previous time step.
This enables the particles to follow the topography.
But because the curvature term was only based on an estimation
(see Sect.~\ref{sec-curvature-acc-term-estimation}),
it can happen that the new particle position is not necessarily on the topography
and the new velocity does not necessarily lie in the tangent plane at this new
position.
Furthermore, in case of a strong convex curvature and high velocities, the particles
can theoretically be in a free fall state (detachment) as mentioned in
Sect.~\ref{sec-pressure-distribution}.
**com1DFA** does not allow detachment of the particles and the particles are
forced to stay on the topography. This consists in a limitation
of the model/method which will lead to nonphysical behaviors in special cases
(material flowing over a cliff).
In both of the previously mentioned cases, the particles positions are projected back
onto the topography and the velocity direction is corrected to be tangential to the
topography.
The position reprojection is done using an iterative method that attempts to conserve
the distance traveled by each particle between :math:`t^n` and :math:`t^{n+1}``.
The velocity reprojection changes the direction of the velocity but its magnitude is
conserved.

Neighbor search
------------------

The lateral pressure forces are computed via the SPH flow thickness gradient.
This method is based on particle interactions within a certain neighborhood, meaning that it
is necessary to keep track of all the particles within the neighborhood of each particle.
Computing the gradient of the flow thickness at a particle location, requires to
find all the particles in its surrounding. Considering the number of particles and
their density, computing the gradient ends up in computing a lot of
interactions and represents the most computationally expensive part of the dense
flow avalanche simulation. It is therefore important that the neighbor search is fast and efficient.
:cite:`IhOrSoKoTe2014` describe different grid neighbor search
methods. In com1DFA, the simplest method is used. The idea is to locate each
particle in a cell, this way, it is possible to keep track of the particles
in each cell. To find the neighbors of a particle, one only needs to read the
cell in which the particle is located (dark blue cell in :numref:`neighborSearch`),
find the direct adjacent cells in all directions (light blue cells) and
simply read all particles within those cells. This is very easily achieved
on grids because locating a particle in a cell is straightforward and
finding the adjacent cells is also easily done.

.. _neighborSearch:

.. figure:: _static/neighborSearch.png
        :width: 90%

        Support grid for neighbor search:
        if the cell side is bigger than the kernel length :math:`r_{kernel}` (red circle in the picture),
        the neighbors for any particle in any given cell (dark blue square)
        can be found in the direct neighborhood of the cell itself (light blue squares)

.. _partInCell:

.. figure:: _static/partInCell.png
        :width: 90%

        The particles are located in the cells using
        two arrays. indPartInCell of size number of cells + 1
        which keeps track of the number of particles in each cell
        and partInCell of size number of particles + 1 which lists
        the particles contained in the cells.

Particle splitting and merging
-------------------------------
There are two different approaches treating splitting of particles in com1DFA.
The first one only deals with splitting of particles with too much mass('split only'). The second approach,
"split/merge" approach aims at keeping a stable amount of particles within a given range. This is done in order to
guaranty a sufficient accuracy of the sph flow thickness gradient computation.

Split (**default**)
~~~~~~~~~~~~~~~~~~~~
If the ``splitOption`` is set to 0, particles are split because of snow entrainment. In this case,
particles that entrain snow grow, i.e. their mass increases. At one point the mass of the particles is considered to be
too big and this particle is split in two. The splitting operation happens if the mass of the
particle exceeds a threshold value (:math:`mPart > massPerPart \times thresholdMassSplit`), where ``thresholdMassSplit``
is specified in the configuration file and ``massPerPart`` depends on the chosen ``massPerParticleDeterminationMethod``
as defined here: :ref:`com1DFAAlgorithm:Initialize particles`.
When a particle is split a new child particle is created with the same properties as the parent apart from
mass and position. Both parent and child get half of the parent mass. The parent and child's position are
adjusted: the first / second is placed forward / backward in the direction of the velocity
vector at a distance :math:`distSplitPart \times rPart` of the initial parent position. Particles are considered to
have a circular basal surface :math:`A = \frac{m}{\rho} = \pi r^2`.

Split and merge
~~~~~~~~~~~~~~~
If the ``splitOption`` is set to 1 particles are split or merged in order to keep the particle count
as constant as possible within the kernel radius.
Assessing the number of particles within one kernel radius is done based on the particle area. Particles
are assumed to be cylindrical, i.e the base is a circle. For particle ``k`` we have :math:`A_k = \frac{m_k}{\rho}`. The area
of the support domain of the sph kernel function is :math:`\pi r_0^2`. The aim is to keep ``nPPK`` particles within
the kernel radius. The particles are split if the estimated number of particles per kernel radius :math:`\frac{\pi r_0^2}{A_k}`
falls below a given value of :math:`n_{PPK}^{min} = C_{n_{PPK}}^{min}n_{PPK}`. Particles are split using the same
method as in :ref:`DFAnumerics:Only split approach`. Similarly, particles are merged if the estimated
number of particles per kernel radius exceeds a given value :math:`n_{PPK}^{max} = C_{n_{PPK}}^{max}n_{PPK}`.
In this case particles are merged with their closest neighbor. The new position and velocity is the mass
averaged one. The new mass is the sum. Here, two coefficients ``C_{n_{PPK}}^{min}`` and ``C_{n_{PPK}}^{max}`` were
introduced. A good balance needs to be found for the coefficients so that the particles are not constantly split or
merged but also not too seldom. The split and merge steps happen only once per time step and per particle.
