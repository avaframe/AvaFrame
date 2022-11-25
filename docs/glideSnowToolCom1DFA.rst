Glide snow tool com1DFA
========================

Glide snow tool com1DFA
----------------------------

The com1DFA computational module provides the option to add an elastic cohesion force between particles. This can
be used for example to simulate small glide snow avalanches.

The cohesion between particles is modeled using elastic bonds with maximum strain rate.
This means that the cohesion force :math:`\mathbf{F}_{lk}` exerted by :math:`l` on :math:`k` is expressed by:

.. math::

  \mathbf{F}_{lk} = -\mathbf{F}_{kl} =  \frac{L-L_0}{L_0} A_{kl} E \mathbf{e_{kl}}

Where :math:`L` and :math:`L_0` are the distances between the particles :math:`k` and :math:`l` at a time t and at rest
(initial time step), :math:`E` is the elastic modulus of the material, :math:`A_{kl}` is the contact area between
particles :math:`k` and :math:`l` and :math:`\mathbf{e_{kl}}` is the unit vector


.. math::

  \epsilon = \frac{L-L_0}{L_0}

represents the strain (:math:`\epsilon` is positive in the case of elongation and negative in compression).
If the strain exceeds a critical value :math:`\epsilon_c` the bond between the particles :math:`k` and :math:`l`
breaks.

The contact area :math:`A_{kl}` between particles :math:`k` and :math:`l` reads:

.. math::

  A_{kl} = h_{kl} d_{kl} = \frac{h_k + h_l}{2} \frac{L}{\sqrt{3}}

:math:`h_{kl}` represents the height of the contact surface and :math:`d_{kl}` represents the
length of the contact surface in the horizontal plane assuming that each particle has 6 neighbors.



Input
-------

The following parameters are required for activating and configuring cohesion:

  * one activates the cohesion by setting ``cohesion`` to 1

  * the maximum strain before breaking of the bond ``cohesionMaxStrain``

  * the Young modulus in N/m² ``cohesiveSurfaceTension``


Initialization of bonds
-------------------------

If the elastic cohesion (``cohesion`` set to 1) is activated, the bonds between particles are initialized.
The construction of the bonds is done by building a triangular mesh on the particles
(used as a point cloud) and the Delaunay `triangulation function <https://matplotlib.org/stable/api/tri_api.html#matplotlib-tri>`_
from matplotlib.
The length :math:`L_0` of the bonds at rest (initial distance between two bonded particles) is also initialized here.


.. Note:: It is recommended to use the triangular (``triangular``) initialization option (:ref:`com1DFAAlgorithm:Initialize particles`)
          rather than the random or the uniform one in the case where cohesion is activated.
