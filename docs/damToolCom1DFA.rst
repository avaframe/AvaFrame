Dam tool com1DFA
=================

The com1DFA computational module provides the option to take the effect of dams into account.
This is done using a ad-hoc method based on particles being reflected/deflected by a dam wall.

Input
-------

The dam is described by the crown line, the slope and the restitution coefficient:


  * crown line as shape file (use the line type and enable the "additional dimensions" option in order
    to specify the z coordinate). The z coordinate corresponds to the absolute height (terrain elevation plus dam height). The dam is then located on the
    left side of the dam. The dam shape files live in the ``avaDir/Inputs/DAM/`` directory (only one file is allowed).

  * the ``slope`` of the dam (in degrees °) between the horizontal plane and the wall to be provided in the shape file
    as an attribute (default value is 60° in the provided examples: avaSlide, avaKot and avaBowl)

  * the restitution coefficient (:math:`\alpha_\text{rest}`), a float between 0 (no reflection in the normal direction) and 1 (full
    reflection) to be specified in the ini file (default value is 0)

Theory
--------

Dam description
~~~~~~~~~~~~~~~~

The dam is described by a crown line, that is to say a series of x, y, z points describing the crown of
the dam (the dam wall is located on the left side of the line), by the slope of the dam wall
(slope measured from the horizontal, :math:`\beta`) and a restitution coefficient (describing if we consider
more elastic or inelastic collisions between the particles and the dam wall, varying between 0 and 1).

The geometrical description of the dam is given on the figure :numref:`fig-DamToolSide`.
The dam crown line (:math:`\mathbf{x_\text{crown}}`) is projected onto the topography, which provides us
with the dam center line (:math:`\mathbf{x_\text{center}}`). We compute the tangent vector to the
center line (:math:`\mathbf{t_f}`). From this tangent vector and the dam slope, it is possible to
compute the wall tangent vector (:math:`\mathbf{t_w}`). Knowing the wall tangent vector and height,
it is possible to determine normal vector to the wall (:math:`\mathbf{n_w}`) and the foot line which
is the intersection between the dam wall and the topography (:math:`\mathbf{x_\text{foot}}`).

When the dam fills up (flow thickness increases), the foot line is modified
(:math:`\mathbf{x_\text{foot}^\text{filled}} = \mathbf{x_\text{foot}} + \frac{h_v}{2} \mathbf{e_z}`).
The normal and tangent vectors to the dam wall are readjusted accordingly.


.. _fig-DamToolSide:

.. figure:: _static/damToolSideView.png
          :width: 90%

          Side view of the dam (cut view). :math:`\mathbf{x_\text{crown}}` describes the crown
          of the dam, :math:`\mathbf{x_\text{center}}` is the vertical projection of the crown
          on the topography (here the light blue line represents the topography).
          The tangent vector to the center line (:math:`\mathbf{t_f}`) is computed from the
          center line points. The tangent vector to the center line with the dam slope angle enable to compute
          the tangent (:math:`\mathbf{t_w}`) and normal (:math:`\mathbf{n_w}`) vector to the dam wall.
          Finally, this normal vector is adjusted depending on the snow thickness at the dam location
          (filling of the dam , :math:`\mathbf{n_w^\text{filled}}`)


In the initialization of the simulation, the dam tangent vector to the center line (:math:`\mathbf{t_f}`),
foot line (:math:`\mathbf{x_\text{foot}}`) and normal vector to the wall (:math:`\mathbf{n_w}`) are computed.
The grid cells crossed by the dam as well as their neighbor cells are memorized (tagged as dam cells).

Dam interaction process
~~~~~~~~~~~~~~~~~~~~~~~~

The particles used in the DFA simulation interact with the dam wall through bouncing (fully elastic up
to nonelastic).

During the computation of the DFA simulation, at each time step, if a particle enters one of the dam cells, the dam
algorithm is called. The first step is to check if the particle crosses the dam during the time step.
If not, the particle position and velocity are updated as if there was no dam. If yes, the intersection point between
the particle trajectory and the dam line is computed. The dam properties are interpolated at the intersection point
(dam tangent and normal vectors). The wall tangent and normal vectors are updated taking the flow thickness
into account. The particle position and velocity are updated taking the dam into account.

Let :math:`\mathbf{x_\text{old}}` and :math:`\mathbf{x_\text{new}}` (respectively :math:`\mathbf{u_\text{old}}`
and :math:`\mathbf{u_\text{new}}`) be the particle position vectors  (respectively velocity vectors) at the beginning
and end of the time step (before the dam interaction). The particle position vector :math:`\mathbf{x_\text{new}^\star}` and
velocity vector :math:`\mathbf{u_\text{new}^\star}` after the dam interaction read:

.. math::

  \mathbf{x^\star_\text{new}} = \mathbf{x_\text{old}} - (1 + \alpha_\text{rest}) \left\{\mathbf{n_w^\text{filled}} \cdot (\mathbf{x_\text{new}} - \mathbf{x_\text{foot}})\right\}\mathbf{n_w^\text{filled}}

.. math::

  \mathbf{u^\star_\text{new}} = \mathbf{u_\text{new}} - (1 + \alpha_\text{rest}) \left\{\mathbf{n_w^\text{filled}} \cdot \mathbf{u_\text{new}}\right\}\mathbf{n_w^\text{filled}}

Finally, the new velocity and position are re-projected onto the topography.
