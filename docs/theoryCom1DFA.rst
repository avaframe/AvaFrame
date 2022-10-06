com1DFA DFA-Kernel theory
============================

.. warning::

   This theory has not been fully reviewed yet. Read its content with a critical mind.

Governing Equations for the Dense Flow Avalanche
------------------------------------------------------

The governing equations of the dense flow avalanche are derived from the
incompressible mass and momentum balance on a Lagrange control volume (:cite:`Zw2000,ZwKlSa2003`).
The detailed derivation can be found in the AvaFrame paper.

Mass balance:
~~~~~~~~~~~~~~~

.. math::
    \frac{\mathrm{d}}{\mathrm{d}t}\underbrace{\int\limits_{V(t)} \rho_0 \,\mathrm{d}V}_{m(t)} = \rho_0 \frac{\mathrm{d}V(t)}{\mathrm{d}t} =
    \oint\limits_{\partial V(t)} q^{\text{ent}} \,\mathrm{d}A,
    :label: mass-balance1

Where :math:`q^{\text{ent}}` represents the snow entrainment rate.

Momentum balance:
~~~~~~~~~~~~~~~~~~~

.. math::
    \frac{\mathrm{d}}{\mathrm{d}t}\int\limits_{V(t)} \rho_0 \, \mathbf{u}\,\mathrm{d}V = \underbrace{\oint\limits_{\partial V(t)}
    \boldsymbol{\sigma}\cdot\mathbf{n}\,\mathrm{d}A}_{\text{surface forces}}
    + \underbrace{\int\limits_{V(t)} \rho_0 \, \mathbf{g}\,\mathrm{d}V}_{\text{body force}} +
    \, \mathbf{F}^\text{ent} + \mathbf{F}^\text{res},
    :label: momentum-balance1


where :math:`\mathbf{u}` is the fluid velocity and :math:`\mathbf{g}` the gravity acceleration.
:math:`\boldsymbol{\sigma} = -pI+\boldsymbol{\mathrm{T}}` represents the
stress tensor, where :math:`I` is the identity tensor, :math:`p` the pressure
and :math:`\boldsymbol{\mathrm{T}}` the deviatoric part of the stress tensor.
:math:`\mathbf{n}` is the normal vector to math:`\partial V(t)`. :math:`\mathbf{F}^{\text{ent}}` represents the force required to break the
entrained snow from the ground and to compress it (since the dense-flow
bulk density is usually larger than the density of the entrained snow,
i.e. :math:`\rho_{\text{ent}}<\rho`) and :math:`\mathbf{F}^{\text{res}}`
represents the resistance force due to obstacles (for example trees).


Boundary conditions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The free surface is defined by :

    .. math:: F_s(\mathbf{x},t) = z-s(x,y,t)=0

The bottom surface is defined by :

    .. math:: F_b(\mathbf{x}) = z-b(x,y)=0

The boundary conditions at the bottom (:math:`\mathcal{S}_b`) and free (:math:`\mathcal{S}_{fs}`) surfaces are introduced.
:math:`\boldsymbol{\sigma^s} = \boldsymbol{\sigma}\cdot\mathbf{n^s}`, respectively :math:`\boldsymbol{\sigma^b} = \boldsymbol{\sigma}\cdot\mathbf{n^b}`, represent the restriction of :math:`\boldsymbol{\sigma}`
to the free surface :math:`\mathcal{S}_{fs}`, respectively the bottom surface :math:`\mathcal{S}_b`:

- traction free free-surface: :math:`\boldsymbol{\sigma_s}\cdot\mathbf{n_s} = \mathbf{0}` on :math:`\mathcal{S}_{fs}`

- impenetrable bottom surface without detachment, :math:`\mathbf{u^b}\cdot\mathbf{n^b} = \mathbf{0}` on :math:`\mathcal{S}_{b}`

- bottom friction law: :math:`\boldsymbol{\tau^b} = \boldsymbol{\sigma^b}\cdot\mathbf{n^b}-
    ((\boldsymbol{\sigma^b}\cdot\mathbf{n^b})\cdot\mathbf{n^b})\mathbf{n^b}
    =\mathbf{f}(\boldsymbol{\sigma^b},\,\overline{\mathbf{u}},\,h,\,\rho_0,\,t,\,\mathbf{x})
    = -f(\boldsymbol{\sigma^b},\,\overline{\mathbf{u}},\,h,\,\rho_0,t,\,\mathbf{x})
    \mathbf{v}_1` on :math:`\mathcal{S}_{b}`

Hypothesis
~~~~~~~~~~~

We consider in the following a shallow flow on moderately curved surfaces. This means
that the aspect ratio, :math:`\qquad\qquad\varepsilon = H / L`, of the characteristic length L
of the avalanche over its characteristic thickness H stays small.

.. _fig-characteristic_size:

.. figure:: _static/characteristic_size.png
        :width: 90%

        Characteristic size of the avalanche along its path (from :cite:`Zw2000`, modified)




Choice of the coordinate system:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The previous equations will be developed in the orthonormal coordinate
system :math:`(B,\mathbf{v_1},\mathbf{v_2},\mathbf{v_3})`, further
referenced as Natural Coordinate System (NCS). In this NCS,
:math:`\mathbf{v_1}` is aligned with the velocity vector at the bottom
and :math:`\mathbf{v_3}` with the normal to the slope, i.e.:

.. math::
   \mathbf{v_1} = \frac{\mathbf{u}}{\left\Vert \mathbf{u}\right\Vert},\quad \mathbf{v_2} = \mathbf{v_3}\wedge\mathbf{v_1},
   \quad \mathbf{v_3} = \mathbf{n^{(b)}}
..   :label: natural-coordinate-system

Thickness averaged equations:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We introduce the following average of a quantity:

.. math::
    \begin{aligned}
    \widetilde{f} &= \frac{1}{V(t)}\int\limits_{V(t)} f\,\mathrm{d}V\\
    \widehat{f} &= \frac{1}{A^b(t)}\int\limits_{\mathcal{S}_b} f\,\mathrm{d}A\\
    \overline{f} &= \frac{1}{h(t)}\int\limits_{0}^{h(t)} f\,\mathrm{d}x_3
    \end{aligned}
    \quad\quad \text{and} \quad \quad
    \begin{aligned}
    \widetilde {f}(x_3) &
    \approx \frac{1}{A^b(\widehat{h}-x_3)}\int\limits_{\mathcal{S}_b}
    \left(\int\limits_{x_3}^{h(t)} f\,\mathrm{d}x_3\right)\,\mathrm{d}A \\
    \overline{f}(x_3) &= \frac{1}{(h-x_3)}\int\limits_{x_3}^{h(t)} f\,\mathrm{d}x_3.
    \end{aligned}

In this NCS and considering a prism-like Control volume, the volume
content :math:`V(t) = A_b(t)h` is obtained by multiplication
of the basal area of the prism, :math:`A_b`.


.. _small-lagrange:

.. figure:: _static/smallLagrange.png
        :width: 90%

        Small Lagrangian prism-like Control volume

Entrainment:
~~~~~~~~~~~~

The snow entrainment is either due to plowing at the front of the avalanche or to erosion
at the bottom. The entrainment rate at the front :math:`q^{\text{plo}}` can be expressed as a function of the
properties of the entrained snow (density :math:`\rho_{\text{ent}}` and
snow thickness :math:`h_{\text{ent}}`), the velocity of the avalanche at the
front :math:`\overline{\mathbf{u}}` and length :math:`w_f` of the front (measured perpendicularly
to the flow velocity :math:`\overline{\mathbf{u}}`). It obviously only happens on the front of
the avalanche:

.. math::
   \oint\limits_{\partial V(t)} q^{\text{plo}}\,\mathrm{d}A = \int\limits_{l_{\text{front}}}\int_b^s q^{\text{plo}}\,
   \mathrm{d}{l}\,\mathrm{d}{z} =  \rho_{\text{ent}}\,w_f\,h_{\text{ent}}\,\left\Vert \overline{\mathbf{u}}\right\Vert
   :label: ploughing

The entrainment rate at the bottom :math:`q^{\text{ero}}` can be expressed as a function of the
bottom area :math:`A_b` of the control volume, the velocity of the avalanche :math:`\overline{\mathbf{u}}`,
the bottom shear stress :math:`\boldsymbol{\tau^b}` and the specific erosion energy :math:`e_b`:

.. math::
    \oint\limits_{\partial V(t)} q^{\text{ero}}\,\mathrm{d}A = \int\limits_{A_b} q^{\text{ero}}\,
    \mathrm{d}A = A_b\,\frac{\boldsymbol{\tau^b}\cdot\mathbf{v}_1}{e_b}\,\left\Vert \overline{\mathbf{u}}\right\Vert
    :label: erosion


This leads in the mass balance :eq:`mass-balance1` to :

.. math::
   \frac{\mathrm{d}V(t)}{\mathrm{d}t} = \frac{\mathrm{d}(A_bh)}{\mathrm{d}t}
   = \frac{\rho_{\text{ent}}}{\rho_0}\,w_f\,h_{\text{ent}}\,\left\Vert \overline{\mathbf{u}}\right\Vert +
   \frac{A_b}{\rho_0}\,\frac{\boldsymbol{\tau^b}\cdot\mathbf{v}_1}{e_b}\,\left\Vert \overline{\mathbf{u}}\right\Vert
   :label: mass-balance2

The force :math:`\mathbf{F}^{\text{ent}}` required to break the entrained snow
from the ground and to compress it is expressed as a function of the required
breaking energy per fracture surface unit :math:`e_s`
(:math:`J.m^{-2}`), the deformation energy per entrained mass element
:math:`e_d` (:math:`J.kg^{-1}`) and the entrained snow thickness
(:cite:`Sa2007,SaFeFr2008,FiFrGaSo2013`):

.. math::
    \mathbf{F}^{\text{ent}} = -w_f\,(e_s+\,q^{\text{ent}}\,e_d)\mathbf{v}_1


Resistance:
~~~~~~~~~~~

The force :math:`\mathbf{F}^{\text{res}}` due to obstacles is expressed
as a function of the characteristic diameter :math:`\overline{d}` and height
:math:`h_{\text{res}}` of the obstacles, the spacing
:math:`s_{\text{res}}` between the obstacles and an empirical
coefficient :math:`c_w` (see :numref:`f-res`). The effective height :math:`h^{\text{eff}}`
is defined as :math:`\min(h, h_{res} )`:

.. math::
   \mathbf{F}^{\text{res}} = -(\frac{1}{2}\,\overline{d}\,c_w/s^2_{\text{res}})\,\rho_0\,A\,
    h^{\text{eff}}\,\overline{u}^2\,\mathbf{v}_1


.. _f-res:

.. figure:: _static/f_res.png
        :width: 90%

        Resistance force due to obstacles (from :cite:`FiKo2013`)



Surface integral forces:
~~~~~~~~~~~~~~~~~~~~~~~~~~

The surface integral is split in three terms, an integral over
:math:`\mathcal{S}_b` the bottom :math:`x_3 = b(x_1,x_2)`, :math:`\mathcal{S}_{fs}` the top
:math:`x_3 = s(x_1,x_2,t)` and :math:`\mathcal{S}_\text{lat}` the lateral surface.
Introducing the boundary conditions :eq:`boundary-conditions` leads to:

.. math::
   \oint\limits_{\partial{V(t)}} \boldsymbol{\sigma}\cdot\mathbf{n}\mathrm{d}A
		= \int\limits_{\mathcal{S}_b} \boldsymbol{\sigma^b}\cdot\mathbf{n^b}\mathrm{d}A
		+ \cancelto{\mathbf{0}}{\int\limits_{\mathcal{S}_{fs}}
		\boldsymbol{\sigma_s}\cdot\mathbf{n_s}\mathrm{d}A}
		+	\int\limits_{\mathcal{S}_{lat}} \boldsymbol{\sigma}\cdot\mathbf{n}\mathrm{d}A
		= \underbrace{\int\limits_{\mathcal{S}_b}
		\boldsymbol{\sigma^b}\cdot\mathbf{n^b}\mathrm{d}A}_{\text{bottom force}}
		+	\underbrace{\oint\limits_{\partial\mathcal{S}_b}
		\left(\int\limits_0^h\boldsymbol{\sigma}\cdot\mathbf{n}\,
		\mathrm{d}x_3\right)\,\mathrm{d}l}_{\text{lateral force}}.
..   :label: surface forces


Mass conservation equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The mass balance :eq:`mass-balance2`
remains unchanged:

.. math::
  \frac{\mathrm{d}V(t)}{\mathrm{d}t} = \frac{\mathrm{d}\left(A_bh\right)}{\mathrm{d}t}
  = \frac{\rho_{\text{ent}}}{\rho_0}\,w_f\,h_{\text{ent}}\,\left\Vert \overline{\mathbf{u}}\right\Vert
  + \frac{A_b}{\rho_0}\,\frac{\boldsymbol{\tau^b}}{e_b}\,\left\Vert \overline{\mathbf{u}}\right\Vert
  :label: mass-balance3


Pressure distribution, thickness integrated pressure and pressure gradient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exploiting the normal component of the momentum equation enables to express the pressure and its gradient:

.. math::
   p^b = p(x_3=0) = \rho_0 h\left(g_{\mathbf{N^b}}
	 + \overline{\mathbf{u}} \cdot \frac{\mathrm{d}\mathbf{N^b}}
   {\mathrm{d}t}\right) = \rho_0 h \, g^\text{eff}

and

.. math::
   h\overline{p} = \int\limits_0^h p(x_3)\,\mathrm{d}x_3
	 = -\rho_0 \frac{h^2}{2}\left(g_{\mathbf{N^b}}
	 + \overline{\mathbf{u}} \cdot \frac{\mathrm{d}\mathbf{N^b}}
   {\mathrm{d}t}\right) = -\rho_0 \frac{h^2}{2} \, g^\text{eff}

and

.. math::
	\label{eq-pressure-gradient}
	 \boldsymbol{\nabla} h\overline{p} = -\rho_0 \, g^\text{eff} \, h \boldsymbol{\nabla} h

Tangential momentum conservation equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the pressure distribution equation and the mass conservation equation, the tangential momentum equations
reads:

.. math::
    \rho_0 V \frac{\mathrm{d}\overline{\mathbf{u}}}{\mathrm{d}t} =
    - A^bp\mathbf{n^b}
    + A^b\boldsymbol{\tau^b}
    - \rho_0 \, g^\text{eff} \, h A^b \boldsymbol{\nabla}_{\!s} h
    + \rho_0 V \mathbf{g}_s
    + \mathbf{F}^{\text{ent}}
    - \mathbf{F}^{\text{res}}
    - \underbrace{\overline{\mathbf{u}}\,\rho_0\,\frac{\mathrm{d}\left(A^b\,h\right)}{\mathrm{d}t}}_{\text{from the mass balance}}
    - \underbrace{\rho_0 V \left( \overline{\mathbf{u}} \cdot \frac{\mathrm{d}\mathbf{v}_3}{\mathrm{d}t} \right)\mathbf{v}_3}_{\text{curvature acceleration}}
    :label: momentum-balance6


The unknown :math:`\overline{\mathbf{u}}` and :math:`h` satisfy :eq:`sigmab`,
:eq:`momentum-balance6` and
:eq:`mass-balance3`. In equation
:eq:`momentum-balance6` the bottom shear
stress :math:`\boldsymbol{\tau^b}` remains unknown and a constitutive equation
has to be introduced in order to completely solve the equations.


Friction Model
~~~~~~~~~~~~~~~~~

The problem can be solved by introducing a constitutive equation which
describes the basal shear stress tensor :math:`\boldsymbol{\tau^b}` as a function
of the flow state of the avalanche.

.. math::
    \boldsymbol{\tau^b} = f(\boldsymbol{\sigma^b},\overline{u},h,\rho_0,t,\mathbf{x})
    :label: samosAT friction model

With

.. math::
   \begin{aligned}
   &\boldsymbol{\sigma^b}\cdot\mathbf{n^b} \qquad &\text{normal component of the stress tensor}\\
   &\overline{u} \qquad &\text{average velocity}\\
   &h \qquad &\text{average flow thickness}\\
   &\rho_0 \qquad &\text{density}\\
   &t \qquad &\text{time}\\
   &\mathbf{x} \qquad &\text{position vector}\end{aligned}

Several friction models already implemented in the simulation tool are
described here.



Mohr-Coulomb friction model
"""""""""""""""""""""""""""""""
The Mohr-Coulomb friction model describes the friction interaction between twos solids.
The bottom shear stress simply reads:

.. math::
 \boldsymbol{\tau^b} = -\left(\tan{\delta}\,\boldsymbol{\sigma^b}\cdot\mathbf{n^b}\right) \mathbf{v}_1

:math:`\tan{\delta}=\mu` is the friction coefficient (and :math:`\delta` the friction angle). The bottom shear stress linearly
increases with the normal stress component :math:`\boldsymbol{\sigma^b}` (:cite:`Zw2000,BaSaGr1999,WaHuPu2004,Sa2007`).

With this friction model, an avalanche starts to flow if the slope inclination is steeper than the
friction angle :math:`\delta`. In the case of an infinite slope of constant inclination,
the avalanche velocity would increase indefinitely. This is unrealistic to model snow
avalanches because it leads to over prediction of the flow velocity.
The Mohr-Coulomb friction model is on the other hand well suited to model
granular flow. Because of its relative simplicity, this friction model is also
very convenient to derive analytic solutions and validate the numerical implementation.

Chezy friction model
""""""""""""""""""""""""
The Chezy friction model describes viscous friction interaction.
The bottom shear stress then reads:

.. math::
 \boldsymbol{\tau^b} = -\left(c_{\text{dyn}}\,\rho_0\,\bar{u}^2\right) \mathbf{v}_1

:math:`c_{\text{dyn}}` is the viscous friction coefficient. The bottom shear stress
is a quadratic function of the velocity. (:cite:`Zw2000,BaSaGr1999,WaHuPu2004,Sa2007`).

This model enables to reach more realistic velocities for avalanche simulations.
The draw back is that the avalanche doesn't stop flowing before the slope inclination approaches zero.
This implies that the avalanche flows to the lowest local point.

Voellmy friction model
"""""""""""""""""""""""""
Anton Voellmy was a Swiss engineer interested in avalanche dynamics :cite:`Vo1955`.
He first had the idea to combine both the Mohr-Coulomb and the Chezy model by summing them up
in order to take advantage of both. This leads to the following friction law:

.. math::
 \boldsymbol{\tau^b} = -\left(\tan{\delta}\,\boldsymbol{\sigma^b}\cdot\mathbf{n^b}
 + c_\text{dyn}\,\rho_0\,\bar{u}^2\right) \mathbf{v}_1


This model is described as Voellmy-Fluid :cite:`Sa2004,Sa2007`, and the turbulent
friction term :math:`\xi` is used instead of :math:`c_{\text{dyn}}`.


SamosAT friction model
""""""""""""""""""""""""

SamosAT friction model is a modification of some more classical models
such as Voellmy model :ref:`Voellmy friction model`. The basal shear stress tensor :math:`\boldsymbol{\tau^b}`
is expressed as (:cite:`Sa2007`):

.. math::
   \boldsymbol{\tau^b} = -\left(\tau_0 + \tan{\delta}\,\left(1+\frac{R_s^0}{R_s^0+R_s}\right)\,\boldsymbol{\sigma^b}\cdot\mathbf{n^b}
    + \frac{\rho_0\,\overline{u}^2}{\left(\frac{1}{\kappa}\,\ln\frac{h}{R} + B\right)^2}\right) \mathbf{v}_1

With

.. math::
   \begin{aligned}
   &\tau_0 \qquad &\text{minimum shear stress}\\
   &R_s \qquad &\text{relation between friction and normal pressure (fluidization factor)}\\
   &R \qquad &\text{empirical  constant}\\
   &R_s^0 \qquad &\text{empirical  constant}\\
   &B \qquad &\text{empirical  constant}\\
   &\kappa \qquad &\text{empirical  constant}\end{aligned}

The minimum shear stress :math:`\tau_0` defines a lower limit below
which no flow takes place with the condition
:math:`\rho_0\,h\,g\,\sin{\alpha} > \tau_0`. :math:`\alpha`
being the slope. :math:`\tau_0` is independent of the flow thickness, which
leeds to a strong avalanche deceleration, especially for avalanches with
low flow heights. :math:`R_s` is expressed as
:math:`R_s = \frac{\rho_0\,\overline{u}^2}{\boldsymbol{\sigma^b}\cdot\mathbf{n^b}}`. Together
with the empirical parameter :math:`R_s^0` the term
:math:`\frac{R_s^0}{R_s^0+R_s}` defines the Coulomb basal friction.
Therefore lower avalanche speeds lead to a higher bed friction, making
avalanche flow stop already at steeper slopes :math:`\alpha`, than
without this effect. This effect is intended to avoid lateral creep of
the avalanche mass (:cite:`SaGr2009`).


Dam 
~~~

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


.. Logarithmic friction model
.. """""""""""""""""""""""""""""
.. Im Gegensatz zum Chezy-Modell, welches eine Näherung der Bodenreibung in turbulenten Strömungen liefert,
.. wird mit der turbulenten Grenzschichttheorie das Verhalten einer solchen Strömung in Bodennähe genauer analysiert.
.. Für die Grundlagen der Grenzschichttheorie.
.. Im Falle einer Lawine kann man von einem vollkommen rauen Untergrund ausgehen. Somit kann in unmittelbarer
.. Bodennäche ein logarithmisches Geschwindigkeitsprofil zugrunde gelegt werden,
.. welches proportional zum Abstand von der Geländeoberfläche $x_3$ ist.
.. Die Verteilung der Geschwindigkeit wird somit über das logarithmisches Geschwindigkeitsprofil $\tilde{u}(x_3)$ ersetzt.
.. Diese Profil wird ''universelles Wandgesetz'' genannt.
..
.. .. math::
..   \begin{aligned}
..     \frac{\tilde{u}}{u_{\tau}} &= \frac{1}{\kappa}\,\ln{\frac{x_3}{R}} + B\\
..     &\text{mit}\\
..     u_{\tau} &= \sqrt{\frac{\boldsymbol{\tau^b}}{\bar{\rho}}},
..   \end{aligned}
..
..
.. .. math::
..   \begin{aligned}
..   &R \qquad &\text{Surface rugosity}\\
..   &B \qquad &\text{empirical  constant}\\
..   &\kappa \qquad &\text{Karman constant (0.4 fur Newtonsche Flussigkeiten)}\\
..   &u_{\tau} \quad &\text{Schubspannungsgeschwindigkeit}\end{aligned}
..
..
..
.. sind. Für Kanalströmungen mit vorgegebener Fließhöhe $\bar{h}$ lässt sich zudem
.. das Geschwindigkeitsprofil auch in größeren Entfernungen von der Wand logarithmisch darstellen.
.. Daraus folgt das ''Mittengesetz''.
..
.. .. math::
..   \begin{align*}
..     \frac{\tilde{u}}{u_{\tau}} &= \frac{\tilde{u}_\text{max}}{u_{\tau}} + \frac{1}{\kappa}\,\ln{\frac{x_3}{\bar{h}}}
..     \intertext{mit}
..     \tilde{u}_\text{max} &= \tilde{u}(x_3 = \bar{h})
..   \end{align*}
..
.. Kombiniert man nun diese beiden Gesetze, indem das Mittengesetz dem Wandgesetz gleichsetzt,
.. erhält man
..
.. .. math::
..     \frac{\tilde{u}_\text{max}}{u_{\tau}} = \frac{1}{\kappa}\,\ln{\frac{\bar{h}}{R}} + B.
..
.. Durch Einsetzen für $u_{\tau}$ und Ersetzen von $\tilde{u}_\text{max}$ durch die in Kapitel \ref{sec:vereinfachtegleichungen}
.. tiefengemittelte Geschwindigkeit $\bar{u}$
.. erhält man nach Umformen schließlich eine Beziehung für die gesuchte Bodenschubspannung $\boldsymbol{\tau^b}$.
..
.. .. math::
..     \boldsymbol{\tau^b} = \frac{\bar{\rho}\,\bar{u}^2}{\left(\frac{1}{\kappa}\,\ln{\frac{\bar{h}}{R}}+B\right)^2}
..
..
.. Dieses Modell lässt sich wie beim Voellmy-Modell mit der Coulomb'schen Reibung kombinieren.
..
.. .. math::
..     \boldsymbol{\tau^b} = \tan{\delta}\,\boldsymbol{\sigma^b} +
..     \frac{\bar{\rho}\,\bar{u}^2}{\left(\frac{1}{\kappa}\,\ln{\frac{\bar{h}}{R}}+B\right)^2}
..
..
.. Im Gegensatz zum klassischen Voellmy-Modell ist die Bodenschubspannung also auch von der Fließmächtigkeit $\bar{h}$
.. und Bodenrauhigkeit $R$ abhängig.
.. Steigende Fließmächtigkeit, bzw. abnehmende Bodenrauhigkeit führt demnach
.. zu einer Reduktion der Reibung \citep[vgl.][]{Sa2007}.



.. .. _fig-infinitesimales_element:
..
.. .. figure:: _static/infinitesimales_element.png
..         :width: 90%
..
..         Infinitesimal volume element and acting forces on it (from [FiKo2013]_)
