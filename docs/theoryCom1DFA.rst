com1DFA DFA-Kernel theory
============================

Governing Equations for the Dense Flow Avalanche
------------------------------------------------------

The derivation is based on the thickness integration of the three-dimensional Navier-Stokes equations,
using a Lagrangian approach with a terrain following coordinate system (:cite:`Zw2000,ZwKlSa2003`).
The equations are simplified using the assumption of shallow flow on a mildly curved topography,
meaning flow thickness is considerably smaller than the length and width of the avalanche and
it is considerably smaller than the topography curvature radius.

We consider snow as the material, however the choice of material does not influence the derivation in the first place.
We assume constant density :math:`\rho_0` and a flow on a surface :math:`\mathcal{S}`, subjected to the gravity force
and friction on the surface :math:`\mathcal{S}`.
If needed, additional processes such as entrainment or other external effects can be taken into account.
The mass conservation equation applied to a Lagrangian volume of material :math:`V(t)` reads:

Mass balance
~~~~~~~~~~~~~~~

.. math::
    \frac{\mathrm{d}}{\mathrm{d}t}\underbrace{\int\limits_{V(t)} \rho_0 \,\mathrm{d}V}_{m(t)} = \rho_0 \frac{\mathrm{d}V(t)}{\mathrm{d}t} =
    \oint\limits_{\partial V(t)} q^{\text{ent}} \,\mathrm{d}A,
    :label: eq-mass-balance1

where the source term :math:`q^{\text{ent}}` represents the snow entrainment rate (incoming mass flux).

Momentum balance
~~~~~~~~~~~~~~~~~~~

.. math::
    \frac{\mathrm{d}}{\mathrm{d}t}\int\limits_{V(t)} \rho_0 \, \mathbf{u}\,\mathrm{d}V = \underbrace{\oint\limits_{\partial V(t)}
    \boldsymbol{\sigma}\cdot\mathbf{n}\,\mathrm{d}A}_{\text{surface forces}}
    + \underbrace{\int\limits_{V(t)} \rho_0 \, \mathbf{g}\,\mathrm{d}V}_{\text{body force}} +
    \, \mathbf{F}^\text{ent} + \mathbf{F}^\text{res},
    :label: eq-momentum-balance1


where :math:`\mathbf{u}` is the fluid velocity and :math:`\mathbf{g}` the gravity acceleration.
:math:`\boldsymbol{\sigma} = -pI+\boldsymbol{\mathrm{T}}` represents the
stress tensor, where :math:`I` is the identity tensor, :math:`p` the pressure
and :math:`\boldsymbol{\mathrm{T}}` the deviatoric part of the stress tensor.
:math:`\mathbf{n}` is the normal vector to math:`\partial V(t)`. :math:`\mathbf{F}^{\text{ent}}` represents the force required to break the
entrained snow from the ground and to compress it (since the dense-flow
bulk density is usually larger than the density of the entrained snow,
i.e. :math:`\rho_{\text{ent}}<\rho`) and :math:`\mathbf{F}^{\text{res}}`
represents the resistance force due to obstacles (for example trees).

Hypothesis
~~~~~~~~~~~

* We consider in the following a shallow flow on moderately curved surfaces. This means
  that the aspect ratio, :math:`\varepsilon = H / L`, of the characteristic length L
  of the avalanche over its characteristic thickness H stays small.

.. _fig-characteristic_size:

.. figure:: _static/characteristic_size.png
        :width: 90%

        Characteristic size of the avalanche along its path (from :cite:`Zw2000`, modified)

* A control volume :math:`V(t)` is assumed to be a small prism shape extending from the bottom surface :math:`\mathcal{S}_b` (lying on the topography
  :math:`\mathcal{S}`) up to the free surface in the surface normal direction :math:`\mathbf{N^b}` as
  illustrated in :numref:`small-lagrange`.
  Note that the bottom surface :math:`\mathcal{S}_b` of area :math:`A^b` has no predefined shape.
  The octagonal shape used in :numref:`small-lagrange` is just one possible option.


.. _small-lagrange:

.. figure:: _static/smallLagrange.png
          :width: 90%

          Small Lagrangian prism-like Control volume

Choice of the coordinate system and thickness averaged quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The normal vector :math:`\mathbf{N^b}` to the bottom surface is pointing upwards
whereas  :math:`\mathbf{n^b}=-\mathbf{N^b}` is the bottom normal vector to the Lagrangian
control volume (pointing out of the volume).

.. math::
   V(t) = \int\limits_{V(t)}\,\mathrm{d}V = \int\limits_{\mathcal{S}_b}
	 \left(\int\limits_b^s\det(\mathbf{J})\,\mathrm{d}x_3\right)\,\mathrm{d}A
  :label: eq-lagrange-volume

where :math:`\mathbf{J}` is the transformation matrix from the Cartesian
coordinate system to the Natural coordinate system (NCS).
The NCS is an orthonormal
coordinate system :math:`(\mathbf{v_1},\mathbf{v_2},\mathbf{v_3})` aligned
with the bottom surface. :math:`\mathbf{v_3}=\mathbf{N^b}=-\mathbf{n^b}` is the normal
vector to the bottom surface pointing upwards.
:math:`\mathbf{v_1}` is pointing in the direction of the thickness integrated fluid velocity
:math:`\overline {\mathbf{u}}` (introduced below).

.. math::
 	 \mathbf{v_1} = \frac{\overline{\mathbf{u}}}
 	 {\left\Vert\overline{\mathbf{u}}\right\Vert},\quad \mathbf{v_2} =
 	 \frac{\mathbf{v_3}\wedge\mathbf{v_1}}{\left\Vert
 	 \mathbf{v_3}\wedge\mathbf{v_1}\right\Vert},
 	 \quad \mathbf{v_3} = \mathbf{N^b}
  :label: eq-natural-coordinate-system

In the case of shallow flow on moderately curved surfaces, :math:`\det(\mathbf{J}) = (1 -
\kappa_1 x_3)(1 - \kappa_2 x_3) \approx1`. :math:`\kappa_{\{1,2\}}` represent
the surface curvature in :math:`\mathbf{v}_{\{1,2\}}` directions and :math:`x_3`
is the elevation from the bottom surface in the direction :math:`\mathbf{N^b}`.
This approximation is valid if the curvature radius is much larger then the
flow thickness :math:`h`. In this case, the control volume reads:

.. math::
	 V(t) \approx \int\limits_{\mathcal{S}_b}\!\!
	 \underbrace{\int\limits_b^s\,\mathrm{d}x_3}_{h(t)}\,\mathrm{d}A
  :label: eq-lagrange-volume2


We introduce the following average of a quantity (where :math:`f` is a scalar or vector function):

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

Note that :math:`\widetilde {f}(0)=\widetilde {f}` and :math:`\overline{f}(0)=\overline{f}`.
When the control volume goes to 0, i.e. basal area goes to a point,
:math:`\widetilde {f}\xrightarrow{A^b\xrightarrow{}0}\overline{f}`
and :math:`{\widehat{f}\xrightarrow{A^b\xrightarrow{}0}f}`.

The NCS has some interesting properties that will be useful for projecting and solving the equations.
Because of the orthogonality of this NCS, we have
:math:`\mathbf{v}_i\cdot\mathbf{v}_j = \delta_{ij},\, \{i,j\}\in \{1,2,3\}^2`
which gives after time derivation:

.. math::
  \frac{\mathrm{d}\mathbf{v}_i\cdot\mathbf{v}_j}{\mathrm{d}t} =
  \mathbf{v}_i\cdot\frac{\mathrm{d}\mathbf{v}_j}{\mathrm{d}t} +
  \frac{\mathrm{d}\mathbf{v}_i}{\mathrm{d}t}\cdot\mathbf{v}_j = 0,
  :label: eq-natural-coordinate-system-identity-base

meaning that:

.. math::
  \left\{
  \begin{aligned}
    \frac{\mathrm{d}\mathbf{v}_i}{\mathrm{d}t}\cdot\mathbf{v}_i &= 0
    \implies \frac{\mathrm{d}\mathbf{v}_i}{\mathrm{d}t} \perp \mathbf{v}_i\\
    \frac{\mathrm{d}\mathbf{v}_i}{\mathrm{d}t}\cdot\mathbf{v}_j &=
    -\mathbf{v}_i\cdot\frac{\mathrm{d}\mathbf{v}_j}{\mathrm{d}t},\quad i \neq j.
  \end{aligned}
  \right.
  :label: eq-natural-coordinate-system-identity-1

It is possible to express :math:`\frac{\mathrm{d}\mathbf{v}_1}{\mathrm{d}t}` in terms
of :math:`(\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3)` and using orthogonality
of :math:`\frac{\mathrm{d}\mathbf{v}_i}{\mathrm{d}t}` and :math:`\mathbf{v}_i`:

.. math::
  \frac{\mathrm{d}\mathbf{v}_1}{\mathrm{d}t} = \alpha_i\mathbf{v}_i
  = \cancelto{\mathbf{0}}{\alpha_1\mathbf{v}_1} + \alpha_2\mathbf{v}_2 + \alpha_3\mathbf{v}_3, \quad \alpha_i =
  \frac{\mathrm{d}\mathbf{v}_1}{\mathrm{d}t}\cdot\mathbf{v}_i
  :label: eq-v1-in-ncs

The derivative of the thickness integrated velocity decomposes to:

.. math::
  \frac{\mathrm{d}\overline{\mathbf{u}}}{\mathrm{d}t} =
  \frac{\mathrm{d}\overline{u}_1\mathbf{v}_1}{\mathrm{d}t} =
  \overline{u}_1\frac{\mathrm{d}\mathbf{v}_1}{\mathrm{d}t} +
  \frac{\mathrm{d}\overline{u}_1}{\mathrm{d}t}\mathbf{v}_1 =
  \overline{u}_1(\alpha_2\mathbf{v}_2 + \alpha_3\mathbf{v}_3) +
  \frac{\mathrm{d}\overline{u}_1}{\mathrm{d}t}\mathbf{v}_1
  :label: eq-du-in-ncs


Boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To complete the conservation :eq:`eq-mass-balance1` and :eq:`eq-momentum-balance1` the following boundary conditions at the bottom (:math:`\mathcal{S}_b`) and free (:math:`\mathcal{S}_{fs}`) surfaces are introduced.
:math:`\boldsymbol{\sigma^s} = \boldsymbol{\sigma}\cdot\mathbf{n^s}`, respectively :math:`\boldsymbol{\sigma^b} = \boldsymbol{\sigma}\cdot\mathbf{n^b}`, represent the restriction of :math:`\boldsymbol{\sigma}`
to the free surface :math:`\mathcal{S}_{fs}`, respectively the bottom surface :math:`\mathcal{S}_b`:

- traction free free-surface: :math:`\boldsymbol{\sigma_s}\cdot\mathbf{n_s} = \mathbf{0}` on :math:`\mathcal{S}_{fs}`

- impenetrable bottom surface without detachment, :math:`\mathbf{u^b}\cdot\mathbf{n^b} = \mathbf{0}` on :math:`\mathcal{S}_{b}`

- bottom friction law: :math:`\boldsymbol{\tau^b} = \boldsymbol{\sigma^b}\cdot\mathbf{n^b}-((\boldsymbol{\sigma^b}\cdot\mathbf{n^b})\cdot\mathbf{n^b})\mathbf{n^b}=\mathbf{f}(\boldsymbol{\sigma^b},\,\overline{\mathbf{u}},\,h,\,\rho_0,\,t,\,\mathbf{x}) = -f(\boldsymbol{\sigma^b},\,\overline{\mathbf{u}},\,h,\,\rho_0,t,\,\mathbf{x})\mathbf{v}_1` on :math:`\mathcal{S}_{b}`


Entrainment
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


Resistance
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


Constitutive relation: friction force
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Up to now in the derivations, the bottom shear stress :math:`\boldsymbol{\tau^b}` is unknown.
To close the momentum equation, a constitutive equation describing the basal shear stress tensor
:math:`\boldsymbol{\tau^b}` as a function of the avalanche flow state is required:

.. math::
  \boldsymbol{\tau^b} =
  \mathbf{f}(\boldsymbol{\sigma^b},\overline{\mathbf{u}},h,\rho_0,t,\mathbf{x})
  :label: eq-bottom-frict-law

..  % where $\boldsymbol{\sigma^b}$ represents the normal component of the stress tensor at the bottom,
.. % $\overline{\mathbf{u}}$ the thickness average velocity, $h$ the flow thickness $\rho_0$ the density of the material,
.. % $t$ and $\mathbf{x}$ the time and position vector.

.. In the following, we employ a Mohr-Coulomb friction model that describes the friction interaction between two solids.
..
.. The bottom shear stress reads:
..
.. .. math::
..   \boldsymbol{\tau^b} =
..   -\tan{\delta}\,\boldsymbol{\sigma^b}\cdot\mathbf{n^b}\, \frac{\mathbf{\overline{u}}}{\Vert\mathbf{\overline{u}}\Vert},
..
..
.. where :math:`\delta` is the friction angle and :math:`\mu=\tan{\delta}` is the friction coefficient.
.. The bottom shear stress linearly increases with the normal stress component :math:`p^b` \citep{BaSaGr1999}.

Different friction models accounting for the influence of flow velocity, flow thickness, etc. have been proposed
(e.g. the Voellmy model :cite:`Vo1955`).
Changing the friction model means changing the :math:`\mathbf{f}` function (:eq:`eq-bottom-frict-law`).
In the **com1DFA** module, three friction models are available.
First a Coulomb one which is used in this paper.
Second a Voellmy friction model (:cite:`Vo1955`) and third the samosAT friction model which is the one used for
hazard mapping by Austrian federal agencies (:cite:`Sa2007)`.

With Mohr-Coulomb friction an avalanche starts to flow if the slope inclination exceeds the friction angle
:math:`\delta`.
In the case of an infinite slope of constant inclination, the avalanche velocity would increase indefinitely.
Using this friction law, flow velocity is overestimated and hence is not suited to model the flow of snow avalanches.
However, because of its relative simplicity, the Mohr-Coulomb friction model is convenient for deriving
analytical solutions and testing numerical implementations.


Expression of surface forces in the NCS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Taking advantage of the NCS and using the boundary conditions, it is possible to split the surface forces into  bottom
(on :math:`\mathcal{S}_b`), lateral (on :math:`\mathcal{S}_\text{lat}`) and free surface forces (on :math:`\mathcal{S}_{fs}`)
and perform further simplifications:

.. math::
		\begin{aligned}
    \oint\limits_{\partial{V(t)}} \boldsymbol{\sigma}\cdot\mathbf{n}\,\mathrm{d}A
		&= \int\limits_{\mathcal{S}_b} \boldsymbol{\sigma^b}\cdot\mathbf{n^b}\,\mathrm{d}A
		+ \cancelto{\mathbf{0}}{\int\limits_{\mathcal{S}_{fs}}
		\boldsymbol{\sigma_s}\cdot\mathbf{n_s}\,\mathrm{d}A}
		+	\int\limits_{\mathcal{S}_{lat}} \boldsymbol{\sigma}\cdot\mathbf{n}\,\mathrm{d}A\\
		&= \underbrace{\int\limits_{\mathcal{S}_b}
		\boldsymbol{\sigma^b}\cdot\mathbf{n^b}\,\mathrm{d}A}_{\text{bottom force}}
		+	\underbrace{\oint\limits_{\partial\mathcal{S}_b}
		\left(\int\limits_0^h\boldsymbol{\sigma}\cdot\mathbf{n}\,
		\mathrm{d}x_3\right)\,\mathrm{d}l}_{\text{lateral force}}
    \end{aligned}.
  :label: surface forces

Using the notations introduced in Sect.~\ref{sec-averaged-quantities} and the decomposition of the stress tensor,
the bottom force can be expressed as a surface normal component and a surface tangential one:

.. math::
		\begin{aligned}
		\int\limits_{\mathcal{S}_b} \boldsymbol{\sigma^b}\cdot\mathbf{n^b}\,\mathrm{d}A
		&= \int\limits_{\mathcal{S}_b} (-p^bI
		+ \boldsymbol{\mathrm{T}})\cdot\mathbf{n^b}\,\mathrm{d}A
		= -\int\limits_{\mathcal{S}_b} p^b\mathbf{n^b}\,\mathrm{d}A +
		\int\limits_{\mathcal{S}_b} \boldsymbol{\mathrm{T}} \cdot \mathbf{n^b}\,\mathrm{d}A\\
		&= -\int\limits_{\mathcal{S}_b} p^b\mathbf{n^b}\,\mathrm{d}A
		+ \int\limits_{\mathcal{S}_b} \boldsymbol{\tau^b}\,\mathrm{d}A
		= -A^b\widehat{p^b\mathbf{n^b}}+A^b\widehat{\boldsymbol{\tau^b}}
    \end{aligned},
  :label: eq-basal-surface-forces

where :math:`\boldsymbol{\tau^b}` is the basal friction term (introduced in Sect.~\ref{sec-boundary-conditions}).
Applying Green's theorem, the lateral force reads:

.. math::
		\begin{aligned}
		\oint\limits_{\partial\mathcal{S}_b} \left(\int\limits_0^h
		\boldsymbol{\sigma}\cdot\mathbf{n}\,\mathrm{d}x_3\right)\,\mathrm{d}l
		&= \oint\limits_{\partial\mathcal{S}_b} \left(\int\limits_0^h
		(-pI	+ \boldsymbol{\mathrm{T}})\,\mathrm{d}x_3\right)\cdot\mathbf{n}
		\,\mathrm{d}l\\
    &= -\oint\limits_{\partial\mathcal{S}_b}
		\left(\int\limits_0^hp\,\mathrm{d}x_3\right)\cdot\mathbf{n}\,\mathrm{d}l
		+ \oint\limits_{\partial\mathcal{S}_b} \left(\int\limits_0^h\boldsymbol{\mathrm{T}}
		\,\mathrm{d}x_3\right)\cdot\mathbf{n}\,\mathrm{d}l\\
		& = -\oint\limits_{\partial\mathcal{S}_b} h\overline{p}\mathbf{n}\,\mathrm{d}l
		+ \oint\limits_{\partial\mathcal{S}_b} h\boldsymbol{\overline{\mathrm{T}}}\cdot\mathbf{n}\,\mathrm{d}l
		= -\int\limits_{\mathcal{S}_b} \boldsymbol{\nabla} h\overline{p}\,\mathrm{d}A
		+ \int\limits_{\mathcal{S}_b} \boldsymbol{\nabla} h\boldsymbol{\overline{\mathrm{T}}}\,\mathrm{d}A\\
		& = -A^b\widehat{\boldsymbol{\nabla} h\overline{p}} + A^b \widehat{\boldsymbol{\nabla} h\boldsymbol{\overline{\mathrm{T}}}}
		\end{aligned}
  :label: eq-lateral-surface-forces

Equations :eq:`eq-basal-surface-forces` and :eq:`eq-lateral-surface-forces` represent the
thickness integrated form of the surface forces and can now be used to write the thickness
integrated momentum equation.

Thickness integrated mass conservation equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The mass balance :eq:`mass-balance2`
remains unchanged:

.. math::
  \frac{\mathrm{d}V(t)}{\mathrm{d}t} = \frac{\mathrm{d}\left(A_bh\right)}{\mathrm{d}t}
  = \frac{\rho_{\text{ent}}}{\rho_0}\,w_f\,h_{\text{ent}}\,\left\Vert \overline{\mathbf{u}}\right\Vert
  + \frac{A_b}{\rho_0}\,\frac{\boldsymbol{\tau^b}}{e_b}\,\left\Vert \overline{\mathbf{u}}\right\Vert
  :label: mass-balance3



Thickness integrated momentum equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the definitions of average values given in Sect.~\ref{sec-averaged-quantities}
and the decomposition of the surface forces given by :eq:`eq-basal-surface-forces`
and :eq:`eq-lateral-surface-forces` combined with the expression of the
entrainment force detailed in Appendix~\ref{ap-entrainment-force}, the momentum equation reads:

.. math::
  \rho_0 \frac{\mathrm{d}V(t) \widetilde {\mathbf{u}}}{\mathrm{d}t} = \rho_0 V
  \frac{\mathrm{d}\widetilde {\mathbf{u}}}{\mathrm{d}t} +
  \rho_0 \widetilde {\mathbf{u}} \frac{\mathrm{d}V}{\mathrm{d}t} = \oint\limits_{\partial V(t)}
  \boldsymbol{\sigma}\cdot\mathbf{n}\,\mathrm{d}A + \rho_0 V \mathbf{g} +
  \mathbf{F}^{\text{ext}}
 :label: eq-momentum-balance2


which leads to:

.. math::
  \rho_0 V \frac{\mathrm{d}\widetilde {\mathbf{u}}}{\mathrm{d}t} =
  \underbrace{-A^b\widehat{p\mathbf{n^b}}}
  _{\substack{\text{bottom} \\ \text{ normal force }}}
  \underbrace{+A^b\widehat{\boldsymbol{\tau^b}}}
  _{\substack{\text{bottom} \\ \text{ shear force }}}
  \underbrace{-A^b\widehat{\boldsymbol{\nabla} h\overline{p}}}
  _{\substack{\text{lateral} \\ \text{ pressure force }}}
  \underbrace{+A^b\cancelto{O(\boldsymbol{\epsilon}^2)}
  {\widehat{\boldsymbol{\nabla} h\overline{\boldsymbol{T}}}}}
  _{\substack{\text{lateral} \\ \text{ shear force }}}
  + \rho_0 V \mathbf{g} + \mathbf{F}^{\text{ext}} \underbrace{-\widetilde {\mathbf{u}}\,\oint\limits_{\partial V(t)}
  q^{\text{ent}}\,\mathrm{d}A}
  _{\substack{\text{ momentum loss } \\ \text{ entrainment }}}
  :label: eq-momentum-balance3

The lateral shear stress term is neglected because of its relative smallness
in comparison to the other terms as shown by the dimensional analysis carried out in
:cite:`GrEd2014`. The mass conservation reads:

.. math::
  \rho_0 \frac{\mathrm{d}V}{\mathrm{d}t} =
  \oint\limits_{\partial V(t)} q^{\text{ent}} \,\mathrm{d}A,
 :label: eq-mass-balance2

Using the approximations from Sect.~\ref{sec-averaged-quantities}, the momentum
equation becomes:

.. math::
  \rho_0 V \frac{\mathrm{d}\overline{\mathbf{u}}}{\mathrm{d}t} = - A^bp\mathbf{n^b}
  + A^b\boldsymbol{\tau^b} - A^b\boldsymbol{\nabla} h\overline{p} + \rho_0 V \mathbf{g}
  + \mathbf{F}^{\text{ext}}
  - \overline{\mathbf{u}}\oint\limits_{\partial V(t)} q^{\text{ent}} \,\mathrm{d}A,
  :label: eq-momentum-balance-approx


where all quantities are evaluated at the center of the basal area (point O in :numref:`small-lagrange`).
This equation is projected in the normal direction :math:`\mathbf{v}_3 = \mathbf{N^b}`
to get the expression of the basal pressure :math:`p^b`. The projection of this same
equation on the tangential plane leads to the differential equations satisfied by
:math:`\overline{\mathbf{u}}`.

Pressure distribution, thickness integrated pressure and pressure gradient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can project the momentum equation (:eq:`eq-momentum-balance-approx`), using the volume between :math:`x_3` and the surface :math:`h`, in the normal
direction (:math:`\mathbf{v_3} = \mathbf{N^b} = -\mathbf{n^b}`). Applying the
properties of the NCS (:eq:`eq-du-in-ncs`) the surface normal
component of :eq:`eq-momentum-balance-approx` reads:

.. math::
	\begin{aligned}
	\rho_0 V(x_3, t) \frac{\mathrm{d}\overline{\mathbf{u}}(x_3)}{\mathrm{d}t} \cdot  \mathbf{v_3}  =&
	\rho_0 A^b (h- x_3) \overline{u}_1(x_3)\frac{\mathrm{d}\mathbf{v_1}}{\mathrm{d}t}
  \cdot\mathbf{v_3}\\
	=& -\rho_0 A^b (h - x_3) \overline{u}_1(x_3) \mathbf{v_1} \cdot \frac{\mathrm{d}\mathbf{v_3}}{\mathrm{d}t}
	= -\rho_0 A^b (h - x_3) \overline{\mathbf{u}}(x_3) \cdot \frac{\mathrm{d}\mathbf{N^b}}{\mathrm{d}t} \\
  =& - A^bp\,\cancelto{-1}{\mathbf{n^b}\cdot\mathbf{N^b}}
	+ A^b\cancelto{0}{\boldsymbol{\tau^b}\cdot\mathbf{N^b}}
	- A^b\boldsymbol{\nabla} \{(h-x_3)\overline{p}\}\cdot\mathbf{N^b}\\
	& + \rho_0 V \cancelto{g_{N^b}}{\mathbf{g}\cdot\mathbf{N^b}}
	+ \cancelto{0}{\mathbf{F}^{\text{ext}}\cdot\mathbf{N^b}}
	- \cancelto{0}{\overline{\mathbf{u}}\oint\limits_{\partial V(t)} q^{\text{ent}} \,\mathrm{d}A\cdot\mathbf{N^b}}
	\end{aligned}
  :label: eq-momentum-balance-x3-projected

Neglecting the normal component of the pressure gradient
gives the expression for pressure. Under the condition that :math:`\overline{\mathbf{u}}_1(x_3)` is
independent of :math:`x_3`, pressure follows a linear profile from the bottom
surface to the free surface:
Exploiting the normal component of the momentum equation enables to express the pressure and its gradient:

.. math::
  p(x_3) = \rho_0 (h - x_3) \left\{-g_{N^b}
  - \overline{\mathbf{u}}\cdot\frac{\mathrm{d}\mathbf{N^b}}{\mathrm{d}t}\right\}
  \quad \mbox{and}  \quad
  p(x_3=0) = p^b
  = \rho_0 h \left\{-g_{N^b}
  - \overline{\mathbf{u}}\cdot\frac{\mathrm{d}\mathbf{N^b}}{\mathrm{d}t}\right\}
  :label: eq-pressure-distribution

Note that the bottom pressure should always be positive.
A negative pressure is nonphysical and means that the material is not in contact with the bottom surface anymore.
This can happen in the case of large velocities on convex topography.
If so, the material should be in a free fall state until it gets back in contact with the topography.
A description on how this is handled within the numerical implementation can be found in Sect.~\ref{sec-curvature-acc-term-estimation}.

Using the previous result of :eq:`eq-pressure-distribution`, it is possible to express the thickness integrated
pressure :math:`\overline{p}`:

.. math::
   h\overline{p} = \int\limits_0^h p(x_3)\,\mathrm{d}x_3
	 = -\rho_0 \frac{h^2}{2}\left(g_{\mathbf{N^b}}
	 + \overline{\mathbf{u}} \cdot \frac{\mathrm{d}\mathbf{N^b}}
   {\mathrm{d}t}\right) = -\rho_0 \frac{h^2}{2} \, g^\text{eff}
   :label: eq-thickness-integrated-pressure

where :math:`g^\text{eff}` is the effective normal acceleration acting on the volume, including the normal component of
gravity and a curvature component.
The expression of the thickness integrated pressure is used to derive the pressure gradient :math:`\boldsymbol{\nabla} h\overline{p}`.
Assuming :math:`g^\text{eff}` to be locally constant (otherwise :math:`g^\text{eff}` would remain inside the gradient
operator), leads to:

.. math::
	\label{eq-pressure-gradient}
	 \boldsymbol{\nabla} h\overline{p} = -\rho_0 \, g^\text{eff} \, h \boldsymbol{\nabla} h
   :label: eq-pressure-gradient

Tangential momentum conservation equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the derived expression of the thickness integrated pressure
(:eq:`eq-pressure-gradient`), we project the momentum
balance (:eq:`eq-momentum-balance-approx`) in the tangent plane leading to the following equation:

.. math::
  \rho_0 V \left(\frac{\mathrm{d}\overline{\mathbf{u}}}{\mathrm{d}t}
  - \left(\frac{\mathrm{d}\overline{\mathbf{u}}}{\mathrm{d}t} \cdot
  \mathbf{v}_3\right)\mathbf{v}_3\right) =
  A^b\boldsymbol{\tau^b} - \rho_0 \, g^\text{eff} \, h A^b \boldsymbol{\nabla}_{s} h
  + \rho_0 V \mathbf{g}_s + \mathbf{F}^{\text{ext}}
	- \overline{\mathbf{u}}\oint\limits_{\partial V(t)} q^{\text{ent}} \,\mathrm{d}A
  :label: eq-momentum-balance-tangent

Where :math:`\boldsymbol{\nabla}_{s} = \boldsymbol{\nabla} - (\boldsymbol{\nabla}\cdot\mathbf{N^b})\mathbf{N^b}` respectively
:math:`\mathbf{g_s} = \mathbf{g} - (\mathbf{g} \cdot \mathbf{N^b})\mathbf{N^b}` is the
tangential component of the gradient operator respectively of the gravity
acceleration.

After replacing the velocity derivative component in the normal direction
by the expression developed in :eq:`eq-momentum-balance-x3-projected`,
:eq:`eq-momentum-balance-tangent` reads:

.. math::
    \begin{aligned}
    \rho_0 V \frac{\mathrm{d}\overline{\mathbf{u}}}{\mathrm{d}t} =&
    - A^bp\mathbf{n^b}
    + A^b\boldsymbol{\tau^b}
    - \rho_0 \, g^\text{eff} \, h A^b \boldsymbol{\nabla}_{\!s} h
    + \rho_0 V \mathbf{g}_s
    + \mathbf{F}^{\text{ent}}
    - \mathbf{F}^{\text{res}}\\
    &- \underbrace{\overline{\mathbf{u}}\,\rho_0\,\frac{\mathrm{d}\left(A^b\,h\right)}{\mathrm{d}t}}_{\text{from the mass balance}}
    - \underbrace{\rho_0 V \left( \overline{\mathbf{u}} \cdot \frac{\mathrm{d}\mathbf{v}_3}{\mathrm{d}t} \right)\mathbf{v}_3}_{\text{curvature acceleration}}
    \end{aligned}
    :label: momentum-balance6

Friction Model
~~~~~~~~~~~~~~~~~

The constitutive equation which
describes the basal shear stress tensor :math:`\boldsymbol{\tau^b}` as a function
of the flow state of the avalanche can be expressed:

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
