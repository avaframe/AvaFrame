Appendix: Module DFA-Kernel
============================

Governing Equations for the Dense Flow Avalanche
------------------------------------------------------

The governing equations of the dense flow avalanche are derived from the
incompressible mass and momentum balance on a Lagrange control volume
([Zw2000]_ [ZwKlSa2003]_).

Mass balance:
~~~~~~~~~~~~~~~

.. math::
    \frac{d}{dt} \int\limits_{V(t)} \rho_0 \,\mathrm{d}V = \rho_0 \frac{dV(t)}{dt} =
    \oint\limits_{\partial V(t)} q_{\text{ent}} \,\mathrm{d}A
    :label: mass-balance1

Where :math:`q_{\text{ent}}` represents the snow entrainment rate.

Momentum balance:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::
    \rho_0 \frac{d}{dt} \int\limits_{V(t)} u_i \,\mathrm{d}V = \oint\limits_{\partial V(t)}
    \sigma^{\text{tot}}_{ij}n_j \,\mathrm{d}A + \rho_0 \int\limits_{V(t)} g_i \,\mathrm{d}V, \quad i=(1,2,3)
    :label: momentum-balance1


We introduce the volume average of a quantity :math:`P(\mathbf{x},t)`:

.. math::
    \overline{P}(\mathbf{x},t) =  \frac{1}{V(t)} \int\limits_{V(t)} P(\mathbf{x},t) \,\mathrm{d}V
    :label: volume-average

and split the area integral into :

.. math::
   \oint\limits_{\partial V(t)} \sigma^{\text{tot}}_{ij}n_j \,\mathrm{d}A =
   \oint\limits_{\partial V(t)} \sigma_{ij}n_j \,\mathrm{d}A + F_i^{\text{ent}} + F_i^{\text{res}}, \quad i=(1,2,3)
   :label: area-integral

:math:`F_i^{\text{ent}}` represents the force required to break the
entrained snow from the ground and to compress it (since the dense-flow
bulk density is usually larger than the density of the entrained snow,
i.e. :math:`\rho_{\text{ent}}<\rho`) and :math:`F_i^{\text{res}}`
represents the resistance force due to obstacles (for example trees).
Which leads to in :eq:`momentum-balance1`:

.. math::
   \rho_0 \frac{dV(t) \overline{u}_i}{dt} = \rho_0 V \frac{d\overline{u}_i}{dt} +
   \rho_0 \overline{u}_i \frac{dV}{dt} = \oint\limits_{\partial V(t)} \sigma_{ij}n_j
   \,\mathrm{d}A + \rho_0 V g_i + F_i^{\text{ent}} + F_i^{\text{res}}, \quad i=(1,2,3)
   :label: momentum-balance2

Using the mass balance equation :eq:`mass-balance1`, we get:

.. math::
   \rho_0 V \frac{d\overline{u}_i}{dt} = \oint\limits_{\partial V(t)} \sigma_{ij}n_j \,\mathrm{d}A
   + \rho_0 V g_i  + F_i^{\text{ent}} + F_i^{\text{res}} - \overline{u}_i \oint\limits_{\partial V(t)} q_{\text{ent}} \,\mathrm{d}A, \quad i=(1,2,3)
   :label: momentum-balance3

Boundary conditions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The free surface is defined by :

    .. math:: F_s(\mathbf{x},t) = z-s(x,y,t)=0

The bottom surface is defined by :

    .. math:: F_b(\mathbf{x}) = z-b(x,y)=0

The boundary conditions at the free surface and bottom of the flow read:

.. math::
   \left\{\begin{aligned}
   &\frac{dF_s}{dt} = \frac{\partial F_s}{\partial t} +  u_i\frac{\partial F_s}{\partial x_i} =0 \quad & \mbox{at  }F_s(\mathbf{x},t) =0 \quad & \mbox{Kinematic BC (Material boundary)}\\
   &\sigma_{ij}n_j = 0 \quad & \mbox{at  }F_s(\mathbf{x},t) =0 \quad & \mbox{Dynamic BC (Traction free surface)}\\
   &u_in_i = 0 \quad & \mbox{at  }F_b(\mathbf{x},t) =0 \quad & \mbox{Kinematic BC (No detachment)}\\
   &\tau^{(b)}_i = f(\sigma^{(b)},\overline{u},\overline{h},\overline{\rho},t,\mathbf{x})\quad & \mbox{at  }F_b(\mathbf{x},t) =0\quad & \mbox{Dynamic BC (Chosen friction law)}
   \end{aligned}
   \right.
   :label: boundary-conditions

:math:`\sigma^{(b)}_i = (\sigma_{kl}n_ln_k)n_i` represents the normal stress at the bottom and
:math:`\tau^{(b)}_i = \sigma_{ij}n_j - \sigma^{(b)}_i` represents the shear stress at the bottom surface.
:math:`f` describes the chosen friction model and are described in `[sec:friction model] <#sec:friction model>`__.
The normals at the free surface (:math:`n_i^{(s)}`) and bottom surface (:math:`n_i^{(b)}`) are:

.. math::
   n_i^{(s,b)} = \frac{\partial F_{s,b}}{\partial x_i}\left(\frac{\partial F_{s,b}}{\partial x_j}
   \frac{\partial F_{s,b}}{\partial x_j}\right)^{-1/2}
   :label: surface-normals

Choice of the coordinate system:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The previous equations will be developed in the orthonormal coordinate
system :math:`(B,\mathbf{e_1},\mathbf{e_2},\mathbf{e_3})`, further
referenced as Natural Coordinate System (NCS). In this NCS,
:math:`\mathbf{e_1}` is aligned with the velocity vector at the bottom
and :math:`\mathbf{e_3}` with the normal to the slope, i.e.:

.. math::
   \mathbf{e_1} = \frac{\mathbf{u}}{\left\Vert \mathbf{u}\right\Vert},\quad \mathbf{e_2} = \mathbf{e_3}\wedge\mathbf{e_1},
   \quad \mathbf{e_3} = \mathbf{n^{(b)}}
   :label: natural-coordinate-system

The origin :math:`B` of the NCS is attached to the slope. This choice
leads to:

.. math::
   n^{(b)}_i = \delta_{i3}, \quad \left.\frac{\partial b}{\partial x_i}\right\rvert_{\mathbf{0}} = 0\quad
   \mbox{for} \quad i=(1,2),\quad \mbox{and} \quad u^{(b)}_2 = u^{(b)}_3 = 0
   :label: NCS-consequence

In this NCS and considering a prism-like Control volume, the volume
content :math:`V(t) = A_b(t)\overline{h}` is obtained by multiplication
of the basal area of the prism, :math:`A_b`, with the averaged value of
the flow-depth,

.. math::
    \overline{h} = \frac{1}{A_b(t)}\int\limits_{A_b(t)} [s(\mathbf{x})-b(\mathbf{x})]\,\mathrm{d}A = \frac{1}{A_b(t)}\int\limits_{A_b(t)} h(\mathbf{x})\,\mathrm{d}A,\qquad
    \overline{u}_i = \frac{1}{V(t)}\int\limits_{V(t)} u_i(\mathbf{x})\,\mathrm{d}V
    :label: hmean-umean

The snow entrainment is either due to plowing at the front of the avalanche or to erosion
at the bottom. The entrainment rate at the front :math:`q_{\text{plo}}` can be expressed function of the
properties of the entrained snow (density :math:`\rho_{\text{ent}}` and
snow depth :math:`h_{\text{ent}}`), the velocity of the avalanche at the
front :math:`\overline{\mathbf{u}}` and length :math:`w_f` of the front cell (measured perpendicularly
to the flow velocity :math:`\overline{\mathbf{u}}`). It obviously only happens on the front cells of
the avalanche (meaning that :math:`w_f` is zero for inner parts of the avalanche):

.. math::
   \oint\limits_{\partial V(t)} q_{\text{plo}}\,\mathrm{d}A = \int\limits_{l_{\text{front}}}\int_b^s q_{\text{plo}}\,
   \mathrm{d}{l}\,\mathrm{d}{z} =  \rho_{\text{ent}}\,w_f\,h_{\text{ent}}\,\left\Vert \overline{\mathbf{u}}\right\Vert
   :label: ploughing

The entrainment rate at the bottom :math:`q_{\text{ero}}` can be expressed function of the
bottom area :math:`A_b` of the control volume, the velocity of the avalanche :math:`\overline{\mathbf{u}}`,
the bottom shear stress :math:`\tau^{(b)}` and the specific erosion energy :math:`e_b`):

.. math::
    \oint\limits_{\partial V(t)} q_{\text{ero}}\,\mathrm{d}A = \int\limits_{A_b} q_{\text{ero}}\,
    \mathrm{d}A = A_b\,\frac{\tau^{(b)}}{e_b}\,\left\Vert \overline{\mathbf{u}}\right\Vert
    :label: erosion


This leads in the mass balance :eq:`mass-balance1` to :

.. math::
   \frac{\mathrm{d}V(t)}{\mathrm{d}t} = \frac{\mathrm{d}(A_b\overline{h})}{\mathrm{d}t}
   = \frac{\rho_{\text{ent}}}{\rho_0}\,w_f\,h_{\text{ent}}\,\left\Vert \overline{\mathbf{u}}\right\Vert +
   A_b\,\frac{\tau^{(b)}}{e_b}\,\left\Vert \overline{\mathbf{u}}\right\Vert
   :label: mass-balance2

The force :math:`F_i^{\text{ent}}` required to break the entrained snow
from the ground and to compress it is expressed function of the required
breaking energy per fracture surface unit :math:`e_s`
(:math:`J.m^{-2}`), the deformation energy per entrained mass element
:math:`e_d` (:math:`J.kg^{-1}`) and the entrained snow depth
[Sa2007]_ [SaFeFr2008]_ [FiFrGaSo2013]_:

.. math:: F_i^{\text{ent}} = -w_f\,(e_s+\,q_{\text{ent}}\,e_d)

The force :math:`F_i^{\text{res}}` due to obstacles is expressed
function of the characteristic diameter :math:`\overline{d}` and height
:math:`h_{\text{res}}` of the obstacles, the spacing
:math:`s_{\text{res}}` between the obstacles and an empirical
coefficient :math:`c_w` (see :numref:`f-res`):

.. math::
   F_i^{\text{res}} = -(\frac{1}{2}\,\overline{d}\,c_w/s^2_{\text{res}})\,\overline{\rho}\,A\,
    \min\left\lbrace\begin{array}{l} \overline{h}\\h_{res}\end{array}\right\rbrace\,\overline{u}^2\,
    \frac{\overline{u}_i}{\|\overline{u}\|}


.. _f-res:

.. figure:: _static/f_res.png
        :width: 90%

        Resistance force due to obstacles (from [FiKo2013]_)

The surface integral is split in three terms, an integral over
:math:`A_b` the bottom :math:`x_3 = b(x_1,x_2)`, :math:`A_s` the top
:math:`x_3 = s(x_1,x_2,t)` and :math:`A_h` the lateral surface.
Introducing the boundary conditions :eq:`boundary-conditions` leads to:

.. math::
   \begin{aligned}
   \oint\limits_{\partial{V(t)}}\sigma_{ij}n_j\,\mathrm{d}A & =
   \int\limits_{A_b}\underbrace{\sigma_{ij}\,n_j^{(b)}}_{-\sigma_{i3}}\,\mathrm{d}A +  \int\limits_{A_s}\underbrace{\sigma_{ij}\,n_j^{(s)}}_{0}\,\mathrm{d}A + \int\limits_{A_h}\sigma_{ij}\,n_j\,\mathrm{d}A\\
   &= -A_b\overline{\sigma}_{i3}^{(b)} + \oint\limits_{\partial A_b}\left(\int_b^s\sigma_{ij}\,n_j\,\mathrm{d}x_3\right)\,\mathrm{d}l
   \end{aligned}
   :label: surface forces

Which simplifies the momentum balance :eq:`momentum-balance3` to:

.. math::
   \rho_0 V \frac{d\overline{u}_i}{dt} = \oint\limits_{\partial A_b}\left(\int_b^s\sigma_{ij}\,n_j\,
   \mathrm{d}x_3\right)\,\mathrm{d}l -A_b\overline{\sigma}_{i3}^{(b)} + \rho_0 V g_i  + F_i^{\text{ent}} +
   F_i^{\text{res}} - \overline{u}_i \oint\limits_{\partial V(t)} q_{\text{ent}} \,\mathrm{d}A, \quad i=(1,2,3)
   :label: momentum-balance5

The momentum balance in direction :math:`x_3` (normal to the slope) is
used to obtain a relation for the vertical distribution of the stress
tensor [Sa2007]_. Due to the choice of
coordinate system and because of the kinematic boundary condition at the
bottom, the left side of :eq:`momentum-balance5` can be
expressed function of the velocity :math:`\overline{u}_1` in direction
:math:`x_1` and the curvature of the terrain in this same direction
:math:`\frac{\partial^2{b}}{\partial{x_1^2}}` [Zw2000]_:

.. math::
   \rho\,A_b\,\overline{h}\,\frac{\,\mathrm{d}\overline{u}_3}{\,\mathrm{d}t} =
   \rho\,A_b\,\overline{h}\,\frac{\partial^2{b}}{\partial{x_1^2}}\,\overline{u}_1^2,

rearranging the terms in the momentum equation leads to:

.. math::
  \overline{\sigma}_{33}(x_3) = \rho_0\,(s-x_3)\left(g_3-\frac{\partial^2{b}}{\partial{x_1^2}}\,\overline{u}_1^2\right)+ \frac{1}{A_b}
  \oint\limits_{\partial A_b}\left(\int_{x_3}^s\sigma_{3j}\,n_j\,\mathrm{d}x_3\right)\,\mathrm{d}l
  :label: sigma33

Non-dimensional Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _fig-charakterisitsche_groessen:

.. figure:: _static/charakterisitsche_groessen.png
        :width: 90%

        Characteristic size of the avalanche along its path (from [Zw2000]_)

The previous equations :eq:`momentum-balance5` and :eq:`sigma33` are can be further simplified by
introducing a scaling based on the characteristic values of the physical
quantities describing the avalanche. The characteristic length L, the
thickness H, the acceleration due to gravity g and the characteristic
radius of curvature of the terrain R are the chosen quantities. From
those values, it is possible to form two non dimensional parameters that
describe the flow:

-  Aspect ratio: :math:`\qquad\qquad\varepsilon = H / L\qquad`

-  Curvature:  :math:`\qquad\lambda = L / R\qquad`

The different properties involved are then expressed in terms of
characteristic quantities :math:`L`, :math:`H`, :math:`g` and :math:`R`
(see :numref:`fig-charakterisitsche_groessen`):

.. math::
   \begin{aligned}
    x_i &= L\, x_i^*\\
    (dx_3,h,\overline{h}) &= H\,(dx_3^*,h^*,\overline{h}^*)\\
    A_b &= L^2\, A_b^*\\
    t &= \sqrt{L/\text{g}}\, t^*\\
    \overline{u_i} &= \sqrt{\text{g}L}\,\overline{u_i}^*\\
    \text{g}_i &= \text{g} \, \text{g}_i^*\\
    \frac{\partial^2{b}}{\partial{x_1}^2} &= \frac{1}{R}\,\frac{\partial^2{b^*}}{\partial{x_1}^{*2}}\end{aligned}

The normal part of the stress tensor is directly related to the
hydro-static pressure:

.. math:: \sigma_{ii} = \rho\,\text{g}\,H\,\sigma_{ii}^*

The dimensionless properties are indicated by a superscripted asterisk.
Introducing those properties in :eq:`sigma33`, leads to
:

.. math::
   \overline{\sigma^*}_{33} = \left(g^*_3-\lambda\frac{\partial^2{b^*}}{\partial{x_1^{*2}}}\,\overline{u}_1^{*2}\right)
   (s^*-x^*_3) + \underbrace{\varepsilon\oint\limits_{\partial A_b^*}\left(\int\limits_{x^*_3}^{s^*}\sigma^*_{31}\,\mathrm{d}x^*_3\right)\,\mathrm{d}l^*}_{O(\varepsilon)}.
   :label: sigma33star

The height, H of dense flow avalanches is assumed to me small compared
to its length, L. Meaning that the equations are examined in the limit
:math:`\varepsilon \ll 1`. It is then possible to neglect the last term
in :eq:`sigma33star` which leads to (after reinserting
the dimensions):

.. math::
    \overline{\sigma}_{33}(x_3) = \rho\,\left(g_3-\overline{u_1}^2\,\frac{\partial^2{b}}{\partial{x_1^2}}\right)
    \left[\overline{h}-x_3\right]
    :label: sigma33dim

And at the bottom of the avalanche, with :math:`x_3 = 0`, the normal
stress can be expressed as:

.. math::
     \overline{\sigma}^{(b)}_{33} = \rho\,\left(g_3-\overline{u_1}^2\,\frac{\partial^2{b}}{\partial{x_1^2}}\right)\,\overline{h}
     :label: sigmab

Calculating the surface integral in equation :eq:`momentum-balance5` requires to
express the other components of the stress tensor. Here again a
magnitude consideration between the shear stresses :math:`\sigma_{12} = \sigma_{21}` and :math:`\sigma_{13}`.
The shear stresses are based on a generalized Newtonian law of materials,
which controls the influence of normal stresse and the rate of deformation through the viscosity.

.. math::
    \tau_{ij} = \eta\left(\frac{\partial{u_i}}{\partial{x_j}}+\frac{\partial{u_j}}{\partial{x_i}}\right), ~ i\neq j

Because :math:`\partial x_1` and :math:`\partial x_2` are of the order of :math:`L`, whereas :math:`\partial x_3`
is of the order of :math:`H`, it follows that:

.. math::
    O\left(\frac{\sigma_{12}}{\sigma_{13}}\right) = \frac{H}{L} = \varepsilon \ll 1

and thus :math:`\sigma_{12} = \sigma_{21}` is negligible compared to :math:`\sigma_{13}`.
:math:`\sigma_{13}` is exressed using the bottom friction law
:math:`\tau^{(b)}_i = f(\sigma^{(b)},\overline{u},\overline{h},\overline{\rho},t,\mathbf{x})`
introduced in :eq:`boundary-conditions`.


In addition, a relation linking the horizontal normal stresses,
:math:`\sigma_{ii}`, :math:`i = (1,2)`, to the vertical pressure distribution given
by :eq:`sigmab` is introduced. In complete analogy to the arguments used by
Savage and Hutter [SaHu1989]_ the horizontal normal stresses are given as:

.. math::
    \sigma_{ii} = K_{(i)}\,\sigma_{33}

Where :math:`K_{(i)}` are the earth pressure coefficients (cf. [ZwKlSa2003]_ [Sa2004]_):

.. math::
   \begin{aligned}
    \sigma_{11} &= K_{x~akt/pass}\,\sigma_{33}\\
    \sigma_{22} &= K_{y~akt/pass}^{(x~akt/pass)}\,\sigma_{33}\end{aligned}

With the above specifications, the integral of the stresses over the
flow height is simplified in equation :eq:`momentum-balance5` to:

.. math::
   \int\limits_b^s\sigma_{ij}\,\mathrm{d}x_3 = \int\limits_b^s K_{(i)}\,\sigma_{33}\,\mathrm{d}x_3 =
    K_{(i)}\,\frac{\overline{h}\,\sigma^{(b)}}{2}

and the momentum balance can be written:

.. math::
    \frac{\,\mathrm{d}\overline{u}_i}{\,\mathrm{d}t} =
    g_i + \frac{K_{(i)}}{\overline{\rho}\,A\,\overline{h}}\,\oint\limits_{\partial{A}}\left(\frac{\overline{h}\,\sigma^{(b)}}{2}\right)n_i\,\mathrm{d}l
    -\delta_{i1}\frac{\tau^{(b)}}{\overline{\rho}\,\overline{h}} - C_{\text{res}}\,\overline{\mathbf{u}}^2\,\frac{\overline{u_i}}{\|\overline{\mathbf{u}}\|}
    -\frac{\overline{u_i}}{A\,\overline{h}}\frac{\,\mathrm{d}(A\,\overline{h})}{\,\mathrm{d}t} + \frac{F_i^{\text{ent}}}{\overline{\rho}\,A\,\overline{h}}
    :label: momentum-balance6

with

.. math:: C_{\text{res}} = \frac{1}{2}\,\overline{d}\,\frac{c_w}{s_{\text{res}}^2}.

The mass balance :eq:`mass-balance2`
remains unchanged:

.. math::
   \frac{\mathrm{d}V(t)}{\mathrm{d}t} = \frac{\mathrm{d}(A_b\overline{h})}{\mathrm{d}t}
   = \frac{\rho_{\text{ent}}}{\rho_0}\,w_f\,h_{\text{ent}}\,\left\Vert \overline{\mathbf{u}}\right\Vert
   + A_b\,\frac{\tau^{(b)}}{e_b}\,\left\Vert \overline{\mathbf{u}}\right\Vert
   :label: mass-balance3

The unknown :math:`\overline{u}_1`, :math:`\overline{u}_2` and
:math:`\overline{h}` satisfy :eq:`sigmab`,
:eq:`momentum-balance6` and
:eq:`mass-balance3`. In equation
:eq:`momentum-balance6` the ground shear
stress :math:`\tau^{(b)}` remains unknown, and and a constitutive equation
has to be introduced in order to completely solve the equations.


Friction Model
~~~~~~~~~~~~~~~~~

The problem can be solved by introducing a constitutive equation which
describes the basal shear stress tensor :math:`\tau^{(b)}` as a function
of the flow state of the avalanche.

.. math::
    \tau^{(b)}_i = f(\sigma^{(b)},\overline{u},\overline{h},\rho_0,t,\mathbf{x})

With

.. math::
   \begin{aligned}
   \quad &\sigma^{(b)} \qquad \text{normal component of the stress tensor}\\
   \quad &\overline{u} \qquad\quad \text{average velocity}\\
   \quad &\overline{h} \qquad\quad \text{average flow depth}\\
   \quad &\rho_0 \qquad\quad \text{density}\\
   \quad &t \qquad\quad\, \text{time}\\
   \quad &\mathbf{x} \qquad\quad \text{position vector}
   \end{aligned}

Several friction models already implemented in the simulation tool are
described here.

SamosAT Model
~~~~~~~~~~~~~

SamosAT friction model is a modification of some more clasical models
such as Voellmy model. The basal shear stress tensor :math:`\tau^{(b)}`
is expressed as [Sa2007]_:

.. math::
   \tau^{(b)} = \tau_0 + \tan{\delta}\,\left(1+\frac{R_s^0}{R_s^0+R_s}\right)\,\sigma^{(b)}
    + \frac{\rho_0\,\overline{u}^2}{\left(\frac{1}{\kappa}\,\ln\frac{\overline{h}}{R} + B\right)^2}

With

.. math::
   \begin{aligned}
   \quad \tau_0 \qquad &\text{minimum shear stress}\\
   \quad R_s \qquad &\text{relation between friction and normal pressure (fluidization factor)}\\
   \quad R_s^0 \qquad &\text{empirical  constant}\\
   \quad B \qquad &\text{empirical  constant}\\
   \quad \kappa \qquad &\text{empirical  constant}\end{aligned}

The minimum shear stress :math:`\tau_0` defines a lower limit below
which no flow takes place with the condition
:math:`\rho_0\,\overline{h}\,g\,\sin{\alpha} > \tau_0`. :math:`\alpha`
being the slope. :math:`\tau_0` is independent of the flow depth, which
leeds to a strong avalanche deceleration, especially for avalanches with
low flow heights. :math:`R_s` is expressed as
:math:`R_s = \frac{\rho_0\,\overline{u}^2}{\sigma^{(b)}}`). Together
with the empirical parameter :math:`R_s^0` the term
:math:`\frac{R_s^0}{R_s^0+R_s}` defines the Coulomb basal friction.
Therefore lower avalanche speeds lead to a higher bed friction, making
avalanche flow stop already at steeper slopes :math:`\alpha`, than
without this effect. This effect is intended to avoid lateral creep of
the avalanche mass [SaGr2009]_


Numerics
========

Mass :eq:`mass-balance3` and momentum
:eq:`momentum-balance6` balance
equations as well as basal normal stress :eq:`sigmab`
are solved numerically using a SPH method (**S**\ moothed **P**\ article
**H**\ ydrodynamis) [Mo1992]_ for the three variables
:math:`\overline{\mathbf{u}}=(\overline{u}_1, \overline{u}_2)` and
:math:`\overline{h}` by discretization of the released avalanche volume
in a large number of mass elements. SPH in general, is a mesh-less
numerical method for solving partial differential equations. The SPH
algorithm discretizes the numerical problem within a domain using
particles [Sa2007]_ [SaGr2009]_, which interact
with each-other in a defined zone of influence. Some of the advantages
of the SPH method are that free surface flows, material boundaries and
moving boundary conditions are considered implicitly. In addition, large
deformations can be modeled due to the fact that the method is not based
on a mesh. From a numerical point of view, the SPH method itself is
relatively robust. Each particle has a mass :math:`m_{p,j}`. Those
particles are projected on a raster and the mass distributed on each
node of the raster (see :numref:`raster`). This leads to the
following expression for the mass :math:`m_{v,i}` of each node on the
raster grid: :math:`m_{v,i} = \sum\limits_j^{n_i}{m_{p,j}}`

.. _raster:

.. figure:: _static/raster.png
        :width: 70%

        Particles in raster grid (from [FiKo2013]_)

Each grid node is also affected with a velocity
:math:`\overline{\mathbf{u}}_{v,i}` expressed as the sum of the momentum
of each raster cell divided by the mass of the same cell:

.. math::
    \overline{\mathbf{u}}_{v,i} = \frac{\sum\limits_j^{n_i}{m_{p,j}}\mathbf{u}_j}{\sum\limits_j^{n_i}{m_{p,j}}}

The flow depth :math:`\overline{h}_{v,i}` can be deduced from the mass
and area of the raster cell:

.. math::
    \overline{h}_{v,i} = \frac{m_{v,i}}{\rho_0\,A_{v,i}}

The bottom area paired to each particle is related to the mass and flow
depth of this one:

.. math::
    A_{p,i} = \frac{m_{p,i}}{\rho_0\,\overline{h}_{p,i}}

The SPH method is introduced when expressing the flow depth of each mass
particle as a weighted sum of its neighbours
[LiLi2010]_ [Sa2007]_:

.. math::
    \overline{h}_{p,j} = \frac{1}{\rho_0}\,\sum\limits_j^{n_i}{m_{p,j}}\,W_{ij}

Where :math:`W` represents the SPH-Kernel function and reads:

.. math::
   W_{ij} = W(\mathbf{r_{ij}},r_0) = \frac{10}{\pi h^5}\left\{
   \begin{aligned}
   & (r_0 - \left\Vert \mathbf{r_{ij}}\right\Vert), \quad &0\leq \left\Vert \mathbf{r_{ij}}\right\Vert \leq  r_0\\
   & 0 , & r_0 <\left\Vert \mathbf{r_{ij}}\right\Vert
   \end{aligned}
   \right.
   :label: kernel function

:math:`\left\Vert \mathbf{r_{ij}}\right\Vert= \left\Vert \mathbf{x_j}-\mathbf{x_i}\right\Vert`
represents the distance between particle :math:`i` and :math:`j` and
:math:`r_0` the smoothing length. Now the lateral pressure forces on
each particle have to be determined. These are calculated from the
compression forces on the boundary of the particle. The boundary is
approximated as a square with the base side length
:math:`\Delta s = \sqrt{A_p}` and the respective flow height. This leads
to:

.. math::
    F_i = K_{(i)}\oint\limits_{\partial{A_p}}\left(\int\limits_{b}^{s}\sigma_{33}\,n_i\,\mathrm{d}x_3\right)\mathrm{d}l

From equation :eq:`momentum-balance6`

.. math::
    F_i = K_{(i)}\,\frac{\Delta{s}}{2}\left((\overline{h}\,\overline{\sigma}^{(b)}_{33})_{x_{p,i}-
    \frac{\Delta{s}}{2}}-(\overline{h}\,\overline{\sigma}^{(b)}_{33})_{x_{p,i}+\frac{\Delta{s}}{2}}\right)

The basal normal pressure :math:`\overline{\sigma}^{(b)}_{33}` from
equation :eq:`sigmab` is determined analogously to the
flow thicknesses :math:`\overline{h}_{p,i}`:

.. math::
   {\overline{\sigma}^{(b)}_{33}}_{p,i} = \sum\limits_j^{n_i}{{\overline{\sigma}^{(b)}_{33}}_{p,i}}\,A_{p,j}\,W_{ij}
     = \frac{1}{\rho_0}\,\sum\limits_j^{n_i}{{\overline{\sigma}^{(b)}_{33}}_{p,i}}\,\frac{m_{p,j}}{\overline{h}_{p,j}}\,W_{ij}

The mass of entrained snow for each particle :math:`p` between :math:`t`
and :math:`t+\Delta{t}` depends on the type of entrainment involved
(ploughing or erosion) and reads:

.. math::
    \Delta{m_p} = A_\text{ent}\,q_{\text{ent}}\,\Delta{t}

with

.. math::
   \begin{aligned}
    A_{\text{plo}} &= w_f\,h_{\text{ent}}= \sqrt{\frac{m_{p}}{\overline{\rho}\,\overline{h}_{p}}}\,h_{\text{ent}}
    \quad &\mbox{and} \quad &q_{\text{plo}} = \rho_{\text{ent}}\,\left\Vert \overline{\mathbf{u}}\right\Vert
    \quad &\mbox{for ploughing}\\
    A_{\text{ero}} &= A_p = \frac{m_{p}}{\overline{\rho}\,\overline{h}_{p}}
    \quad &\mbox{and} \quad &q_{\text{ero}} = \frac{\tau^{(b)}}{e_b}\,\left\Vert \overline{\mathbf{u}}\right\Vert
    \quad &\mbox{for erosion}
    \end{aligned}

Descretizing the momentum balance
:eq:`momentum-balance6` in time enables
to write the the velocity of the particle at the next time step:

.. math::
   u_i^{k+1} = \frac{u_i^k + \Delta{t}\,\left(g_i+\frac{F_i+F_i^\text{ent}}{m_p}\right)}
    {1 + \Delta{t}\left(\frac{\tau^{(b)}}{\bar\rho\,\overline{h}\,\|\overline{u}\|^k}+C_\text{res}\,\|\overline{u}\|^k\right)}
    -u_i^k\,\frac{m_p}{m_p+\Delta{m}_p}

With

.. math::
   \begin{aligned}
   \quad k \qquad &\text{current time step,}\\
   \quad k+1 \qquad &\text{next time step,}\\
   \quad \Delta{t} \qquad &\text{time step size.}
   \end{aligned}

The new position of the particle (in the next time step :math:`k+1`)
reads:

.. math::
    X_i^{k+1} = X_i^k + \frac{\Delta{t}}{2}(u_i^k + u_i^{k+1})


.. _fig-infinitesimales_element:

.. figure:: _static/infinitesimales_element.png
        :width: 90%

        Infinitesimal volume element and acting forces on it (from [FiKo2013]_)


References
----------

.. [FiFrGaSo2013] J. T. Fischer and R. Fromm and P. Gauer and B. Sovilla. (2013)
  Evaluation of probabilistic snow avalanche simulation ensembles with Doppler radar observations. Cold Regions Science and Technology.

.. [FiKo2013] J. T. Fischer and A. Kofler. (2013)
    SamosAT CoSiCa. Concepts for enhanced Simulation and multivariate Calibration. BFW

.. [LiLi2010] M.B. Liu and G.R. Liu. (2010).
    Smoothed Particle Hydrodynamics (SPH): an Overview and Recent Developments. Arch Computat Methods Eng 17, 25--76.

.. [Mo1992] J.J. Monaghan. (1992).
      Smoothed particle hydrodynamics. Annual review of astronomy and astrophysics. Vol 30. 543--574.

.. [Sa2004] B. Salm. (2004).
    A short and personal history of snow avalanche dynamics. Cold Regions Science and Technology. Vol. 39. 83--92.

.. [Sa2007] P. Sampl. (2007).
    SamosAT Modelltheorie und Numerik. AVL List GMBH.

.. [SaFeFr2008] R. Sailer and W. Fellin and R. Fromm and P. J{\"o}rg and L. Rammer and P. Sampl and A. Schaffhauser. (2008).
    Snow avalanche mass-balance calculation and simulation-model verification. Annals of Glaciology. Vol. 48, 183--192.

.. [SaGr2009] P. Sampl and M. Granig. (2009).
    Avalanche simulation with SAMOS-AT. Proceedings of the International Snow Science Workshop, Davos.

.. [SaHu1989] S.Savage and K. Hutter. (1989).
    The motion of a finite mass of granular material down a rough incline. Journal of Fluid Mechanics, 199, 177-215.

.. [Zw2000] T. Zwinger. (2000).
    Dynamik einer Trockenschneelawine auf beliebig geformten Berghangen, Technischen Universitaet Wien.

.. [ZwKlSa2003] T. Zwinger and A. Kluwick and P. Sampl. (2003).
    Numerical simulation of dry-snow avalanche flow over natural terrain.
    In: Hutter K., Kirchner N. (eds) Dynamic Response of Granular and Porous Materials under Large and Catastrophic Deformations.
    Lecture Notes in Applied and Computational Mechanics, vol 11. Springer, Berlin, Heidelberg.
