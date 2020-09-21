com1DFA: Module DFA-Kernel
===========================

This module is a simulation tool for dense flow (snow) avalanches (DFA). It is based on the
samos-AT (Snow Avalanche Modeling and  Simulation- Advanced Technologies) developed by the Austrian government
in cooperation with the company AVL List GmbH in Graz.
The calculation of the DFA is based on the depth integrated governing equations and 
solved numerically using the smoothed particle hydrodynamic (sph) method.

Input
---------

Output
---------

To run
--------

Theory
--------


The governing equations of the dense flow avalanche are derived from the
incompressible mass and momentum balance on a Lagrange control volume
([Zw2000]_ [ZwKlSa2003]_). Assuming the avalanche is much longer and larger
than thick, it is possible to integrate the governing equations over the thickness
of the avalanche and operate some simplifications due to the shape of the avalanche.
This leads, after some calculation steps described in details in Appendix
:ref:`appendixCom1DFA:Governing Equations for the Dense Flow Avalanche` to:

.. math::
    :label: mass-balance3

    \frac{\mathrm{d}V(t)}{\mathrm{d}t} = \frac{\mathrm{d}(A_b\overline{h})}{\mathrm{d}t}
    = \frac{\rho_{\text{ent}}}{\rho_0}\,w_f\,h_{\text{ent}}\,\left\Vert \overline{\mathbf{u}}\right\Vert


.. math::
    \frac{\,\mathrm{d}\overline{u}_i}{\,\mathrm{d}t} =
    g_i + \frac{K_{(i)}}{\overline{\rho}\,A\,\overline{h}}\,\oint\limits_{\partial{A}}\left(\frac{\overline{h}\,\sigma^{(b)}}{2}\right)n_i\,\mathrm{d}l
    -\delta_{i1}\frac{\tau^{(b)}}{\overline{\rho}\,\overline{h}} - C_{\text{res}}\,\overline{\mathbf{u}}^2\,\frac{\overline{u_i}}{\|\overline{\mathbf{u}}\|}
    -\frac{\overline{u_i}}{A\,\overline{h}}\frac{\,\mathrm{d}(A\,\overline{h})}{\,\mathrm{d}t} + \frac{F_i^{\text{ent}}}{\overline{\rho}\,A\,\overline{h}}
    :label: momentum-balance6


.. math::
    \overline{\sigma}^{(b)}_{33} = \rho\,\left(g_3-\overline{u_1}^2\,\frac{\partial^2{b}}{\partial{x_1^2}}\right)\,\overline{h}
    :label: sigmab
