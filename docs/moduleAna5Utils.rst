
#####################################################################
ana5Utils
#####################################################################


distanceTimeAnalysis: Visualizing the temporal evolution of flow variables
------------------------------------------------------------------------------

With the functions gathered in this module, flow variables of avalanche simulation results can be
visualized in a distance versus time diagram, the so called **thalweg-time diagram**.
The **tt-diagram** provides a way to identify main features of the temporal evolution of
flow variables along the :term:`avalanche path`.
This is based on the ideas presented in :cite:`FiFrGaSo2013` and :cite:`RaKo2020`, where
avalanche simulation results have been transformed into the radar coordinate system to facilitate
direct comparison, combined with the attempt to analyze simulation results in an avalanche path
dependent coordinate system (:cite:`Fi2013`).
In addition to the **tt-diagram**, :py:mod:`ana5RangeTime` also offers the possibility to
produce simulated **range-time diagrams** of the flow variables with respect to a radar's field
of view. With this, simulation results can be directly compared to radar measurements in terms
of front position and inferred approach velocity. The colorcoding of the simulated
**range-time** diagrams refers to average values of the chosen flow parameter
(e.g. flow depth, flow velocity) at specified range gates. This colorcoding is not directly
comparable to the intensity given in the range-time diagram from radar measurements.

.. Note::
  The **tt-diagram** and the **range-time diagram** can be produced *during run time* of
  :py:mod:`com1DFA`, or as a postprocessing step. However, the second option requires saving all the
  required time steps of the flow parameter fields, which is much more computationally expensive
  compared to the first option.

To run
~~~~~~~

During run-time of :py:mod:`com1DFA`:

* in your local copy of ``com1DFA/com1DFACfg.ini`` in [VISUALISATION] set `createRangeTimeDiagram`
  to True and chose if you want a *TTdiagram* by setting this flag to True or in the case of a
  simulated range-time diagram to False

* in your local copy of ``ana5Utils/distanceTimeAnalysisCfg.ini`` you can adjust the default settings
  for the generation of the diagrams

As a postprocessing step:

* have a look at :py:mod:`runScripts.runThalwegTimeDiagram.py` and :py:mod:`runScripts.runRangeTimeDiagram.py`

* in your local copy of ``ana5Utils/distanceTimeAnalysisCfg.ini`` you can adjust the default settings
  for the generation of the diagrams

The resulting figures can be found in ``avalancheDirectory/Outputs/ana5Utils``.


.. figure:: _static/thalwegTime_FD.png
    :width: 90%

    Thalweg-time diagram example.


.. Note::
  The **tt-diagram** requires info on an avalanche path (see :ref:`moduleAna3AIMEC:ana3AIMEC: Aimec`).
  The simulated **range-time diagram** requires info on the x, y coordinate of the radar location, a point
  in the direction of the field of view, the aperture angle and the width of the range gates.
  The maximum approach velocity indicatd in the diagrams is computed as the ratio of the
  distance travelled by the front and the respective time needed for two increments along the path/
  radar's field of view.


Theory
~~~~~~~~~

Thalweg-time diagram
~~~~~~~~~~~~~~~~~~~~~~

First, the flow variable result field is transformed into a path-following coordinate system.
For this step, functions from :py:mod:`ana3AIMEC` are used.
The distance to the *start of runout area point* is determined using a user defined threshold
of the flow variable, giving the avalanche front position. The front positions defined with this
method for all the time steps are shown as black dots in the **tt-diagram**.
The mean values of the flow variable are computed at cross profiles along the avalanche path for
each time step and included in the **tt-diagram** as colored field.
For this analysis, all available flow variables can be chosen, but the interpretation of the
tt-diagram structures and the corresponding meaning of avalanche front changes comparing flow thickness or flow velocity.

Simulated Range-Time diagram
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The radar's field of view is determined using its location, a point in the direction of the field of
view and the aperture angle. The distance of every point in the simulation results to the radar
location is computed. The simulation results which lie outside the radar's field of view are masked.
The distance of the avalanche front with respect to the radar location is determined for a user
defined threshold in the flow variable and the average values of range gates along the radar's line
of sight are computed from the result field.
This data is plotted in a range-time diagram, where the black dots indicate the avalanche front,
and the colored field indicates the mean values of the flow parameter for the range gates for each
time step.





-
