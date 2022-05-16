##################################
out1Peak: Peak plots
##################################



Plot all peak fields
====================

:py:func:`out1Peak.outPlotAllPeak.plotAllPeakFields` generates one plot for each
peak field in ``Outputs/modName/peakFiles``, constrained to existing data
including a buffer specified in ``out3Plot/plotUtilsCfg.ini``. These peak fields
represent the peak values of the simulation result parameters (*dynamic peak
pressure*, *peak flow thickness*, *peak velocity*), and modName corresponds to the
name of the computational module that has been used to perform the simualtions.
Details on this function, as for example required inputs, can be found in
:py:func:`out1Peak.outPlotAllPeak.plotAllPeakFields`.


Plot all fields
====================

:py:func:`out1Peak.outPlotAllPeak.plotAllFields` generates one plot for each
simulation result field provided in the specified input directory. This function
is designed to work for result fields that follow a certain naming convention in
order to provide full functionality:
*releaseAreaName_simulationType_modelType_simulationIdentifier_resultType.asc*

One plot for each field is saved using the name of the input field. By default,
the plot is constrained to where there is data, however this can be turned off
by setting ``constrainData=False`` in the inputs. Details on this function, as
for example required inputs, can be found in
:py:func:`out1Peak.outPlotAllPeak.plotAllField`.
