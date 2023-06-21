AvaFrameConnector (GUI)
=======================

The AvaFrameConnector allows QGis users to access certain base workflows directly from QGis. The connector 
only provides the interface to functions within the AvaFrame python package, and is developed separately, see the
`QGisAF github repository <https://github.com/avaframe/QGisAF>`_. 
It makes use of the QGis processing plugin, which is included in all current QGis releases. 

Here is a quick overview of the processing scripts the AvaFrameConnector provides:

Operational
-----------

.. glossary::
   :sorted:
  
   Dense Flow Standard (com1)
      Runs dense flow avalanche module com1DFA. For more specific info about the inputs (shapefile attributes etc), see 
      :ref:`moduleCom1DFA:Input`. You can select multiple shapefiles for release areas, each will be calculated as one scenario. 
      Each shapefile can contain multiple polygons. If you provide entrainment/resistance areas, each scenario will be calculated 
      once without them (i.e null simulation) and once with them. 
   
   Alpha Beta (com2) 
      Runs the alpha beta calculation via module com2AB. For more specific info about the inputs (shapefile attributes etc), see 
      :ref:`moduleCom2AB:Input`. 

   Full Operational Run
      Runs com1DFA and additionally com2AB if profile and splitpoint are set. 
      Additional info see description of runCom1DFA. 


Experimental
------------

These either require deeper knowledge about the AvaFrame python package, i.e. configuration files etc., or are not yet 
fully tested and might produce unwanted results. 

.. glossary::
   :sorted:
  
   AvaFrameLayerRename   
      Renames com1DFA result layers by adding the values of the given variable (from the configuration file) 

   Probability run (ana5, com1)
      Runs probability simulations via module com1DFA. The release shape HAS TO HAVE a ci95 field containing the 
      95 percentile confidence interval (same unit as the release thickness). Multiple scenarios can be provided, 
      final map includes variations from all scenarios combinded. Release thickness and SamosAT friction mu 
      are being varied, 40 variations per scenario. 

   Release area stats(in1, com1)
      Returns info for release area statistics in a csv file. See output of QGis processing for the location of
      the file. 

   Snow slide (com5)
      Runs snow slide simulations via module com1DFA. For more info see com5SnowSlide section in the documentation.
      The resistance layer is meant for building outlines.


Admin
-----

.. glossary::
   :sorted:
  
   GetVersion 
      Displays the current version of the AvaFrame python package (NOT the AvaFrameConnector)

   UpdateAvaFrame
      Updates the AvaFrame python package
