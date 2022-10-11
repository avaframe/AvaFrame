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
  
   Dense Flow Standard 
      Runs dense flow avalanche module com1DFA. For more specific info about the inputs (shapefile attributes etc), see 
      :ref:`moduleCom1DFA:Input`. You can select multiple shapefiles for release areas, each will be calculated as one scenario. 
      Each shapefile can contain multiple polygons. If you provide entrainment/resistance areas, each scenario will be calculated 
      once without them (i.e null simulation) and once with them. 

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

Admin
-----

.. glossary::
   :sorted:
  
   GetVersion 
      Displays the current version of the AvaFrame python package (NOT the AvaFrameConnector)

   UpdateAvaFrame
      Updates the AvaFrame python package
