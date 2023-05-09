Glossary
========

.. Note:: **This is still under construction**

.. Note:: This glossary has the main purpose to introduce and consolidate an
   unambiguous language and definitions with respect to avalanche dynamics,
   avalanche simulations and modelling for the digital toolbox AvaFrame and
   particularly the Avalanche Modelling Atlas (AMA). Terms are organized in
   alphabetical order and the terminology builds on the UNESCO Avalanche atlas
   :cite:`DQ1981` and relates to existing guidelines provided and used by
   Avalanche Associations and Warning Services (Canadian Avalanche Association 
   :cite:`CAA_OGRS2016` , American Avalanche Association (AAA) 
   :cite:`AAA_SWAG2016` or the European Avalanche Warning Services (EAWS)). 
   The glossary has the goal to follow common conventions used in publications
   and related projects.

.. glossary::
    :sorted:

    alpha point
    alpha angle
        The alpha point is the point of the furthest reach of an avalanche
        (visible deposition or affected area). The alpha angle is the angle between the
        line connecting the release point and the alpha point and the
        horizontally projected :term:`travel length` along the :term:`thalweg`
        (see :term:`travel length` and :term:`runout angle`). The alpha angle is also
        termed the angle of reach or the :term:`runout angle`. 

    beta point
    beta angle
        The boundary between the :term:`transition` and :term:`deposition` zone (start of
        deposition) is often called beta point and associated to a certain
        slope threshold of the :term:`thalweg`. The corresponding :term:`travel angle` is
        referred  to as beta angle. For major avalanche :term:`path` (s) that may produce
        avalanches of :term:`size` 4-5, the beta point is associated to the 10° Point
        (or beta_10), i.e. the point where the slope angle of the :term:`thalweg`
        decreases below 10°. For avalanche :term:`path` (s) with :term:`event` (s) of avalanche :term:`size`
        1-3 the beta point may accordingly be associated to larger :term:`thalweg`
        slope angles (i.e. up to 30°, see start of :term:`transition`).

    coordinate transformation
        Coordinate transformations refer to the operation of changing
        coordinates, e.g. between a fixed, Eularian, global coordinate system
        with east and north orientation, to an avalanche :term:`path` dependent
        coordinate system along the :term:`thalweg`, or even a Langrangian coordinate
        system, moving along particle trajectories. 

    projection
        The projection refers to
        the projection of coordinates within a coordinate system, I.e. the
        projection of the :term:`runout point` (as furthest reach of the avalalanche)
        to the :term:`thalweg`. Or projecting a 3D travel length (xyz) to a 2D
        travel length measured only along xy.

    cycle
        An avalanche cycle describes a series of :term:`event` (s) that occur across a region
        over a relatively short time span.

    danger scale
        The avalanche danger scale refers to the avalanche hazard and is an
        inherrent part of avalanche warning.

    dense flow
        Dense flow is a :term:`form of movement` in the :term:`transition` zone of the
        avalanche. Dense flow avalanches (DFA) are flowing along the ground and are rather
        associated to warm flow. Mixed types of movement are often observed,
        combining different flow regimes and their partial or complete
        transitions, e.g. ‘mixed flow and powder avalanches’ or ‘flow avalanche
        with powder component’, towards the evolution of a fluidized layer in
        the avalanche flow.

    depth
        Release, :term:`entrainment`, flow, or deposition depth refers to the extent of
        the avalanche measured in the direction of gravity.

    density
        Release, :term:`entrainment`, flow, or deposition density. Important quantity 
        relating mass and volume, influencing impact pressure and particular friction relations.

    entrainment
        Entrainment describes the process of mass intake during the avalanche flow.

    event
        see :term:`scenario`
    
    scenario
        One or multiple avalanche events or corresponding simulation scenarios
        are associated to a certain avalanche :term:`path` and have distinct criteria
        and characteristics such as avalanche :term:`size`, :term:`release area` (s),
        or :term:`runout area`. These properties are
        morphologically connected to the different zones (:term:`origin`, :term:`transition`,
        :term:`deposition`) of an avalanche :term:`path` and allow to define other
        associated properties, such as :term:`alpha angle` or :term:`runout length` that are
        defined in combination with the avalanche :term:`thalweg`. Besides observed and
        documented avalanche events, design events of particular :term:`return period` (s)
        are of particular interest for engineering applications.

    flow variables
        Flow variables include :term:`thickness`, :term:`velocity`, or :term:`density` and are
        determined by the form of movement. Flow models that are implemented
        usually calculate the spatio-temporal evolution of these variables and
        where the maximum over the whole flow or computational duration, I.e.,
        their peak values are the most used results. The flow variables are
        used to derive other variables such as impact pressure or kinetic
        energy of the flow.

    form of movement
        Is an avalanche criterion in the zone of :term:`transition` and has :term:`powder snow`
        or :term:`dense flow` as characteristics.

    manner of starting
        Is an avalanche criterion in the zone of :term:`origin` and has the possible
        characteristics loose, slab, or gliding.

    terrain classification
        Terrain may be classified according to the Avalanche Terrain Exposure
        Scale (ATES) into simple (low angle or primarily forested terrain with
        some openings that may involve the :term:`deposition` zone of infrequent avalanche
        :term:`path` (s)), challenging (well defined avalanche :term:`path` (s), starting zones, or
        terrain traps), and complex (exposure to multiple overlapping avalanche
        :term:`path` (s), large expanses of steep, open terrain, multiple starting zones,
        and terrain traps below).

    thalweg
        The thalweg is defined by the main flow direction of an avalanche
        :term:`path` of one or multiple,  i.e. not regarding a specific
        :term:`event` or :term:`scenario`. avalanche :term:`event` (s).
        Technically it is the two-dimensional
        terrain representation, displaying the terrain altitude along the
        horizontally projected :term:`travel length`.

    thickness
        Release, :term:`entrainment`, flow, or deposition thickness refers to the
        extent (distance) of the avalanche measured perpendicular to the slope.

    travel angle
    travel length
        Travel lengths are measured as horizontally projected travel length
        (:math:`s_{XY}`) along the :term:`thalweg` and are associated with the
        corresponding travel angle, measured between the line connecting the current
        location with the uppermost point of the release and the horizontal
        plane. Alternatively, the surface parallel travel length
        (:math:`s_{XYZ}`) may be defined as the three-dimensional length
        travelled by the avalanche.

    path
        The avalanche path summarizes the total catchment and is divided into
        different zones (zone of :term:`origin`, :term:`transition`, :term:`deposition`) 
        with different criteria and characteristics. An inherent property of the avalanche path is the
        :term:`thalweg` and the associated avalanche :term:`event`.

    powder snow
        Powder snow avalanches (PSA) refer to the :term:`form of movement` in the zone of
        :term:`transition`, referring to the dust or suspension cloud in avalanches.
        PSA are associated with cold, dry cohesionless snow. Mixed types of
        movement are often observed, combining different flow regimes and their
        transitions, e.g., ‘powder avalanche with flow component’.

    return period
        Return periods are related to return levels describing the :term:`size` or
        magnitude of design or recorded :term:`event` (s) on a respective avalanche :term:`path`.
        The return level is often determined by the run out length of
        historically documented avalanche :term:`event` (s) accompanied with return period
        estimates, which are associated to the occurrence probability.

    release area
        Potential release areas are located in the zone of :term:`origin`. Each documented
        :term:`event` or simulation :term:`scenario` is associated to one or more primary 
        and/or secondary release areas, that can further be described by the :term:`manner of starting`.

    runout area
    runout angle
    runout length
    runout point
        Runout lengths and angles are intricately linked to the :term:`alpha point`, utilizing the
        :term:`projection` to the :term:`thalweg`. In the same manner as :term:`travel length`, run out lengths
        are measured as horizontally projected lengths along the :term:`thalweg`, from
        the uppermost point of the :term:`release area` to furthest reach of the runout
        area. The runout may refer to visible deposition (associated to dense
        flow), damages or the impacted and affected area (associated to air
        blast or :term:`powder snow`) in the zone of :term:`deposition` and is usually defined
        via flow :term:`thickness`, velocity, kinetic energy or impact pressure
        thresholds.

    size
        Avalanche size refers to the magnitude or intensity of an :term:`event`, classified by destructive
        potential, :term:`runout length` and dimension according to the EAWS size
        classification, which is closely related to the CAA destructive size.

    velocity
        Flow velocities are usually measured in a surface parallel direction.
        Alternatively approach velocities are measured along the line of sight.

    wet snow
        A wet snow avalanche (WSA) implies the presence of liquid water within
        the avalanche and is usually associated to :term:`dense flow` type of movement
        in the :term:`transition` zone of the avalanche.

    origin
        see :term:`zone of origin`

    zone of origin
        The zone of origin delineates the area, in which typical :term:`release area` (s) are located, and an avalanche's
        appearance is characterized by the :term:`manner of starting`. The uppermost
        possible point is referred to as start of origin.
    
    transition
        see :term:`zone of transition`

    zone of transition    
        The zone of transition is the area between zone of :term:`origin` and
        zone of :term:`deposition` along the :term:`thalweg`. The :term:`form
        of movement` is linked to the :term:`flow variables`. The start of
        transition links the zone of :term:`origin` and transition and is usually
        associated with a slope inclination of about 28-30°.
    
    deposition
        see :term:`zone of deposition`

    zone of deposition
        The zone of deposition is where the :term:`runout area` of the avalanche is located and where the
        avalanche stops due to frictional energy dissipation. The boundary with
        the :term:`transition` zone (start of deposition) is often called :term:`beta point`.
