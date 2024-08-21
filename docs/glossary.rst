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
   :cite:`AAA_SWAG2016` or the European Avalanche Warning Services (:cite:`EAWS.2024`)).
   The glossary has the goal to follow common conventions used in publications
   and related projects.

.. glossary::
    :sorted:

    alpha point
        see :term:`runout point`

        The alpha point describes the furthest point of the runout of an avalanche :term:`event` 
        and marks the end of the avalanche deposition.

    alpha angle
        see :term:`runout angle`

        The alpha angle is the angle between the :term:`origin point` and the :term:`runout
        point`. It describes the general inclination between the release of the avalanche to the
        furthest point of the deposition. The alpha angle is based on the :term:`runout length`
        and the altitude difference between the :term:`origin point` and the :term:`runout point`.
        It is an important element in avalanche simulation and modeling (e.g. alpha-beta model,
        FlowPy simulations). The alpha angle is also termed the angle of reach or the :term:`runout angle`.

    beta point
        see :term:`deposition point`

        The boundary between the :term:`transition` and :term:`deposition` zone (start of
        deposition) is often called beta point or :term:`deposition point` and associated
        with a certain slope threshold of the :term:`thalweg`. For major avalanche :term:`path` that may produce
        avalanches of :term:`size` 4-5, the beta point is associated with the 10° Point
        (or beta_10), i.e., the point where the slope angle of the :term:`thalweg`
        decreases below 10°. For avalanche :term:`path` with :term:`event` of avalanche :term:`size`
        1-3, the beta point may accordingly be associated with larger :term:`thalweg`
        slope angles (i.e., up to 30°, see start of :term:`transition`).

    beta angle
        see :term:`travel angle`

        The beta angle describes the general incliniation along the :term:`thalweg` between the 
        :term:`origin point` and the :term:`deposition point`. The corresponding :term:`travel angle` is
        referred to as beta angle. It is an important concept in avalanche modeling and simulation. 
        The travel angle increases with decreasing avalanche :term:`size`.

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
        or :term:`runout area`. These properties are morphologically connected to the 
        different zones (:term:`zone of origin`, :term:`zone of transition`, :term:`zone of deposition`)
        of an avalanche :term:`path` and allow to define other associated properties, such as :term:`alpha angle`
        or :term:`runout length` that are defined in combination with the avalanche :term:`thalweg`.
        Besides observed and documented avalanche events, design events of particular 
        :term:`return period` (s) are of particular interest for engineering applications.

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
        The thalweg is defined as the line representing the main flow direction of all potential avalanche
        events within a specific avalanche :term:`path`. The thalweg is delineated according to the terrain 
        characteristics and is independent of a specific :term:`event`. Technically it is the two
        dimensional terrain representation along the horizontally projected :term:`travel length` and the 
        altitude difference. 

    thickness
        Release, :term:`entrainment`, flow, or deposition thickness refers to the
        extent (distance) of the avalanche measured perpendicular to the slope.

    trajectory length
        Used in com1DFA particle dictionaries, where the trajectory length is computed as the
        distance traveled by a particle from one time step to the next and then accumulated over
        time. Three different trajectory lengths are computed (1) trajectoryLengthXY - computed in the
        x, y plane, (2) trajectoryLengthXYZ - also taking the slope of the topography into account and
        (3) trajectoryLengthXYCor - same as trajectoryLengthXY but corrected for the potential
        angle difference of the slope and the normal.

    travel angle
        The travel angle describes the general inclination between the :term:`origin point` and the current location of interest. 
        It is calculated based on the :term:`travel lenght` and the altitude difference between the :term:`origin point` and the
        point of interest. Important travel angles are the :term:`beta angle` and the :term:`alpha angle`. 

    travel length
        Travel lengths are measured as horizontally projected travel length
        (:math:`s_{XY}`) along the :term:`thalweg`, measured between the current
        location with the uppermost point of the release, :term:`origin point`.
        Alternatively, the surface parallel travel length (:math:`s_{XYZ}`) may be
        defined as the three-dimensional length travelled by the avalanche.

    path
        The avalanche path summarizes the total catchment and is divided into
        different zones (:term:`zone of origin`, :term:`zone of transition`, :term:`zone of deposition`)
        with different criteria and characteristics. An inherent property of the avalanche path is the
        :term:`thalweg` and the associated avalanche :term:`event`.

    powder snow
        Powder snow avalanches (PSA) refer to the :term:`form of movement` in the 
        :term:`zone of transition`, referring to the dust or suspension cloud in avalanches.
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
        Potential release areas are located in the :term:`zone of origin`. Each documented
        :term:`event` or simulation :term:`scenario` is associated to one or more primary
        and/or secondary release areas, that can further be described by the :term:`manner of starting`.

    runout area
        see :term:`zone of deposition`

    runout angle
        The runout angle, also referred to as :term:`alpha angle`, describes the general inclination of the
        the avalanche event from the :term:`origin point` to the :term:`runout point`.

    runout length
        Runout lengths is intricately linked to the :term:`alpha point`, utilizing the
        :term:`projection` to the :term:`thalweg`. In the same manner as :term:`travel length`, run out lengths
        are measured as horizontally projected lengths along the :term:`thalweg`, from
        the uppermost point of the :term:`release area` (:term:`origin point`) to furthest reach of the runout
        area (:term:`runout point`).

    runout point
        The runout point is also referred to as :term:`alpha point`. It describes the furthest point of the runout,
        and marks the outer most end of the avalanche deposit. The runout may refer to visible deposition
        (associated to dense flow), damages or the impacted and affected area (associated to air
        blast or :term:`powder snow`) in the zone of :term:`deposition` and is usually defined
        via flow :term:`thickness`, velocity, kinetic energy or impact pressure thresholds.

    destructive size
        The destructive size refers to the magnitude or intensity of an :term:`event`, according to the EAWS size classification, which is closely
        related to the CAA destructive size. Thereby the size refers to the destructive potential, :term:`runout length`,  volumn or mass 
        of each avalanche event. 

    relative size
        The relative size refers to the size of an avalanche in relation to the :term:`thalweg`. The classification scheme was developed 
        by the AAA and relys on the horizontal extent, vertical depth of the fracture, volumn and mass of the deposited snow as well as the
        runout length of the avalanche event. Thereby the largest possible event along a :term:`thalweg` is classified with R5.

    Dmax
        Dmax is the size of the maximum potential avalanche event size in an infinite time series of avalanches 
        within the terrain of a given avalanche path. This maximum potential avalanche event corresponds with the 
        relative size R5 (major or maximum, relative to path) of the relative avalanche size classification 
        (American Avalanche Association, 2016). Therefore Dmax is defined accordingly to the destructive size
        (accoridng to the avalanche size classifications of the Canadian Avalanche Association and the European
        Avalanche Warning Services) of the R5 event. 
        

    velocity
        Flow velocities are usually measured in a surface parallel direction.
        Alternatively approach velocities are measured along the line of sight.

    wet snow
        A wet snow avalanche (WSA) implies the presence of liquid water within
        the avalanche and is usually associated to :term:`dense flow` type of movement
        in the :term:`transition` zone of the avalanche.

    origin point
        see also :term:`zone of origin`

        The origin point describes the highest possible release point along a :term:`thalweg` within the terrain.

    zone of origin
        The zone of origin delineates the area, in which typical :term:`release area` (s) are located, and an avalanche's
        appearance is characterized by the :term:`manner of starting`. The uppermost
        possible point is referred to as :term:`origin point`.

    transition point
        see also :term:`zone of transition`

        The transition point describes the transition between the zone of origin and the zone of transition.
        It marks the lowest point of the release area. It is generally assumed that the transition occurs at
        a slope angle of approximately 30°.

    zone of transition
        The zone of transition is the area between :term:`zone of origin` and
        :term:`zone of deposition` along the :term:`thalweg`. The :term:`form
        of movement` is linked to the :term:`flow variables`. The start of
        transition links the :term:`zone of origin` and transition and is usually
        associated with a slope inclination of about 28-30°.

    deposition point
        see also :term:`zone of deposition`

        The deposition point describes the point along the thalweg where an avalanche event starts to decelerates
        due to increased friction and terrain characteristics. The point of deposition marks the end of the travel
        track and the beginning of the deposition area. It is located between the :term:`zone of transit` and the 
        :term:`zone of deposition`. It is also known as the :term:`beta point`.

    zone of deposition
        The zone of deposition is where the :term:`runout area` of the avalanche is located and where the
        avalanche stops due to frictional energy dissipation. Whithin this zone the avalanche deposit is to be found.
        The boundary with the :term:`zone of transition` is often called :term:`beta point`.
