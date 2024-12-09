Glossary
========

.. Note:: **This is still under construction**

.. Note:: The main purpose of this glossary is to introduce and consolidate an
   unambiguous language and definitions related to avalanche dynamics,
   avalanche simulations, and modeling for the digital toolbox AvaFrame and
   in particular the Avalanche Modelling Atlas (AMA). Terms are organized in
   alphabetical order, and the terminology builds on the UNESCO Avalanche Atlas
   :cite:`DQ1981`. It also relates to existing guidelines provided by
   avalanche associations and warning services (e.g. the European Avalanche Warning
   Services (EAWS) :cite:`EAWS.2024`, the Canadian Avalanche Association (CAA)
   :cite:`CAA_OGRS2016`, and the American Avalanche Association (AAA) :cite:`AAA_SWAG2016`).
   The glossary follows common conventions used in related publications
   and projects (see our :ref:`zreferences:References`).

.. glossary::
    :sorted:

    alpha point
        The alpha point refers to the furthest point of the runout of an avalanche :term:`event`
        and marks the terminus of the avalanche deposition.

        see also :term:`runout point`

    alpha angle
        The alpha angle is the angle between the :term:`origin point` and the :term:`runout
        point`. It describes the inclination of the line connecting the top of the avalanche release and the
        furthest point of the deposition. The alpha angle is based on the :term:`runout length`
        and the altitude difference between the :term:`origin point` and the :term:`runout point` along the :term:`thalweg`.
        It is an important element in avalanche simulation and modeling (e.g. alpha-beta model :cite:`LiBa1980`,
        FlowPy simulations). The alpha angle is also referred to as the 'angle of reach' or
        the :term:`runout angle`.

        see also :term:`runout angle`

    beta point
        The boundary between the :term:`zone of transition` and :term:`zone of deposition` (start of
        deposition) is often referred to as the 'beta point' or :term:`deposition point` and is associated
        with a certain slope threshold of the :term:`thalweg`. For major avalanche :term:`path` s that may produce
        avalanches of :term:`size` 4-5, the beta point is associated with the 10° point
        (or beta_10), i.e. the point where the slope angle of the :term:`thalweg`
        decreases below 10°. For an avalanche :term:`path` with a :term:`size` 1-3 :term:`event`,
        the beta point may accordingly be associated with larger :term:`thalweg`
        slope angles (up to 30°, see :term:`zone of transition`).

        see also :term:`deposition point`

    beta angle
        The beta angle describes the inclination along the :term:`thalweg` between the
        :term:`origin point` and the :term:`deposition point`. The corresponding :term:`travel angle` is
        referred to as 'beta angle'. It is an important concept in avalanche modeling and simulation. 
        The beta angle increases with decreasing avalanche :term:`size`.

        see also :term:`travel angle`

    coordinate transformation
        Coordinate transformations refer to the operation of changing
        coordinates, e.g. between a fixed, Eularian, global coordinate system
        with geographical orientation, to an avalanche :term:`path` dependent
        coordinate system along the :term:`thalweg`, or even a Langrangian coordinate
        system, moving along particle trajectories.

    projection
        Projection refers to
        the transition of a defined set of coordinates from one coordinate system to another, i.e. the
        projection of the :term:`runout point` (as furthest reach of the avalanche)
        to the :term:`thalweg`. Or projecting a 3D travel length (xyz) to a 2D
        travel length measured only along xy.

    avalanche cycle
        An avalanche cycle describes a series of natural :term:`event` (s) that occur across a region
        over a relatively short time span (hours to days).

    danger scale
        The avalanche danger scale refers to the avalanche hazard and is an
        inherent part of avalanche warning.

    dense flow
        Dense flow is a :term:`form of movement` in the :term:`zone of transition` of the
        avalanche. Dense Flow Avalanches (DFA) flow along the ground. Mixed types of movement are often observed,
        combining different flow regimes and their partial or complete
        transitions, e.g. ‘mixed flow and powder avalanches’ or ‘flow avalanche
        with powder component’, towards the evolution of a fluidized layer in
        the avalanche flow (see :term:`powder snow`).

    depth
        Release, :term:`entrainment`, flow, or deposition depth refers to the extent of
        the avalanche measured in the direction of gravity.

        see also :term:`thickness`

    thickness
        Release, :term:`entrainment`, flow, or deposition thickness refers to the
        extent of the avalanche, measured perpendicular to the slope.

        see also :term:`depth`

    density
        Release, :term:`entrainment`, flow, or deposition density. Important quantity
        relating to mass and volume, influencing impact pressure and in particular friction relations.

    entrainment
        Entrainment describes the process of mass intake during the avalanche flow.

    event
        An avalanche associated with a certain avalanche :term:`path`
        Can be observed and documented avalanches, but also design events
        (eg. for certain :term:`return period` (s)), which are of particular interest for engineering applications.
        It has properties that are morphologically connected to
        different zones (i.e. :term:`zone of origin`, :term:`zone of transition`, :term:`zone of deposition`)
        of an avalanche :term:`path` and allow defining other associated properties, such as :term:`alpha angle`
        or :term:`runout length` that are defined in combination with the avalanche :term:`thalweg`.

    release scenario
        One or more release areas, which are associated with a certain avalanche :term:`path` and release a
        the same time.

    flow variables
        Flow variables include the flow :term:`thickness`, flow :term:`velocity`, and flow :term:`density` and are
        determined by the :term:`form of movement`. The spatio-temporal evolution of these variables
        are usually calculated by the implemented flow models. These models output
        the maximum values over the whole flow, or peak values thoughout the computational duration.
        The flow variables are used to derive other variables such as impact pressure or kinetic
        energy of the flow.

    form of movement
        Is an avalanche criterion in the :term:`zone of transition` and can have :term:`dense flow` and/or
        :term:`powder snow` characteristics.

    manner of starting
        Is an avalanche criterion in the :term:`zone of origin` and has the
        characteristics 'loose', 'slab', or 'gliding'.

    terrain classification
        Terrain may be classified according to the Avalanche Terrain Exposure
        Scale (ATES) into 'simple' (low angle or primarily forested terrain with
        some openings that may involve the :term:`zone of deposition` of infrequent avalanche
        :term:`path` (s)), 'challenging' (well defined avalanche :term:`path` (s), starting zones, or
        terrain traps), 'complex' (exposure to multiple overlapping avalanche
        :term:`path` (s), large expanses of steep, open terrain, multiple starting zones,
        and terrain traps below) and 'extreme' (exposure to very steep faces with cliffs, spines, couloirs, crevasses
        or sustained overhead hazard).

    thalweg
        The thalweg is defined as the line representing the main flow direction of all potential avalanche
        events within a specific avalanche :term:`path`. The thalweg is delineated according to the terrain 
        characteristics and is independent of a specific :term:`event`. Technically it is the two-
        dimensional terrain representation along the horizontally projected :term:`travel length` and the 
        altitude difference.

    trajectory length
        Used in com1DFA particle dictionaries, where the trajectory length is computed as the
        distance traveled by a particle from one time step to the next and then accumulated over
        time. Three different trajectory lengths are computed (1) trajectoryLengthXY - computed in the
        x, y plane, (2) trajectoryLengthXYZ - also taking the slope of the topography into account, and
        (3) trajectoryLengthXYCor - same as trajectoryLengthXY but corrected for the potential
        angle difference of the slope and the normal.

    travel angle
        The travel angle describes the inclination between the :term:`origin point` and the current location of interest. 
        It is calculated based on the :term:`travel length` and the altitude difference between the :term:`origin point` and the
        point of interest. Important travel angles are the :term:`beta angle` and the :term:`alpha angle`. 

    travel length
        Travel lengths are measured as horizontally projected travel length
        (:math:`s_{XY}`) along the :term:`thalweg`, between the current
        location with the uppermost point of the release, :term:`origin point`.
        Alternatively, the surface parallel travel length (:math:`s_{XYZ}`) may be
        defined as the three-dimensional length travelled by the avalanche.

    path
        The avalanche path summarizes the total catchment and can be divided into
        different zones (:term:`zone of origin`, :term:`zone of transition`, :term:`zone of deposition`)
        with different criteria and characteristics. An inherent property of the avalanche path is the
        :term:`thalweg` and the associated avalanche :term:`event`.

    powder snow
        Powder snow is a :term:`form of movement`,  referring to the dust or suspension
        cloud of a powder snow avalanche (PSA). PSAs are associated with cold, dry cohesionless snow. In reality, mixed
        types of movement are often observed, combining different flow regimes and their transitions, e.g. ‘powder
        avalanche with :term:`dense flow` component’.

    return period
        The return period is the average time interval between events of a certain intensity or :term:`size`.
        The return period is often determined by the :term:`runout length` of historically documented avalanche :term:`event` (s).
        It is a time-based measure (how often an event occurs).
        Return periods are often used for planning and risk management.

    return level
        The return level is the magnitude or size of an event that is expected to be exceeded with a certain probability
        over a specified time period. It is a magnitude-based measure (the size of the event).
        Return levels are used for designing infrastructure (like dams) to withstand specific events.



    release area
        Potential release areas are located in the :term:`zone of origin`. Each documented
        :term:`event` or simulation :term:`scenario` is associated to one or more primary
        and/or secondary release areas, that can further be described by the :term:`manner of starting`.

    runout area
       The runout area is associated to a specific avalanche event and usually located in the :term:`zone of deposition`.

    runout angle
        The runout angle, also referred to as :term:`alpha angle`. It describes the inclination of the
        the avalanche event from the point of release or :term:`origin point` to the :term:`runout point`.

    runout length
        Runout length is intrinsically linked to the :term:`alpha point`, utilizing the
        :term:`projection` to the :term:`thalweg`. In the same manner as :term:`travel length`, runout lengths
        are measured as horizontally projected lengths along the :term:`thalweg`, from
        the uppermost point of the :term:`release area` (:term:`origin point`) to furthest reach of the runout
        area (:term:`runout point`).

    runout point
        The runout point is also referred to as the :term:`alpha point`. It describes the furthest point of the runout,
        and marks the outer most end of the avalanche deposit. The runout may refer to visible deposition
        (associated to :term:`dense flow`), damages or the impacted and affected area (associated to air
        blast or :term:`powder snow`) in the :term:`zone of deposition`. It is usually defined
        via flow :term:`thickness`, velocity, kinetic energy or impact pressure thresholds.

    size
        Two prevalent avalanche size classifications exist.

        see :term:`destructive size` and :term:`relative size`

    destructive size
        The destructive size refers to the magnitude or intensity of an avalanche :term:`event`, according to the EAWS size
        classification (:cite:`EAWS.2024`), which is closely related to the CAA destructive size (:cite:`CAA_OGRS2016`).
        Thereby the size refers to the destructive potential, :term:`runout length`, volume or mass of each avalanche event.

    relative size
        The relative size refers to the size of an avalanche in relation to the :term:`thalweg`. The classification
        scheme was developed by the AAA (:cite:`AAA_SWAG2016`) and relies on the horizontal extent, vertical depth of
        the fracture, volume and mass of the deposited snow as well as the runout length of the avalanche event.
        The largest possible event along a :term:`thalweg` is classified with R5.

    Dmax
        Dmax is the size of the maximum potential avalanche event size in an infinite time series of avalanches 
        within the terrain of a given avalanche path. This maximum potential avalanche event corresponds to the 
        relative size R5 (major or maximum, relative to path) of the relative avalanche size classification 
        (:cite:`AAA_SWAG2016`). Therefore Dmax is defined according to the destructive size
        (:cite:`CAA_OGRS2016`, :cite:`EAWS.2024`) of the R5 event. 

    velocity
        Flow velocities are usually measured in a surface parallel direction. Alternatively, approach velocities are
        measured along the line-of-sight.

    wet snow
        The term wet snow avalanche (WSA) implies the presence of liquid water within
        an avalanche and is usually associated to :term:`dense flow` type of movement
        of the avalanches :term:`zone of transition`.

    origin point
        The origin point refers to the highest possible release point along a :term:`thalweg` within the terrain.

        see also :term:`zone of origin`

    zone of origin
        The zone of origin delineates the area, in which typical :term:`release area` (s) are located, and an avalanche's
        appearance is characterized by the :term:`manner of starting`. The uppermost
        possible point is referred to as :term:`origin point`.

    transition point
        The transition point describes the transition between the :term:`zone of origin` and the :term: 'zone of transition'.
        It marks the lowest point of the release area. It is often assumed that the transition occurs at
        a slope angle of approximately 30°.

        see also :term:`zone of transition`

    zone of transition
        The zone of transition is the area between the :term:`zone of origin` and the
        :term:`zone of deposition` along the :term:`thalweg`. The :term:`form
        of movement` is linked to the :term:`flow variables`. The start of the
        zone of transition is usually associated with a slope inclination of about 30°.

    deposition point
        The deposition point describes the point along the :term:`thalweg`, where an avalanche event starts to decelerate,
        due to increased friction and terrain characteristics. The point of deposition marks the end of the
        :term:`zone of transition` and the beginning of the :term:`zone of deposition`. It is also known as the :term:`beta point`.

        see also :term:`zone of deposition`

    zone of deposition
        The zone of deposition is where the :term:`runout area` of the avalanche is located and the
        avalanche stops due to frictional energy dissipation. The majority of the final avalanche deposit is located within this zone.
        The boundary with the :term:`zone of transition` is often referred to as the :term:`beta point`.
