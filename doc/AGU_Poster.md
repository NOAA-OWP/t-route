# Efficient Routing Computations with a Graph-Based Routing Framework

<img src=https://s3.amazonaws.com/ipostersessions-agu/3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37.png width="200" height="100">

__James S. Halgren, Dong Ha Kim, Juzer F. Dhondia, Nels J. Frazier,
Nicholas Chadwick, Ryan D. Grout, Alexander A. Maestre, Jacob E. Hreha, Adam N. Wlostowski, Graeme R.
Aggett__

_Office of Water Prediction, NOAA/NWS, Tuscaloosa, AL, USA;_<br>
_Lynker Technologies, Leesburg, VA, USA;_<br>
_ERT, Inc., Laurel, MD, USA;_<br>
_University of Alabama Remote Sensing Center, Tuscaloosa, AL, USA_<br>

\
<img src=https://s3.amazonaws.com/ipostersessions-agu/6f1e8bbb-6ff0-459d-8f73-25ee1cc5f134.jpg alt="@jameshalgren" width="100" height="100"> 
<img src=https://s3.amazonaws.com/ipostersessions-agu/6e0a282b-3631-4874-9533-3a9e58fc4e12.png alt="@kumdonoaa" width="100" height="100"> 
<img src=https://s3.amazonaws.com/ipostersessions-agu/2d5dd5fa-5ea3-4e0e-ae3d-bad1f075608c.jpg alt="@jdhondia-noaa" width="100" height="100"> 
<img src=https://s3.amazonaws.com/ipostersessions-agu/d22d5f15-69ca-48ea-b5df-420ec3e2eb60.jpeg alt="@groutr" width="100" height="100"> 
<img src=https://s3.amazonaws.com/ipostersessions-agu/6358cd22-3645-46ca-a31b-1ea2cc498f80.png alt="@awlostowski-noaa" width="100" height="100"> 
<img src=https://s3.amazonaws.com/ipostersessions-agu/ede3d422-0af6-4ea4-baf1-44eb027d75ae.jpg alt="gaggett" width="100" height="100"> 

__Originally a Poster Presented at:__
![](https://agu2020fallmeeting-agu.ipostersessions.com/GetFile.ashx?file=020_40144_FM_LogInPage_Header_1920x240-01.jpg)


# The CONUS routing challenge
Under direction of the NOAA-NWS Office of Water Prediction (OWP) we have
created **a new routing framework for the National Water Model
(NWM).** This new framework permits use of advanced routing methods but
implies additional compute burden. Our challenge is to implement a parallel execution 
method to permits efficient calculation of the advanced routing while still satisfying 
the operational execution constrations for use with the National Water Model.

## The National Water Model (NWM) Channel Routing Network
Based on the NHD+ V2.0 Medium Resolution data set, the National Water Model CONUS river network is composed of:

-   4.3M stream/river miles
-   2,729,077 individual segments
-   2,102,010 reaches 
-   1,029,217 junctions, and 
-   14,713 independent drainage basins (disjoint networks). 

### Key Challenges for CONUS-scale routing
#### Computation Volume
Continental scale routing in the NWM is an enormous computational challenge. The NWM routing computation includes:
-   multiple forecast realizations (analysis, short, medium-range
    deterministic, medium-range ensemble, and long-range ensemble),
-   carried out on over 4.3M river miles consisting of 2.7M+ segments
    ...
-   representing over 600 billion routing computations daily, 
-   or about 7 million routing calculations per second on average.

*Table: Daily Volume of Operational NWM Routing Calculations*

<!--
_NWM Forecast_ | _Time Steps_ | _Segments (approx.)_ | _Number of Times per day_ | _Number of Ensemble Members_ | _Total Daily Routing Calculations_
-- | -- | -- | -- | -- | --
Analysis and Assimilation | 3 | 2,700,000 | 24 | 1 | 194,400,000
Short Range | 216 | 2,700,000 | 24 | 1 | 13,996,800,000
Medium Range deterministic | 2880 | 2,700,000 | 4 | 1 | 31,104,000,000
Medium Range ensembles | 2448 | 2,700,000 | 4 | 7 | 185,068,800,000
Long Range | 8640 | 2,700,000 | 4 | 4 | 373,248,000,000
  |   |   |   | TOTAL | 603,612,000,000
-->

<img src=https://s3.amazonaws.com/ipostersessions-agu/536eeae6-5b42-45ae-8f72-ff80fd3335bd.png height="200">

#### Topological dependencies
OWP is developing the new framework to support diffusive- and
dynamic-wave hydraulic routing simulations. These higher-order 
routing methods require explicit topological tracking throughout 
network. 

Tracking the topological connectivity and enforcing dependence of the
calculations permits use of these routing methods, but also adds an
additional challenge to NWM routing. Instead of being
completely embarrassingly parallel within each timestep, some calculations
must wait for others to be completed. (The current
method uses simplifying assumptions, incompatible with higher-order
routing solutions, to allow for simplified, fully parallel routing
execution within each timestep.)
Our challenge was to introduce this topological dependency in the NWM
routing framework while still managing the required hourly calculation
volume. 
Here, we present our new routing framework and a graph-based computing scheme 
to drive rapid parallel computation while maintaining the topological dependence
of the routing network. 

## Graph Representation of NWM Channels
We represent the CONUS river network as a series of directed acyclic
graphs, each consisting of a hydraulically independent drainage basin
exiting to the ocean or to an inland sink.

<img src=https://s3.amazonaws.com/ipostersessions-agu/d3725fe6-d7ae-40cb-aae2-93f3b29da588.png width="200">

*Caption: Elementary components of the CONUS river network graph. The
smallest elements, denoted by discrete colors, are individual stream
segments. Linear combinations of segments between junctions form
reaches. Junctions exist at the confluence of two or more reaches.*

Network complexity, expressed as a number of junctions, is a useful
measure of the level of dependence of the graph, and gives an idea of
the computational burden for each independent network.

<img src=https://s3.amazonaws.com/ipostersessions-agu/dfdf0e2c-b9f3-4746-978d-347d2e457e06.png height="400">

*Caption: Size distribution of river networks with greater than 2k junctions in the CONUS NWM dataset.*

### Independent Networks in the CONUS Dataset
<img src=https://s3.amazonaws.com/ipostersessions-agu/658cf3ec-95dc-4fb5-a4c5-528c69917e67.png height=200>

*Caption: Distribution of the 213 independent networks in the CONUS NWM
dataset of river segments below all existing national weather service
forecast points called the “Mainstem” network. Colors are used to
indicate independent networks, each draining to unique tailwater
segments at the CONUS boundary.*

Terminal segments, which represent river outlets emptying into the ocean
or an inland sink, define independent river networks within the larger
NWM routing dataset. There are more than 5,000 terminal segments of
order 1, meaning that they define an independent network consisting of a
single reach with no tributaries emptying directly into the ocean. Only
the Mississippi outlet reaches order 10

<img src=https://s3.amazonaws.com/ipostersessions-agu/90e08f09-c4d5-4ef9-bc1d-21185f9d2750.png height=300>

*Caption: Tabulation of orders represented in the 14k+ independent
networks of the NWM channel dataset network. Distribution of the 213
independent networks in the CONUS NWM dataset of river segments below
all existing national weather service forecast points.*

## Parallel Methods: Network-based or Just-in-Time?

Two parallelization approaches, By-network, and Just-in-time (JIT)
as detailed below, were tested in comparison to pure serial computation:

<img src=https://s3.amazonaws.com/ipostersessions-agu/1ae5ea14-48a9-4ee3-a6fd-7762207afa81.gif height=150>

*Caption: Animation of Just-In-Time network traversal, with calculations
on separate portions of the tree coalescing to finish simultaneously at
the outlet. At any moment in the animation, red sparks highlight reaches
of a common reverse network order that may be computed in parallel.*

Serial computation 
-   Starting upstream, proceeding downstream,
-   One network at a time.
-   Computationally inefficient, but a useful benchmark.

Independent network parallelization
-   Divide computation according to separate networks, e.g., the
    Colorado River, Mississippi, etc. are computed independently. 
-   Performance limited by the size of the largest basin, i.e., the
    Mississippi. 

Just-In-Time parallelization

-   Calculate first the headwater reaches of edges of the longest
    network...
-   Followed by all reaches below headwaters, etc.
-   Orchestrating the computation so each dependency is computed just
    before it is needed downstream (hence "Just-in-time") provides best
    theoretical potential speedup.
-   Practical efficiencies are obtained by grouping the reaches into
    cascading orders of subnetworks.

<img src=https://s3.amazonaws.com/ipostersessions-agu/110c901d-0bc5-499b-bcbf-75a2893255c6.png height=400>

*Caption: Grouping the reaches into subnetworks balances the practical
impact of many parallel calls. An optimal subnetwork size is small
enough to permit sufficient parallelism and large enough to ensure that
parallel overhead is not burdensome. Colors denote subnetworks of common
reverse network order. Higher order subnetworks are computed prior to
lower order subnetworks in order to maintain topological dependencies.*

## Parallel scaling results

By representing the NWM routing network as a graph, **we have achieved up
to 20% of the theoretical 40x speedup** possible with carefully
orchestrated graph-based parallelization.

**We estimate an approximately 1000x theoretical potential speedup for
the full resolution NWM routing dataset.**

We estimated speedup with additional CPU cores for four cases: 

-   theoretical network-based, 
-   theoretical by reaches (which is the ideal JIT, assuming no parallel
    overhead), 
-   real network-based, 
-   real Just-in-time (JIT).

Tests were conducted primarily on the subset of channels below existing
National Weather Service forecast points, referred to as the
“Mainstems”.

<img src=https://s3.amazonaws.com/ipostersessions-agu/fac0a954-5cf6-48e7-baed-8a8d97bfaa84.png height=200>

*Caption: Performance improvement with additional parallel cores for
Mainstems network domain. Note that the network-based performance is
very near the theoretical maximum; the JIT performance is significantly
better, but falls short of the theoretical maximum.*
  

Theoretical potential parallel computational speedup is calculated
as the ratio of the total segment count to the parallel method limiting
size. 

Network-based parallel execution limited by total segment count of
largest network. 

-   Estimated 2x speedup in both Mainstems and Full Resolution datasets
    (Mississippi basin accounts for roughly half of the total segment
    count.)
-   Experimental results yielded close to the theoretical maximum. 

JIT method speedup limited by largest network depth, i.e., length of
path from the furthest headwater to the outlet. 

-   Depth of Mainstems network mississippi basin is 73 reaches compared
    to 27k reaches overall so JIT theoretical maximum improvement is 27k
    / 73 (i.e., \~40x) 
-   Experimental results yield about 10% efficiency (i.e. \~4x speedup).
    Careful re-distribution of subnetworks to threads has yielded a 20%
    efficiency in a smoke test (i.e. 8x speedup).
-   Theoretical potential performance improvement grows with network
    size -- the depth of the Mississippi Basin in the Full resolution
    dataset is 2218 reaches compared to 2.1M reaches overall so JIT
    theoretical maximum improvement is 2.1M / 2218 (i.e., \~1000x).

The actual speedup will be affected by numerous computational realities
including: i/o overhead, parallel thread or process pool spin-up time,
array access efficiency (i.e. cache misses), etc.

---
# Learn more: T-route on GitHub

The new routing framework is publicly developed and we encourage
interested community members to access...

**[http://github.com/NOAA-OWP/t-route](http://github.com/NOAA-OWP/t-route)**

...to try the approach, provide feedback, and contribute to further
development. Our goal is to significantly reduce barriers to efficient
application of higher order routing solutions in the National Water
Model, enabling more useful forecasts that help communities prepare
effectively for hydrologic hazards.


-----
# Author Information

The authors are part of the development team for the National Water
Model at the Office of Water Prediction for NOAA's National Weather
Service. We wish to gratefully acknowledge the excellent work of the OWP
and NCAR teams (and others) who have prepared the National Water Model
as it stands today. Author affiliations are as follows: James S.
Halgren<sup>1,2</sup>, Dong Ha Kim<sup>1,2</sup>, Juzer F. Dhondia<sup>1,4</sup>, Nels J.
Frazier<sup>1,3</sup>, Nicholas Chadwick,<sup>1,3</sup>, Ryan D. Grout<sup>1,2</sup>, Alexander A.
Maestre<sup>1,4</sup>, Jacob E. Hreha<sup>1,2</sup>, Adam N. Wlostowski<sup>1,2</sup>, Graeme R.
Aggett<sup>2 </sup>

<sup>1</sup>Office of Water Prediction, NOAA/NWS, Tuscaloosa, AL, USA <br>
<sup>2</sup>Lynker Technologies, Leesburg, VA, USA <br>
<sup>3</sup>ERT, Inc., Laurel, MD, USA <br>
<sup>4</sup>University of Alabama Remote Sensing Center, Tuscaloosa, AL, USA <br>

---
# Abstract
To resolve non-uniform and unsteady flows in the National Water Model
(NWM), the Office of Water Prediction is developing additional routing
engines to power simulations with the dynamic and diffusive
approximations of the St Venant equations. This gives rise to two major
computational challenges. First, the presence of both upstream and
downstream boundary conditions requires tracking topological
connectivity of the entire network within the computation. Second, all
solution methods, whether explicit or implicit, become computationally
expensive when scaled to continental domains. To be viable for
operational modeling as an element of the National Water Model, the
computational framework for dynamic routing must address these
challenges.

We present a continental-scale flow routing framework that represents
the flow network as a collection of directed acyclic graphs where edges
point in the direction of downstream flow. We use information from this
graph representation to efficiently drive a parallelized computation of
flow from headwaters downstream to the tailwaters. This approach has
achieved modest performance gains in terms of overall compute time and
resources for the routing cases we have tested. The framework is
publicly developed and we encourage interested community members to use
our approach and provide feedback.

Initial results show that we can simulate 5 days of continental scale
flow routing below all existing national weather service forecast points
in approximately 10 minutes using only 4 processors. Also, the new
framework permits computation using upstream dependencies in all
timesteps, which is not possible in the present NWM routing framework.
We will continue our work with the goal of significantly reducing
barriers to efficient application of higher order routing solutions in
the National Water Model, enabling more useful forecasts that help
communities prepare for hydrologic hazards.

[![](./AGU%20-%20iPosterSessions.com_files/Paper_758298_abstract_728584_0.gif)](./AGU%20-%20iPosterSessions.com_files/Paper_758298_abstract_728584_0.gif)

To resolve non-uniform and unsteady flows in the National Water Model
(NWM), the Office of Water Prediction is developing additional routing
engines to power simulations with the dynamic and diffusive
approximations of the St Venant equations. This gives rise to two major
computational challenges. First, the presence of both upstream and
downstream boundary conditions requires tracking topological
connectivity of the entire network within the computation. Second, all
solution methods, whether explicit or implicit, become computationally
expensive when scaled to continental domains. To be viable for
operational modeling as an element of the National Water Model, the
computational framework for dynamic routing must address these
challenges.

We present a continental-scale flow routing framework that represents
the flow network as a collection of directed acyclic graphs where edges
point in the direction of downstream flow. We use information from this
graph representation to efficiently drive a parallelized computation of
flow from headwaters downstream to the tailwaters. This approach has
achieved modest performance gains in terms of overall compute time and
resources for the routing cases we have tested. The framework is
publicly developed and we encourage interested community members to use
our approach and provide feedback.

Initial results show that we can simulate 5 days of continental scale
flow routing below all existing national weather service forecast points
in approximately 10 minutes using only 4 processors. Also, the new
framework permits computation using upstream dependencies in all
timesteps, which is not possible in the present NWM routing framework.
We will continue our work with the goal of significantly reducing
barriers to efficient application of higher order routing solutions in
the National Water Model, enabling more useful forecasts that help
communities prepare for hydrologic hazards.

![](https://agu.confex.com/data/abstract/agu/fm20/8/9/Paper_758298_abstract_728584_0.gif)
