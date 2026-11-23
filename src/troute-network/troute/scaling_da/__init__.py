"""Area-scaled streamflow data assimilation (off by default).

What it does
------------
Streamflow nudging corrects the flow at each observed gage and lets that
correction travel downstream. Area-scaled DA also corrects the river network
upstream of the gage, where nudging changes nothing.

After each routing loop, the correction that nudging applied at a gage
(dQ_gage, one value per timestep) is spread over the gage's upstream tree:

    along a channel:    dQ_up = dQ_gage * (A_up / A_gage) ** theta
    at a confluence:    each branch gets the parent's correction times
                        Q_branch / max(Q_parent, sum of Q over all branches)

A is the total drainage area of a segment (hydrofabric total_dasqkm) and Q its
modeled flow. theta defaults to 0.77 and can be set per VPU or per gage. The
multiplier is clipped to [0, 1], so a correction never grows away from the gage.

A tree stops at the next observed gage, at lakes and reservoirs, at the
diffusive domain, and max_reach_km (default 200 km) upstream along the network,
so no segment is corrected by two gages. Corrected flows are floored at 0 and
depth follows the flow (h ~ Q ** 0.6, ratio clipped to [1/2, 2]); velocity is
not changed.

The corrected flow is written to the loop's output every loop, but it enters
the routing state only after the final loop (the hand-off to a following
forecast). It is enabled by streamflow_da.streamflow_scaling, needs
streamflow_nudging, and runs only with HYFeatures networks in the -V4 driver.

Modules: trees (build the trees), spread (the multiplier arithmetic),
runtime (setup and the per-loop step called from main_v04).

Limitations
-----------
- theta = 0.77 comes from flood-peak scaling in a small (0.17-21.2 km2)
  Hortonian basin; using it for per-timestep corrections in larger basins
  goes beyond that study.
- The confluence flow split and the distance cap are heuristics.
- The correction is an analysis increment, not routed water: mass is not
  conserved (the injected volume is logged each loop).

Reference: Ogden, F. L., and Dawdy, D. R. (2003). Peak discharge scaling in
small Hortonian watershed. J. Hydrol. Eng. 8(2):64-73. The method follows the
NGWPC t-route implementation of area-scaled streamflow DA.
"""
