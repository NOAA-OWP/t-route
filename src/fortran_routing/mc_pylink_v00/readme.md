# Muskingum Cunge Methods
## Method 1
* Contain all loops (segment, timestep) in Python, Call Only a single segment Fortran computation
* qd_curr, vel_curr, depth_curr = f(qd_prev, qup_prev, qup_curr, bw, tw, s0, n, twcc, ncc, cs, qlat) 

## Method 2
* Timestep loop python, Segment loop (by reach) fortran 
### Method 2a
* The current method

### Method 2b
* Remove the 'iseg' from the mcNWM call

## Method 3
* Multi-segment fortran

## Method 4
* Multi-timestep, multi-segment fortran

# TODO:
* Solidify test dataset and canonical output (Dong Ha)

# MSH-based dynamic Methods
## MSH rectangular

## MSH compound channel 

## MSH KDD

