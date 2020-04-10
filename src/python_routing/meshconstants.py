PHI = 1.0
THETA = 1.0
THETAS = THETA
# thetas is the weighting coefficient for the source term
# (to determine if we want to compute the source term in an
# implicit or an explicit manner.)
THESINV = 1.0
ALFA2 = 0.5
ALFA4 = 0.1
CELERITY_EPSILON = 0.00001
DX_TOLERANCE = 0.00001
DEPTH_TOLERANCE = 0.00001
AREA_TOLERANCE = 0.0001
'''
1.0 =: phi           temporal weighting coefficient
1.0 =: theta         spatial weighting coefficient for all terms except the source term
1.0 =: thetas        spatial weighting coefficient for the source term (0:explicit, 1:implicit)
1.0 =: thesinv       ?
0.5 =: alfa2         emp parameter for artificial diffusion (lit)
0.1 =: alfa4         maximum value for artificial diffusion
'''
