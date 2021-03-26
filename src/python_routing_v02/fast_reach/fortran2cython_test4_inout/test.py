import numpy as np
import diffusive

print("This is a test to make sure data passing and recieving happen correctly.")

mxncomp_g = 3
nrch_g = 6
z_ar_g = np.arange(mxncomp_g*nrch_g, dtype = np.double).reshape(mxncomp_g,nrch_g)
ntss_ev_g = 7

print("Here is a 2D array that we will pass to Fortran:")
print(z_ar_g)
print("--------------------------------------------")
print("Now, Fortran will crate a 3D array where the 2D array above is repeated across",ntss_ev_g,"planes")
print("--------------------------------------------")

q_ev_g, elv_ev_g  = diffusive.compute_diffusive(
                                mxncomp_g,
                                nrch_g,
                                np.asfortranarray(z_ar_g),
                                ntss_ev_g,
                               )

print("Here is the 3D array returned from Fortran:")
print(q_ev_g)