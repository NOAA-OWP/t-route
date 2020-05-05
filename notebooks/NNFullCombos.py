#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from pylab import *
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense,Dropout
import numpy as np
import json
import time
import sys;sys.path.append(r'../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS')
import mc_sseg_stime_NOLOOP as mc
import numba
from numba import jit

depthp_min = 0.010; depthp_max = .011
qlat_min = 35; qlat_max = 45
qdp_min = .01; qdp_max = 1
quc_min = .01; quc_max = 1
qup_min = .01; qup_max = 1
s0_min = .00001; s0_max = .002; 
cs_min = .085; cs_max = 2.254;
tw_min = 150; tw_max = 500;
bw_min = 112; bw_max = 150;

# singlesegment():
array_length = 6   

dt = 60 # Time step
dx = 1800 # segment length
# bw = np.linspace(0.135, 230.035, array_length, endpoint=True) # Trapezoidal bottom width
bw = np.linspace(bw_min, bw_max, array_length, endpoint=True) # Trapezoidal bottom width
tw = np.linspace(tw_min, tw_max, array_length, endpoint=True) # Channel top width (at bankfull)
# twcc = np.linspace(0.67, 1150.17, array_length, endpoint=True) # Flood plain width
# twcc = tw*  # Flood plain width tw*3
n_manning = .028   # manning roughness of channel
n_manning_cc = .028 # manning roughness of floodplain
cs = np.linspace(cs_min,cs_max, array_length, endpoint=True)# channel trapezoidal sideslope
s0 = np.linspace(s0_min, s0_max, array_length, endpoint=True) # Lateral inflow in this time step
qup = np.linspace(qup_min, qup_max, array_length, endpoint=True) # Flow from the upstream neighbor in the previous timestep
# quc = np.linspace(10, 1000, array_length, endpoint=True) # Flow from the upstream neighbor in the current timestep 
quc = np.linspace(quc_min, quc_max, array_length, endpoint=True)  # Flow from the upstream neighbor in the current timestep 
# qdp = np.linspace(10, 1000, array_length, endpoint=True) # Flow at the current segment in the previous timestep
qdp = np.linspace(qdp_min, qdp_max, array_length, endpoint=True)  # Flow at the current segment in the previous timestep
qlat = np.linspace(qlat_min, qlat_max, array_length, endpoint=True) # lat inflow into current segment in the current timestep
velp = .5  # Velocity in the current segment in the previous timestep NOT USED AS AN INPUT!!!
depthp = np.linspace(depthp_min ,depthp_max , array_length, endpoint=True) # D


# @jit
def normalize(val,max,min,target_max=1,target_min=0):
    return (val-min)/(max-min)*(target_max-target_min)+target_min

def singlesegment(
        dt # dt
        , qup = None # qup can be anything
        , quc = None # quc will not be more than a 10 percent diff than qup
        , qlat = None # ql can be anything - key 
        , qdp = None # qdp will not be more than 20 percent diff than qup+qlat
        , dx = None # dx fully variable 
        , bw = None # bw correlated to tw, tw always > bw
        , tw = None # tw correlated to bw, bw always < tw
        , twcc = None # twcc always > than tw, tw of broader floodplain
        , n_manning = None # usually smaller than n_manning_cc
        , n_manning_cc = None # ncc usually greater than n_manning
        , cs = None # cs correlated to bw and tw
        , s0 = None # s0 variable 
        , velp = None # velocity at previous time step not rel
        , depthp = None # depth at previous time step starting point for iteration depthp = approx(y_direct(bw,n_manning,s0,avg(qup,qdp)))
    ):

    
    

    
    # call Fortran routine
    return mc.muskingcungenwm(
        dt, qup, quc, qdp, qlat, dx, bw, tw, twcc
        ,n_manning, n_manning_cc, cs, s0, velp, depthp
    )
    #return qdc, vel, depth
    
    

Y = []
M = []

for qup_i in range(0,(array_length)-1):
    for quc_j in range(0,(array_length)-1):
        for qlat_k in range(0,(array_length)-1):
            for qdp_l in range(0,(array_length)-1):
                for bw_n in range(0,(array_length)-1):
                    for tw_o in range(0,(array_length)-1):
                        for cs_p in range(0,(array_length)-1):
                            for s0_q in range(0,(array_length)-1):
                                for depthp_s in range(0,(array_length)-1):
                                    # M.append([dt, qup[qup_i], quc[quc_j], qlat[qlat_k],qdp[qdp_l],dx,  bw[bw_n], tw[tw_o], tw[tw_o]*3,n_manning, n_manning_cc, cs[cs_p], s0[s0_q], velp, depthp[depthp_s]])
                                    M.append([
                                        # dt, 
                                        normalize(qup[qup_i],qup_max,qup_min), 
                                        normalize(quc[quc_j],quc_max,quc_min), 
                                        normalize(qlat[qlat_k],qlat_max,qlat_min),
                                        normalize(qdp[qdp_l],qdp_max,qdp_min),
                                        # dx,  
                                        normalize(bw[bw_n],bw_max,bw_min),
                                        normalize(tw[tw_o],tw_max,tw_min),
                                        # normalize(tw[tw_o]*3,tw_max,tw_min),
                                        # n_manning, 
                                        # n_manning_cc, 
                                        normalize(cs[cs_p],cs_max,cs_min),
                                        normalize(s0[s0_q], s0_max, s0_min),
                                        # velp, 
                                        normalize(depthp[depthp_s],depthp_max,depthp_min)])

                                    S = singlesegment(
                                    dt=dt,
                                    qup=qup[qup_i],
                                    quc=quc[quc_j],
                                    qlat=qlat[qlat_k],
                                    qdp=qdp[qdp_l],
                                    
                                    dx=dx ,
                                    bw=bw[bw_n],
                                    tw=tw[tw_o],
                                    twcc=tw[tw_o]*3,
                                    n_manning=n_manning,
                                    n_manning_cc=n_manning_cc,
                                    cs=cs[cs_p],
                                    s0=s0[s0_q],
                                    velp=velp,
                                    depthp=depthp[depthp_s])
                                    Y.append((S[0]))
                                    if len(Y)%100000 == 0:
                                        print(len(Y))

dt = 60.0 # diffxxxxx
dx = 1800.0 # diffxxxxx
bw = 112.0 #small diffxxxxx
tw = 448.0# small diffxxxxx
twcc = 623.5999755859375 # no differencexxxxx
n_manning = .02800000086426735 #diffxxxxxxx
n_manning_cc = .03136000037193298 # no differencexxxxxxx
cs = 1.399999976158142 # tiny diffxxxxx
s0 = .0017999999690800905 # big diffxxxxxxxxx
qlat = 40.0 # diffxxxx
qup = .04598825052380562 # almost 1 to 1 with qucxxxx
quc = .04598825052380562#xxxxxx
qdp = .21487340331077576 # same as qup qucxxxxx
velp = .070480190217494964 # no differencedepthp = 0.010033470578491688) # large diff
depthp = 0.010033470578491688

M.append([ normalize(qup,qup_max,qup_min), 
    normalize(quc,quc_max,quc_min), 
    normalize(qlat,qlat_max,qlat_min),
    normalize(qdp,qdp_max,qdp_min),
    # dx,  
    normalize(bw,bw_max,bw_min),
    normalize(tw,tw_max,tw_min),
    # normalize(tw[tw_o]*3,tw_max,tw_min),
    # n_manning, 
    # n_manning_cc, 
    normalize(cs,cs_max,cs_min),
    normalize(s0, s0_max, s0_min),
    # velp, 
    normalize(depthp,depthp_max,depthp_min)])
Y.append(0.7570106983184814)

# M = np.array(M)
# for i in range(0,len(M),1):
#     S = singlesegment(*M[i])
#     Y.append(S[0])
Y = np.array(Y)    
M = np.array(M)

print(Y[-1])
print(M[-1])


num_samp = 10000

dt = 60 # Time step
dx = 1800 # segment length
# bw = np.linspace(0.135, 230.035, array_length, endpoint=True) # Trapezoidal bottom width
bw = np.random.uniform(bw_min, bw_max, num_samp) # Trapezoidal bottom width
tw = np.random.uniform(tw_min, tw_max, num_samp) # Channel top width (at bankfull)
# twcc = np.linspace(0.67, 1150.17, array_length, endpoint=True) # Flood plain width
# twcc = tw*  # Flood plain width tw*3
n_manning = .028   # manning roughness of channel
n_manning_cc = .028 # manning roughness of floodplain
cs = np.random.uniform(cs_min,cs_max, num_samp)# channel trapezoidal sideslope
s0 = np.random.uniform(s0_min, s0_max, num_samp) # Lateral inflow in this time step
qup = np.random.uniform(qup_min, qup_max, num_samp) # Flow from the upstream neighbor in the previous timestep
# quc = np.linsprandom.uniformace(10, 1000, array_length, endpoint=True) # Flow from the upstream neighbor in the current timestep 
quc = np.random.uniform(quc_min, quc_max, num_samp)  # Flow from the upstream neighbor in the current timestep 
# qdp = np.linspace(10, 1000, array_length, endpoint=True) # Flow at the current segment in the previous timestep
qdp = np.random.uniform(qdp_min, qdp_max, num_samp)  # Flow at the current segment in the previous timestep
qlat = np.random.uniform(qlat_min, qlat_max, num_samp) # lat inflow into current segment in the current timestep
velp = .5  # Velocity in the current segment in the previous timestep NOT USED AS AN INPUT!!!
depthp = np.random.uniform(depthp_min ,depthp_max , num_samp) # D

VAL_x = []
VAL_y = []
for i in range(num_samp):
    VAL_x.append( [normalize(qup[i],qup_max,qup_min), 
    normalize(quc[i],quc_max,quc_min), 
    normalize(qlat[i],qlat_max,qlat_min),
    normalize(qdp[i],qdp_max,qdp_min),
    # dx,  
    normalize(bw[i],bw_max,bw_min),
    normalize(tw[i],tw_max,tw_min),
    # normalize(tw[tw_o]*3,tw_max,tw_min),
    # n_manning, 
    # n_manning_cc, 
    normalize(cs[i],cs_max,cs_min),
    normalize(s0[i], s0_max, s0_min),
    # velp, 
    normalize(depthp[i],depthp_max,depthp_min)])
    S = singlesegment(
                                dt=dt,
                                qup=qup[i],
                                quc=quc[i],
                                qlat=qlat[i],
                                qdp=qdp[i],
                                
                                dx=dx ,
                                bw=bw[i],
                                tw=tw[i],
                                twcc=tw[i]*3,
                                n_manning=n_manning,
                                n_manning_cc=n_manning_cc,
                                cs=cs[i],
                                s0=s0[i],
                                velp=velp,
                                depthp=depthp[i])
    VAL_y.append(S[0])
VAL_x = np.array(VAL_x)
VAL_y = np.array(VAL_y)

import matplotlib as plt
from pylab import *
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense,Dropout,BatchNormalization
# import time
# start_time = time.time()

print(M[:100])
#Define the model
def baseline_model():
    model = tf.keras.Sequential()
    model.add(Dense(1028, activation=tf.nn.relu, input_shape=[9],use_bias=False))
#     model.add(BatchNormalization())
    model.add(Dense(512, activation='relu'))
#     model.add(Dense(256, activation='relu'))
#     model.add(Dense(128, activation='relu'))
#     model.add(Dense(64, activation='relu'))
#     model.add(Dense(32, activation='relu'))
#     model.add(Dense(16, activation='relu'))
    
    model.add(Dense(1, activation = 'linear'))

    model.compile(optimizer = 'adam', loss = 'mean_squared_error', metrics = ['mse','mean_absolute_error'])
    return model
# model.summary()
#Use the model
regr = baseline_model()

history = regr.fit(M[:], Y[:], epochs=10, batch_size=1000,  validation_data=(VAL_x,VAL_y))
# plt.plot(history.history['mse'])
# plt.plot(history.history['val_mse'])

# plt.title('Model accuracy')
# plt.ylabel('MSE')
# plt.xlabel('Epoch')
# plt.legend(['MSE', 'val_mse'], loc='upper left')
# plt.show()
# print("--- %s seconds ---" % (time.time() - start_time))
print(regr.predict(M[-1:]))
print(Y[-1])



import random 
from random import randint



Y = []
M = []

num_samp = 100

dt = 60 # Time step
dx = 1800 # segment length
# bw = np.linspace(0.135, 230.035, array_length, endpoint=True) # Trapezoidal bottom width
bw = np.random.uniform(bw_min, bw_max, num_samp) # Trapezoidal bottom width
tw = np.random.uniform(tw_min, tw_max, num_samp) # Channel top width (at bankfull)
# twcc = np.linspace(0.67, 1150.17, array_length, endpoint=True) # Flood plain width
# twcc = tw*  # Flood plain width tw*3
n_manning = .028   # manning roughness of channel
n_manning_cc = .028 # manning roughness of floodplain
cs = np.random.uniform(cs_min,cs_max, num_samp)# channel trapezoidal sideslope
s0 = np.random.uniform(s0_min, s0_max, num_samp) # Lateral inflow in this time step
qup = np.random.uniform(qup_min, qup_max, num_samp) # Flow from the upstream neighbor in the previous timestep
# quc = np.linsprandom.uniformace(10, 1000, array_length, endpoint=True) # Flow from the upstream neighbor in the current timestep 
quc = np.random.uniform(quc_min, quc_max, num_samp)  # Flow from the upstream neighbor in the current timestep 
# qdp = np.linspace(10, 1000, array_length, endpoint=True) # Flow at the current segment in the previous timestep
qdp = np.random.uniform(qdp_min, qdp_max, num_samp)  # Flow at the current segment in the previous timestep
qlat = np.random.uniform(qlat_min, qlat_max, num_samp) # lat inflow into current segment in the current timestep
velp = .5  # Velocity in the current segment in the previous timestep NOT USED AS AN INPUT!!!
depthp = np.random.uniform(depthp_min ,depthp_max , num_samp) # D


errors = np.array([])
for i in range(num_samp):
    temp_y = singlesegment(
                                dt=dt,
                                qup=qup[i],
                                quc=quc[i],
                                qlat=qlat[i],
                                qdp=qdp[i],
                                
                                dx=dx ,
                                bw=bw[i],
                                tw=tw[i],
                                twcc=tw[i]*3,
                                n_manning=n_manning,
                                n_manning_cc=n_manning_cc,
                                cs=cs[i],
                                s0=s0[i],
                                velp=velp,
                                depthp=depthp[i])

    temp_y_interp = regr.predict(np.array( [[normalize(qup[i],qup_max,qup_min), 
                                    normalize(quc[i],quc_max,quc_min), 
                                    normalize(qlat[i],qlat_max,qlat_min),
                                    normalize(qdp[i],qdp_max,qdp_min),
                                    # dx,  
                                    normalize(bw[i],bw_max,bw_min),
                                    normalize(tw[i],tw_max,tw_min),
                                    # normalize(tw[tw_o]*3,tw_max,tw_min),
                                    # n_manning, 
                                    # n_manning_cc, 
                                    normalize(cs[i],cs_max,cs_min),
                                    normalize(s0[i], s0_max, s0_min),
                                    # velp, 
                                    normalize(depthp[i],depthp_max,depthp_min)]]))
    # print(temp_y,temp_y_interp)
    error = abs(temp_y_interp - temp_y[0])
    # print(error, temp_y_interp, temp_y)
    errors = np.append(errors, error)

print(f"Average error is {np.mean(errors)}")
print(f"Max error is {np.max(errors)}")


num_speed_samp = 2700000

dt = 60 # Time step
dx = 1800 # segment length
# bw = np.linspace(0.135, 230.035, array_length, endpoint=True) # Trapezoidal bottom width
bw = np.random.uniform(bw_min, bw_max, num_speed_samp) # Trapezoidal bottom width
tw = np.random.uniform(tw_min, tw_max, num_speed_samp) # Channel top width (at bankfull)
# twcc = np.linspace(0.67, 1150.17, array_length, endpoint=True) # Flood plain width
# twcc = tw*  # Flood plain width tw*3
n_manning = .028   # manning roughness of channel
n_manning_cc = .028 # manning roughness of floodplain
cs = np.random.uniform(cs_min,cs_max, num_speed_samp)# channel trapezoidal sideslope
s0 = np.random.uniform(s0_min, s0_max, num_speed_samp) # Lateral inflow in this time step
qup = np.random.uniform(qup_min, qup_max, num_speed_samp) # Flow from the upstream neighbor in the previous timestep
# quc = np.linsprandom.uniformace(10, 1000, array_length, endpoint=True) # Flow from the upstream neighbor in the current timestep 
quc = np.random.uniform(quc_min, quc_max, num_speed_samp)  # Flow from the upstream neighbor in the current timestep 
# qdp = np.linspace(10, 1000, array_length, endpoint=True) # Flow at the current segment in the previous timestep
qdp = np.random.uniform(qdp_min, qdp_max, num_speed_samp)  # Flow at the current segment in the previous timestep
qlat = np.random.uniform(qlat_min, qlat_max, num_speed_samp) # lat inflow into current segment in the current timestep
velp = .5  # Velocity in the current segment in the previous timestep NOT USED AS AN INPUT!!!
depthp = np.random.uniform(depthp_min ,depthp_max , num_speed_samp) # D

speed_x = []

for i in range(num_samp):
    speed_x.append( [normalize(qup[i],qup_max,qup_min),
    normalize(quc[i],quc_max,quc_min), 
    normalize(qlat[i],qlat_max,qlat_min),
    normalize(qdp[i],qdp_max,qdp_min),
    # dx,  
    normalize(bw[i],bw_max,bw_min),
    normalize(tw[i],tw_max,tw_min),
    # normalize(tw[tw_o]*3,tw_max,tw_min),
    # n_manning, 
    # n_manning_cc, 
    normalize(cs[i],cs_max,cs_min),
    normalize(s0[i], s0_max, s0_min),
    # velp, 
    normalize(depthp[i],depthp_max,depthp_min)])
speed_x = np.array(speed_x)

import time 

starttime=time.time()
regr.predict(speed_x)
print(time.time()-starttime)

# In[ ]:




