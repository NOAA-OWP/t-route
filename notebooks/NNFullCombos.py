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
# singlesegment():
array_length = 11    

dt = 60 # Time step
dx = 1800 # segment length
# bw = np.linspace(0.135, 230.035, array_length, endpoint=True) # Trapezoidal bottom width
bw = np.linspace(112.000000000000000000, 150.000000000000000000, array_length, endpoint=True) # Trapezoidal bottom width
tw = np.linspace(150.000000000000000, 500.000000000000000, array_length, endpoint=True) # Channel top width (at bankfull)
# twcc = np.linspace(0.67, 1150.17, array_length, endpoint=True) # Flood plain width
# twcc = tw*  # Flood plain width tw*3
n_manning = .028   # manning roughness of channel
n_manning_cc = .028 # manning roughness of floodplain
cs = np.linspace(0.085000000000000000, 2.254000000000000000, array_length, endpoint=True)# channel trapezoidal sideslope
s0 = np.linspace(0.000010000000000, .002000000000000000, array_length, endpoint=True) # Lateral inflow in this time step
qup = np.linspace(.010000000000000000, 1.0000000000000000, array_length, endpoint=True) # Flow from the upstream neighbor in the previous timestep
# quc = np.linspace(10, 1000, array_length, endpoint=True) # Flow from the upstream neighbor in the current timestep 
quc = np.linspace(.010000000000000000, 1.0000000000000000, array_length, endpoint=True)  # Flow from the upstream neighbor in the current timestep 
# qdp = np.linspace(10, 1000, array_length, endpoint=True) # Flow at the current segment in the previous timestep
qdp = np.linspace(.010000000000000000, 1.000000000000000, array_length, endpoint=True)  # Flow at the current segment in the previous timestep
qlat = np.linspace(35.0000000000000000, 45.00000000000001, array_length, endpoint=True) # lat inflow into current segment in the current timestep
velp = np.linspace(0.050000000000000000, .10000000000000000, array_length, endpoint=True) # Velocity in the current segment in the previous timestep NOT USED AS AN INPUT!!!
depthp = np.linspace(0.01000000000000000 ,.0110000000000000000 , array_length, endpoint=True) # D



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

for i in range(0,(array_length)-1):
    for j in range(0,(array_length)-1):
        for k in range(0,(array_length)-1):
            for l in range(0,(array_length)-1):
                for n in range(0,(array_length)-1):
                    for o in range(0,(array_length)-1):
                        for p in range(0,(array_length)-1):
                            for q in range(0,(array_length)-1):
                                for r in range(0,(array_length)-1):
                                    for s in range(0,(array_length)-1):
                                        M.append([dt, qup[i], quc[j], qlat[k],qdp[l],dx,  bw[n], tw[o], tw[o]*3,n_manning, n_manning_cc, cs[p], s0[q], velp[r], depthp[s]])
                S = singlesegment(
                dt=dt,
                qup=qup[i],
                quc=quc[j],
                qlat=qlat[k],
                qdp=qdp[l],
                
                dx=dx ,
                bw=bw[n],
                tw=tw[o],
                twcc=tw[o]*3,
                 n_manning=n_manning,
                n_manning_cc=n_manning_cc,
                cs=cs[p],
                s0=s0[q],
                velp=velp[r],
                depthp=depthp[s])
                Y.append(round(S[0],4))
                

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

M.append([dt, qup, quc, qlat,qdp,dx,  bw, tw, twcc,n_manning, n_manning_cc, cs, s0, velp, depthp])


M = np.array(M)
for i in range(0,len(M),1):
    S = singlesegment(*M[i])
    Y.append(S[0])
Y = np.array(Y)    


print(Y[-1])
print(M[-1])

from pylab import *
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense,Dropout,BatchNormalization
# import time
# start_time = time.time()



#Define the model
def baseline_model():
    model = tf.keras.Sequential()
    model.add(Dense(512, activation=tf.nn.relu, input_shape=[15],use_bias=False))
#     model.add(BatchNormalization())
    model.add(Dense(256, activation='relu'))
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

history = regr.fit(M[:TOP], Y[:TOP], epochs=10, batch_size=10,  validation_split=0.1)
plt.plot(history.history['mse'])
plt.plot(history.history['val_mse'])

plt.title('Model accuracy')
plt.ylabel('MSE')
plt.xlabel('Epoch')
plt.legend(['MSE', 'val_mse'], loc='upper left')
plt.show()
# print("--- %s seconds ---" % (time.time() - start_time))
regr.predict(M[-1:])
Y[-1]


# In[ ]:




