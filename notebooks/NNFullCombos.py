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
import matplotlib as plt
from pylab import *
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense,Dropout,BatchNormalization
import itertools
import random 
from random import randint
# import time
# start_time = time.time()



# can set AL to a list to run multiple tests on different sized arrays
AL = [4,5]

#output lists if multiple runs were perfermed
mean_errors_list = []
max_errors_list = []
mean_rel_errors_list = []
max_rel_errors_list = []

for size in AL:
# these are the min a max ranges for each input variable. Based on array length specified this will slice these ranges up to be used in our combinations of test data.
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
    array_length = size 

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


    
    #this is a normalization function that is used to scale the input variables between 0-1 to make the network train more efficiently. If not normalized some values may influence the NN too much
    def normalize(val,max,min,target_max=1,target_min=0):
        return (val-min)/(max-min)*(target_max-target_min)+target_min


#MC function that allows us to calculate the expected outputs used as our Y output. Can be used to generate real calculations from the MC code.
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
        
   
#for loops used to find every combination of our variables possible based on our original array size between the mins and maxs of each variable.
    Y = []
    M = []

    # for qup_i in range(0,(array_length)-1):
    #     for quc_j in range(0,(array_length)-1):
    #         for qlat_k in range(0,(array_length)-1):
    #             for qdp_l in range(0,(array_length)-1):
    #                 for bw_n in range(0,(array_length)-1):
    #                     for tw_o in range(0,(array_length)-1):
    #                         for cs_p in range(0,(array_length)-1):
    #                             for s0_q in range(0,(array_length)-1):
                                    # for depthp_s in range(0,(array_length)-1):
    temp = list(itertools.product(qup,quc,qlat,qdp,bw,tw,cs,s0,depthp))
    for i in temp:
        M.append([
            # dt, 
            normalize(i[0],qup_max,qup_min), 
            normalize(i[1],quc_max,quc_min), 
            normalize(i[2],qlat_max,qlat_min),
            normalize(i[3],qdp_max,qdp_min),
            # dx,  
            normalize(i[4],bw_max,bw_min),
            normalize(i[5],tw_max,tw_min),
            # normalize(tw[tw_o]*3,tw_max,tw_min),
            # n_manning, 
            # n_manning_cc, 
            normalize(i[6],cs_max,cs_min),
            normalize(i[7], s0_max, s0_min),
            # velp, 
            normalize(i[8],depthp_max,depthp_min)])
                                        # M.append([dt, qup[qup_i], quc[quc_j], qlat[qlat_k],qdp[qdp_l],dx,  bw[bw_n], tw[tw_o], tw[tw_o]*3,n_manning, n_manning_cc, cs[cs_p], s0[s0_q], velp, depthp[depthp_s]])
                                        # M.append([
                                        #     # dt, 
                                        #     normalize(qup[qup_i],qup_max,qup_min), 
                                        #     normalize(quc[quc_j],quc_max,quc_min), 
                                        #     normalize(qlat[qlat_k],qlat_max,qlat_min),
                                        #     normalize(qdp[qdp_l],qdp_max,qdp_min),
                                        #     # dx,  
                                        #     normalize(bw[bw_n],bw_max,bw_min),
                                        #     normalize(tw[tw_o],tw_max,tw_min),
                                        #     # normalize(tw[tw_o]*3,tw_max,tw_min),
                                        #     # n_manning, 
                                        #     # n_manning_cc, 
                                        #     normalize(cs[cs_p],cs_max,cs_min),
                                        #     normalize(s0[s0_q], s0_max, s0_min),
                                        #     # velp, 
                                        #     normalize(depthp[depthp_s],depthp_max,depthp_min)])

        S = singlesegment(
        dt=dt,
        qup=i[0],
        quc=i[1],
        qlat=i[2],
        qdp=i[3],
        dx=dx ,
        bw=i[4],
        tw=i[5],
        twcc=i[5]*3,
        n_manning=n_manning,
        n_manning_cc=n_manning_cc,
        cs=i[6],
        s0=i[7],
        velp=velp,
        depthp=i[8])
        Y.append((S[0]))
        if len(Y)%100000 == 0:
            print(len(Y))
#this section just adds a test sample with the expected output of .75 to the end of our data in case we would like to compare it
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

#this section randomly generates validation data between the min and the max for us to train against. The model uses its training x and y data to improve itself, but will compare on unseen validation 
#data to make sure it is not overfitting as much 
#takes random number samples between the min/max for each variable
    num_samp = 1000000

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



#creation of the model

    

    # print(M[:100])
    #Define the model
    def baseline_model():
        #creates a sequential model
        model = tf.keras.Sequential()
        #select input node size #, activation function(relu is standard), and input shape being passed from out X/M training array
        model.add(Dense(1028, activation=tf.nn.relu, input_shape=[9],use_bias=False))
    #     model.add(BatchNormalization())
        model.add(Dense(512, activation='relu'))
        # model.add(Dense(128, activation='relu'))
        model.add(Dense(256, activation='relu'))
        model.add(Dense(128, activation='relu'))
    #     model.add(Dense(64, activation='relu'))
    #     model.add(Dense(32, activation='relu'))
    #     model.add(Dense(16, activation='relu'))
        #outputs to a single number with a linear activation function so output can be between -inf,inf
        model.add(Dense(1, activation = 'linear'))
#optimize on mean squared error. can experiment with other loss functions
        model.compile(optimizer = 'adam', loss = 'mse', metrics = ['mse','mean_absolute_error'])
        return model

    regr = baseline_model()
# here is where we input our X/M along with out Y training data, the number of iterations to run through the NN, the number of points to pass the NN at once, and our validation data set
    history = regr.fit(M[:], Y[:], epochs=10, batch_size=100,  validation_data=(VAL_x,VAL_y))

    print(regr.predict(M[-1:]))
    print(Y[-1])



    

#creates random samples to be used in predictions so we can calculate the average error etc.

    Y = []
    M = []
#select the number of points you'd like to sample
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


#you could take real values and insert them here from historical data is you wish (Alex)
    rel_errors = []
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
#calculates errors
        if i%1000 == 0:
            print(i)
        # print(temp_y,temp_y_interp)
        error = abs(temp_y_interp - temp_y[0])
        rel_error = error/temp_y
        rel_errors.append(rel_error)
        # print(error, temp_y_interp, temp_y)
        errors = np.append(errors, error)
    print(f"For a sample of {num_samp}")
    print(f"Average error is {np.mean(errors)}")
    print(f"Max error is {np.max(errors)}")
    print(f"Average relative error is {np.mean(rel_errors)}")
    print(f"Max relative error is {np.max(rel_errors)}")

    mean_errors_list.append(np.mean(errors))
    max_errors_list.append(np.max(errors))
    mean_rel_errors_list.append(np.mean(rel_errors))
    max_rel_errors_list.append(np.max(rel_errors))

print(AL)
print(mean_errors_list)
print(max_errors_list)
print(mean_rel_errors_list)
print(max_rel_errors_list)

    # num_speed_samp = 2700000

    # dt = 60 # Time step
    # dx = 1800 # segment length
    # # bw = np.linspace(0.135, 230.035, array_length, endpoint=True) # Trapezoidal bottom width
    # bw = np.random.uniform(bw_min, bw_max, num_speed_samp) # Trapezoidal bottom width
    # tw = np.random.uniform(tw_min, tw_max, num_speed_samp) # Channel top width (at bankfull)
    # # twcc = np.linspace(0.67, 1150.17, array_length, endpoint=True) # Flood plain width
    # # twcc = tw*  # Flood plain width tw*3
    # n_manning = .028   # manning roughness of channel
    # n_manning_cc = .028 # manning roughness of floodplain
    # cs = np.random.uniform(cs_min,cs_max, num_speed_samp)# channel trapezoidal sideslope
    # s0 = np.random.uniform(s0_min, s0_max, num_speed_samp) # Lateral inflow in this time step
    # qup = np.random.uniform(qup_min, qup_max, num_speed_samp) # Flow from the upstream neighbor in the previous timestep
    # # quc = np.linsprandom.uniformace(10, 1000, array_length, endpoint=True) # Flow from the upstream neighbor in the current timestep 
    # quc = np.random.uniform(quc_min, quc_max, num_speed_samp)  # Flow from the upstream neighbor in the current timestep 
    # # qdp = np.linspace(10, 1000, array_length, endpoint=True) # Flow at the current segment in the previous timestep
    # qdp = np.random.uniform(qdp_min, qdp_max, num_speed_samp)  # Flow at the current segment in the previous timestep
    # qlat = np.random.uniform(qlat_min, qlat_max, num_speed_samp) # lat inflow into current segment in the current timestep
    # velp = .5  # Velocity in the current segment in the previous timestep NOT USED AS AN INPUT!!!
    # depthp = np.random.uniform(depthp_min ,depthp_max , num_speed_samp) # D

    # speed_x = []

    # for i in range(num_samp):
    #     speed_x.append( [normalize(qup[i],qup_max,qup_min),
    #     normalize(quc[i],quc_max,quc_min), 
    #     normalize(qlat[i],qlat_max,qlat_min),
    #     normalize(qdp[i],qdp_max,qdp_min),
    #     # dx,  
    #     normalize(bw[i],bw_max,bw_min),
    #     normalize(tw[i],tw_max,tw_min),
    #     # normalize(tw[tw_o]*3,tw_max,tw_min),
    #     # n_manning, 
    #     # n_manning_cc, 
    #     normalize(cs[i],cs_max,cs_min),
    #     normalize(s0[i], s0_max, s0_min),
    #     # velp, 
    #     normalize(depthp[i],depthp_max,depthp_min)])
    # speed_x = np.array(speed_x)

    # import time 

    # starttime=time.time()
    # regr.predict(speed_x)
    # print(time.time()-starttime)

    # In[ ]:




