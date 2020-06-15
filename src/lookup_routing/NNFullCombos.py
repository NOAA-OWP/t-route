
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
import NN_gen_training_data
import NN_normalization 
import nwm_single_segment
import NN_gen_val_data
import NN_gen_predicted_values
# import time
# start_time = time.time()


def main():
    # AL = [5]

    # #output lists if multiple runs were perfermed
    # can set AL to a list to run multiple tests on different sized arrays
    depthp_min = 0.010; depthp_max = .011
    qlat_min = 35; qlat_max = 45
    qdp_min = .01; qdp_max = 1
    quc_min = .01; quc_max = 1
    qup_min = .01; qup_max = 1
    s0_min = .00001; s0_max = .002; 
    cs_min = .085; cs_max = 2.254;
    tw_min = 150; tw_max = 500;
    bw_min = 112; bw_max = 150;
    
    (M,
    Y,
    mean_errors_list,
    max_errors_list,
    mean_rel_errors_list,
    max_rel_errors_list,
    AL)= NN_gen_training_data.main(
    depthp_min,
    depthp_max,
    qlat_min,
    qlat_max,
    qdp_min,
    qdp_max,
    quc_min,
    quc_max,
    qup_min,
    qup_max,
    s0_min,
    s0_max,
    cs_min,
    cs_max,
    tw_min,
    tw_max,
    bw_min,
    bw_max
    )
    
    (VAL_x,
    VAL_y) = NN_gen_val_data.main(
    depthp_min,
    depthp_max,
    qlat_min,
    qlat_max,
    qdp_min,
    qdp_max,
    quc_min,
    quc_max,
    qup_min,
    qup_max,
    s0_min,
    s0_max,
    cs_min,
    cs_max,
    tw_min,
    tw_max,
    bw_min,
    bw_max
    )
    #this section randomly generates validation data between the min and the max for us to train against. The model uses its training x and y data to improve itself, but will compare on unseen validation 
    #data to make sure it is not overfitting as much 
    #takes random number samples between the min/max for each variable
    

    


    


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
    history = regr.fit(M[:], Y[:], epochs=300, batch_size=1000,  validation_data=(VAL_x,VAL_y))

    print(regr.predict(M[-1:]))
    print(Y[-1])

    NN_gen_predicted_values.main(
    depthp_min,
    depthp_max,
    qlat_min,
    qlat_max,
    qdp_min,
    qdp_max,
    quc_min,
    quc_max,
    qup_min,
    qup_max,
    s0_min,
    s0_max,
    cs_min,
    cs_max,
    tw_min,
    tw_max,
    bw_min,
    bw_max,
    regr,
    AL)

main()

