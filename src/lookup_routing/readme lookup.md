# Neural Network Substitution-based routing for the National Water Model
version 3.0

# Description
The National Water Model (NWM) is a hydrologic model that simulates observed and forecast
streamflow of approximately 2.7 million stream segments in the continental US (CONUS). The
channel routing computation of the current operational NWM is based on the Muskingum-Cunge
method (MC), which solves its two parameters using a numerical trial-and-error method called
the secant method. These iterations can be computationally expensive, requiring approximately one
core hour to complete five minutes of simulation-time the channel routing calculations in the model..
We propose the use of a machine learning method to solve nonlinear partial differential routing
equations. By training a neural network to reproduce the hydraulic routing computation, we
expect to generate ensembles of predictions that are interchangeable with the output of the
traditional method in a fraction of the time required. For example, elements of the network that
are currently in a receding condition are expected to quickly reach a solution due predictable
patterns of relaxation. By accelerating the rate in which the NWM’s hydraulic routing
computation is performed, we explore the impact to the overall computational burden of the
model.
We have created a recurrent neural network model that imitates the output of a formulated MC
single segment and Manning’s equation. Additionally, we have confirmed parity of the single
segment MC calculation with the native WRF-Hydro version. Evaluations of the impacts of each
input parameter to the MC single segment calculation were also assessed during this process.
We developed modifications to the model to make the prototyping of parameters more efficient.
Finally, we explored if the application of neural networks could be applied to other
computationally heavy processes in the NWM.

# Installation
Install Python3
Install Tensorflow
Install Keras
Run the file located at  https://github.com/NOAA-OWP/t-route/blob/master/src/fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS/mc_sseg_stime_NOLOOP_demo.py to generate a required file.

# Usage
Substitute for MC single segment

# Components
Main - NNFullCombos.py :
Run this file to generate the model and predictions.

**NN_argpase.py** :
List of default prototyping values for command line input. Ex.1) `python3 NNFullCombos.py --epochs 5` will change the default number of epochs from 3 to 5.
Ex.2) `python3 NNFullCombos.py --batch_size 5000` will change the batch size from a default of 1000 to 5000 for this run. 

**NN_gen_training_data.py** : 
Generates training data based on the min to max of each parameter and a specified number of slices between the two. For example, if the min=1 and max=3 for
a specific input paramater A and the array size is 3 then it will create a list [1,2,3] for input A. If input B has a min=4 and max=6 it will create a list of [4,5,6]. Once all arrays 
are created for each input every possible combination of inputs is created from these lists. [1,4] [1,5], [1,6], [2,4] ... and so on. Once all the possible combinations of inputs that the model may see in the 
real quadratic formula the model is mimicking all of these inputs are run through the SingleSeg.py function. The goal is to replace the SingleSeg.py function by obtaining each of the related inputs and outputs to train on. This creates a large
search space of inputs and their outputs for the model to train upon.

**NN_gen_val_data.py** :
Follows the same methodology as NN_gen_training_data.py, but instead creates a random selection of inputs from within the search min/max rather than equal search combinations based on the array size selected.

**NN_normalization.py** :
This file normalizes the data used to train the model.

def normalize(val, max, min, target_max=1, target_min=0):
    return (val - min) / (max - min) * (target_max - target_min) + target_min

**nwm_single_segment.py** :
Contains the SingleSegment function used to produce the outputs the model is trying to mimick and train against.

**NN_gen_predicted_values.py** :
Creates a new set of random inputs and tests them against real outputs run through the SingleSeg.py function. Produces error rates between predicted 
and expect once complete. 

Prototypable inputs :
parser.add_argument(
        "-e",
        "--epochs",
        help="Set the number of epochs (e > 0)",
        dest="epochs",
        default=3,
    )
    parser.add_argument(
        "-b",
        "--batch_size",
        help="Set the batch size (b > 0)",
        dest="batch_size",
        default=1000,
    )
    parser.add_argument(
        "-v",
        "--validation_samples",
        help="Set the validation sample size",
        dest="num_samp_val",
        default=100000,
    )
    parser.add_argument(
        "-p",
        "--prediction_samples",
        help="Set the prediction sample size",
        dest="num_samp_pred",
        default=1000,
    )
    parser.add_argument(
        "-a",
        "--array_length",
        help="Set the array length or number of slices between the min and max of each input parameter",
        dest="AL",
        default=4,
    )

# Support
https://github.com/NOAA-OWP/t-route/

# Contributing
https://github.com/NOAA-OWP/t-route/

Open to the public 

# Authors and acknowledgment
Jacob E. Hreha, James S. Halgren,  Adam N. Wlostowski,  Alexander A. Maestre,  kumdonoaa

# Project status
Looking to improve the model accuracy beyond its current threshold. Model is reaching a limit of about +-.01 and we would like to see that improved by at least one order of magnitude. 
Not sure if the model is rounding heavily or incapable of obtaining that level of accuracy because of a backend issue. 

Author:      jhrehanoaa, Lynker Technologies/NOAA
Created:     7/28/2020

National Water Center 
Office of Water Prediction
205 Hackberry Ln, Tuscaloosa, AL 35401

