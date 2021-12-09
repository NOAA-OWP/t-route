import pandas as pd
import numpy as np
import hydroeval as he
import math 
# import matplotlib.pyplot as plt

s1 = pd.read_pickle("/glade/u/home/jhreha/t-route/test/jobs/parity_check/chanobs_output/example1")
s2 = pd.read_pickle("/glade/u/home/jhreha/t-route/test/jobs/parity_check/chanobs_output/example2")

s1_q = s1.iloc[:,::3].sort_index()
s2_q = s2.iloc[:,::3].sort_index()

print("\n")
print("****************************************")
print("****** Are simulated flows equal? ******")
print("****************************************")
print("\n")
test = np.array_equal(s1_q.to_numpy().astype(float), s2_q.to_numpy().astype(float), equal_nan = True)


if test:
    print("YES! - Simulations produce the same results")
    
else:
    print("NO. - Simulations produce different results")
    diff = abs(s1_q - s2_q)
    col_idx = diff.max().argmax()
    row_idx = diff.iloc[:,col_idx].argmax()
    
    print("The largest flow difference in flow is: ", diff.max().max(), "cms")
    print("at link ID", s1_q.index[col_idx], "and timestep", row_idx)
    print("simulation 1 reports a flow value of", s1_q.iloc[row_idx, col_idx])
    print("simulation 2 reports a flow value of", s2_q.iloc[row_idx, col_idx])
    print("\n")
       
    print("****************************************")
    print("**** Are simulated flows similar? ******")
    print(" ** (... within 0.1% of one another) ***")
    print("****************************************")
    print("\n")
    diff_frac = abs(s1_q - s2_q)/s1_q
    col_idx = diff_frac.max().argmax()
    row_idx = diff_frac.iloc[:,col_idx].argmax()

    if diff_frac.max().max() < 0.001:
        print("YES! - Simulation results are very close")
    else:
        print("NO. Simulation results are unacceptably different.")
        print("The largest proportional flow difference is: ", diff_frac.max().max() * 100, "%") 
        print("at link ID", s1_q.index[col_idx], "and timestep", row_idx)
        print("simulation 1 reports a flow value of", s1_q.iloc[row_idx, col_idx])
        print("simulation 2 reports a flow value of", s2_q.iloc[row_idx, col_idx])

    overall_nse = he.evaluator(he.nse,s2_q.to_numpy().flatten(), s1_q.to_numpy().flatten())
    nse_list = []
    for id in s2_q.index: nse_list.append(he.evaluator(he.nse,s2_q.loc[id].dropna().to_numpy().flatten(), s1_q.loc[id].dropna().to_numpy().flatten())[0])
    nse_list_filtered = [x for x in nse_list if math.isnan(x) == False]
    print("Overall nse:", overall_nse)
    print("Minimum nse:", min(nse_list_filtered))
    # plt.hist(nse_list_filtered, density=True, bins=30)
    # plt.ylabel('NSE')
    # plt.xlabel('Data')
    # plt.suptitle("Florence Serial vs Parallel NSE by ID")
    # plt.show()