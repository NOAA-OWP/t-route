import sys
import numpy as np

debuglevel = 0
COMPILE = True
if COMPILE:
    try:
        #Assume fortran module is compiled in same folder as calling python script
        #If not, a sys.path.append() call is needed to add the module path
        import subprocess
        fortran_compile_call = []
        fortran_compile_call.append(r'f2py3')
        fortran_compile_call.append(r'-c')
        fortran_compile_call.append(r'varSingleChStime_f2py.f90')
        fortran_compile_call.append(r'MCsingleChStime_f2py_clean.f90')
        fortran_compile_call.append(r'-m')
        fortran_compile_call.append(r'mc_sc_stime')
        #subprocess.run(fortran_compile_call, cwd=fortran_source_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if debuglevel <= -2: 
            subprocess.run(fortran_compile_call)
        else:
            subprocess.run(fortran_compile_call, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        import mc_sc_stime as mc
    except Exception as e:
        print (e)

#TODO: generalize with a direction flag
def compute_mc_reach_up2down(
        head_segment = None
        , reach = None
        , reach_connections = None
        , reach_flowdepthvel = None
        , upstream_inflow = 0
        , supernetwork_data = None
        , ts = 0
        , verbose = False
        , debuglevel = 0
        ):

    if verbose: print(f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments'])})")
    
    # filename = f'./logs/{head_segment}_{ts}.log'
    # file = open(filename, 'w+') 
    
    ntim=2;       #the number of time steps necessary for variables passed to mc module to compute successfully
    nlinks=2;     #the number of links needed to define varialbe qd. ** nlinks is not used in fortran source code.

    mc.var.uslinkid=1
    mc.var.linkid=2
    ncomp0=1; mc.var.ncomp0=ncomp0  #the number of segments of a reach upstream of the current reach
    ncomp = len(reach['segments']) ;  mc.var.ncomp= ncomp  #the number of segments of the current reach 
    #mxseg=max(ncomp0,ncomp)    
    #MC model outputs
    mc.var.qd=np.zeros((ntim,ncomp,nlinks))  #will store MC output qdc (flow downstream current timestep) 
    mc.var.vela=np.zeros((ntim,ncomp)) 
    mc.var.deptha=np.zeros((ntim,ncomp))
    #lateral flow
    mc.var.qlat=np.zeros((ncomp))
    mc.var.dx=np.zeros((ncomp))
    mc.var.bw=np.zeros((ncomp))
    mc.var.tw=np.zeros((ncomp))
    mc.var.twcc=np.zeros((ncomp))
    mc.var.ncc=np.zeros((ncomp))
    mc.var.cs=np.zeros((ncomp))
    mc.var.so=np.zeros((ncomp))
    mc.var.n =np.zeros((ncomp))
    
    
    mc.var.qd[0,0,0]= upstream_inflow
            
    current_segment = reach['reach_head']
    
    # writeString = f'timestep: {ts} cur : {current_segment}  {upstream_inflow}'
    # writetoFile(file, writeString)
    # writeString = f'variables'
    # writetoFile(file, writeString)

    next_segment = reach_connections[current_segment]['downstream'] 
    i = 0
    #input flow to upstream reach of current reach     
    while True:
        # Thanks to SO post for a reminder of this "Loop-and-a-half" construct
        # https://stackoverflow.com/questions/1662161/is-there-a-do-until-in-python

        # for now treating as constant per reach 
        dt=300.0 ;      mc.var.dt= dt  #60.0;

        #TODO: James Ask Dong Ha/Juzer: Does i == 0 have any meaning for these fortran arrays?
        mc.var.dx[i]=reach_connections[current_segment]['data'][supernetwork_data['length_col']] 
        mc.var.bw[i]=reach_connections[current_segment]['data'][supernetwork_data['bottomwidth_col']] # ;        mc.var.bw= bw #50
        mc.var.tw[i]=0.01*mc.var.bw[i] # ;  mc.var.tw= tw
        mc.var.twcc[i]=mc.var.tw[i] # ;      mc.var.twcc=twcc
        mc.var.n[i]=reach_connections[current_segment]['data'][supernetwork_data['manningn_col']]  # ;      mc.var.n=n #0.03
        mc.var.ncc[i]=mc.var.n[i] # ;      mc.var.ncc=ncc
        mc.var.cs[i]=reach_connections[current_segment]['data'][supernetwork_data['ChSlp_col']]  # ;    mc.var.cs=cs #1.0e6
        mc.var.so[i]=reach_connections[current_segment]['data'][supernetwork_data['slope_col']] # ;    mc.var.so=so #0.002
        #so = max(0.002, so)
        #ck= current_segment['data'][supernetwork['MusK_col']];   mc.var.ck = ck 
        #cx= current_segment['data'][supernetwork['MusX_col']];   mc.var.cx = cx
        #print (f'{current_segment}')
   
        #tmp allocation
        dx = mc.var.dx[i]
        bw = mc.var.bw[i]
        tw = mc.var.tw[i]
        n = mc.var.n[i]
        cs = mc.var.cs[i]
        so = mc.var.so[i]
        
        #print (f'counter = {i}')
        #if current_segment == 5559368 or i == 100:

        mc.var.qlat[i]= reach_flowdepthvel[current_segment]['qlat']['curr']  # temporary assigned qlat 
        mc.var.qd[0,i,1]= reach_flowdepthvel[current_segment]['flow']['prev']  # temporary assigned qd
        mc.var.vela[0,i] = reach_flowdepthvel[current_segment]['vel']['prev']
        mc.var.deptha[0,i] = reach_flowdepthvel[current_segment]['depth']['prev']
        
        # writeString = f'timestep: {ts} cur : {current_segment}  {dx} {bw} {tw} {n} {cs} {so} {dt}'
        # writetoFile(file, writeString)
        # writeString =  f"Previous: {reach_flowdepthvel[current_segment]['flow']['prev']} "
        # writeString =  writeString + f"{reach_flowdepthvel[current_segment]['depth']['prev']} "
        # writeString =  writeString + f"{reach_flowdepthvel[current_segment]['vel']['prev']} "
        # writeString =  writeString + f"{reach_flowdepthvel[current_segment]['qlat']['curr']}"
        # writetoFile(file, writeString)

        i += 1
        
        if current_segment == reach['reach_tail']:
            if verbose: print(f'{current_segment} (tail)')
            break
        if verbose: print(f'{current_segment} --> {next_segment}\n')
        current_segment = next_segment
        next_segment = reach_connections[current_segment]['downstream'] 
    #end loop initialized the MC vars 
    
    # writeString = f'data'
    # writetoFile(file, writeString)
  
    mc.mc.main()

    #print (f'{ts} end mc')
    
    current_segment = reach['reach_head']
    next_segment = reach_connections[current_segment]['downstream'] 
    i = 0
    while True:
        reach_flowdepthvel[current_segment]['flow']['curr'] = mc.var.qd[1,i,1] 
        reach_flowdepthvel[current_segment]['depth']['curr'] = mc.var.deptha[1,i]
        reach_flowdepthvel[current_segment]['vel']['curr'] = mc.var.vela[1,i]
        d = reach_flowdepthvel[current_segment]['depth']['curr'] 
        q = reach_flowdepthvel[current_segment]['flow']['curr']
        v = reach_flowdepthvel[current_segment]['vel']['curr']
        ql = reach_flowdepthvel[current_segment]['qlat']['curr']
        # writeString = f'timestep: {ts} cur : {current_segment}  {q} {d} {v} {ql} '
        # writetoFile(file, writeString)
        i += 1
        
        if current_segment == reach['reach_tail']:
            if verbose: print(f'{current_segment} (tail)')
            break
        if verbose: print(f'{current_segment} --> {next_segment}\n')
        current_segment = next_segment
        next_segment = reach_connections[current_segment]['downstream'] 
    #end loop collect MC output 
    return {head_segment:reach_flowdepthvel}

def main():
    pass
#TODO: PUT 4 link test Implemenation Here
'''
Call Musk Cunge with test parameters. 
When importing, call with parameters from NHD network
'''

if __name__ == '__main__':
    main()
