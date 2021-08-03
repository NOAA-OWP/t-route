
_Question:_
How do you perform the domain partitioning results for parallelization?

_Answer:_
To identify the independent networks, it is a simple DFS started at all points in the network that have a "to" label not contained within the network. I.e., we start at the outlet segments of the network and search upstream until we find all members of the collection and tie them to their corresponding outlet.
That establishes the first level by-network parallelization units.
To break things down further, we use a threshold limited BFS search on the largest networks. The BFS is started on the outlet segments (like the first-level DFS) and the total number of segments encountered inside the search front is tracked until it exceeds a set threshold. Once the threshold is reached, the search is terminated at the next junction along each search path and then the search is recursively restarted at that search fringe. The closest thing to a name that I came up with for that is based on the way we sequence the calculation of those subnetworks -- JIT -- just in time, vs. Opportunistic parallelization. Opportunistic parallelization would be to do any of the headwaters first, then to do any networks for which all precedents have been calculated. The JIT means that we look forward one level and anticipate what will be needed next along the longest dependency tree and compute now what will be needed as a dependency next -- hence just-in-time. This saves some of the headwaters that might connect into a larger mainstem further down the network for a little later, allowing some load balancing.

_Question:_
Fred had some misunderstanding on my AGU abstract, so we were working to clarify that we want to apply previous decomposition methods in the routing domain to the catchment domain, so I'm just looking for something to fill in that (application of what method) here.  So in synopsis, (application of what method) == just in time decomposition, a combination of depth first search and a threshold limited breadth first search
So the sentence would read:Prior studies of hydrologic (e.g. Muskingum) stream flow routing models show that an optimal domain partitioning results from just in time decomposition, a combination of depth first search and a threshold limited breadth first searchDoes that hold up?
Or to view it in full context: https://docs.google.com/document/d/1TwsHQK2m4jXhiBt_3kzS7ofqHMbfX7jMUHd535c7yP4/edit

_Question:_
How does the subdividing the largest elements via a threshold-limited breadth first search relate to the just-in-time parallel approach?

_Answer:_
The parallelization via just-in-time traversal of the ordered subnetworks is what makes the decomposition optimal.
If you parallelized another way, it becomes less optimal -- maybe even worse than serial.

_Question:_
So the subnetworks identified by the BFS are not executed by the same resource?

_Answer:_
No, they are not. The BFS is a single core process (though the potential exists to parallelize it) which performs the domain decomposition in a single pass and stores the information. The output of the BFS gives an ordering (related to the subnetwork) in which to execute
and things in the entire network of the same order can be run in parallel.
Touching each segment to provide the decomposition is super low cost, so it doesn't make sense to parallelize it.

_Question:_
Then what happens in execution?

_Answer:_
We run through the collection of subnetworks from highest to lowest order (as Identified by the recursions of the BFS) with the guarantee that 1) if we do it in that order, the precedents of every calculation will be satisfied before the calculation is performed; and 2) that the precedents will be satisfied just before they are needed, i.e., in the prior order.
As far as the "optimality" of the solution, there is definitely a threshold size that optimizes the tradeoff between computation burden per subnetwork + parallel overhead AND the opportunity for parallel execution with increased number of parts.

_Question:_
How did you prove that?

_Answer:_
Well, this is applied math, somewhere in here, but there's not a math proof. 

But one more thing (Keep your pencil out) as far as optimality, there is a longest path length that has to be traversed which represents an absolute minimum number of computation cycles -- i.e., if we could assign a parallel node to every calculation, there would still have to 2219 sequential segments to get through the longest path in the CONUS network. The JIT solution with the decomposition as described tends to keep the overall number of required calculation cycles close to this theoretical minimum without requiring any complex message passing or other external management.
Or rather, it tends to use each cycle very efficiently, so that a fewer parallel nodes can perform the computation as easily as a greater number with less efficient utilization.

_Question:_
What was your process for figuring this all out?

_Answer:_
The main challenge of our project was to re-introduce the topological dependence without breaking the computational bank. 

To think about how to do that, we looked at a couple things we know about the problem. We know there is a longest path that has to be run in sequence -- that's one hard limit. And we know there is a set number of segments -- getting through all of those is another hard limit. With those, we can set some absolute limits on wall time.

For a given calculation speed, we can establish a curve that represents the minimum theoretical calculation wall time for different numbers of parallel cores -- the longest path sets the absolute floor. 

Then I did two spreadsheet experiments. 
1) Assuming no parallel overhead and assuming we figure out how to get each parallel unit to be the same computational size, what is the length of computation for different number of computation cores based on different ways of farming the parallel units to the cores? This was the comparison of Just-In-TIme, Opportunistic, and partial opportunistic approaches. This showed me that how we choose to order the calculation can either optimize or not the use of the parallel compute resources. For a given topological decomposition, the calculation ordering that tended to ensure consistent maximum core usage was the JIT. Other methods had a tendency to saturate the cores early but to leave many cores idle at the end of the computation sequence. The JIT isn't 100% optimal, but I was surprised at how well it did without requiring any complex forward evaluation or message passing. 
2) The absolute minimum assumes no parallel overhead. The second spreadsheet exercise was to assume a certain amount of overhead per cycle and an overhead for each new parallel call. That is what gave us the guaranteed optimum somewhere between doing everything all at once and splitting things up to the tiniest possible decomposition. On one side, the lack of parallelism dominates. On the other side, the wall time pushes back up because of compounding inefficiencies in parallelization. The true optimum comes from interaction of the potential parallel scaling (see experiment 1 -- if you have a very narrow, long network, you won't get much of a speedup no matter how many cores you add) and the ratio of the overhead amounts to the calculation time (e.g., if you have really big parallel costs, the solution pushes closer to the serial); and to one another (e.g. a big per-cycle overhead doesn't impact much if the network depth is small.) 
