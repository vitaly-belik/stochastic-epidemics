# stochastic-epidemics

The repository contains a prototyte code for stochastic agent-based simulation of an epidemic spread on a metapopulation network. 
It uses a Gillespie algorithm with a priority queue (similar to e.g. https://doi.org/10.1021/jp993732q). Movements of individuals 
between locations are given by a temporal network. 

The code was used to perform simulation for the spread of nosocomial pathogens (such as a resistant MRSA) on referral network of 
hospitals: 

Belik, V., Karch, A., Hövel, P., Mikolajczyk, R. (2017). Leveraging Topological and Temporal Structure of Hospital Referral Networks
for Epidemic Control. In: Masuda, N., Holme, P. (eds) Temporal Network Epidemiology. Theoretical Biology. Springer, Singapore. 
https://doi.org/10.1007/978-981-10-5287-3_9
