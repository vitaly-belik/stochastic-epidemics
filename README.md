# stochastic-epidemics

The code was used to perform simulation for the spread of nosocomial pathogens (such as a resistant MRSA) on referral network of 
hospitals: 

Belik, V., Karch, A., Hövel, P., Mikolajczyk, R. (2017). Leveraging Topological and Temporal Structure of Hospital Referral Networks
for Epidemic Control. In: Masuda, N., Holme, P. (eds) Temporal Network Epidemiology. Theoretical Biology. Springer, Singapore. 
https://doi.org/10.1007/978-981-10-5287-3_9

@incollection{Belik2017,

  author    = {Belik, Vitaly and Karch, Andr{\'e} and H{\"o}vel, Philipp and Mikolajczyk, Rafael},
  
  title     = {Leveraging Topological and Temporal Structure of Hospital Referral Networks for Epidemic Control},
  
  booktitle = {Temporal Network Epidemiology},
  
  editor    = {Masuda, Naoki and Holme, Petter},
  
  series    = {Theoretical Biology},
  
  year      = {2017},
  
  publisher = {Springer},
  
  address   = {Singapore},
  
  doi       = {10.1007/978-981-10-5287-3_9},
  
  url       = {https://doi.org/10.1007/978-981-10-5287-3_9}
  
}

The repository contains a prototyte code for stochastic agent-based simulation of an epidemic spread on a metapopulation network. 
It uses the Gillespie algorithm with a priority queue. Movements of individuals 
between locations are pre-scheduled and given by a temporal network. 
