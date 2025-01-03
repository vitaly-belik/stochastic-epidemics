#!/usr/bin/python

"""
Agent-based simulation for epidemics on a set of locations(hospitals). Travel
events are read from an input file.
(Vitaly Belik)

"""

from utils import *
import sys
import random
from random import choice
import math
import numpy as np
import Queue
import re
import datetime

"""
    21.10.2016 screening fraction now in the pars file in util.py.

    22.09.2016 this version is for screening every time upon admission
    i.e. the status of a screened ("-1") patient changes from "-1" (if "-1" then
    infected) to "1" upon admission.
    Also note there is no difference between where from a patient is coming, from
    home or another hospital for screening and "forgetting" about screening.


    25.07.2014
    - version, where the health state of incoming patients is checked and their
    prevalence is calculated just as fraction of infected guys.


    22.07.2014
    - version with newly arriving into the system patients infected with some
    probability (0.2%).


    21.07.2014

    - patients have a small prevalence in the community. Every time a patient
    comes to the hospitals with the prob. of apprx. 2% of being infected

    - now the appearence and disappearenc of patients in the system is implemented
    to lower the memory requirements (could this reduce the run time?)



    17.07.2014

    - now the increase in clearance rate is included. When the recovery time is
    assigned with probability of xi it could multiplied by either 1 or 0.5 etc.

    16.07.2014

    - now screening is implemented using the "recovery" state. Guys who are screened
    upon admission as carrier of the disease, become "recovered" and thus do not
    contribute to the dynamics anymore. The only issue could be their readission
    and thus a renewed carriage of the disease. One can monitor this.


    15.07.2014

    - screening (100% effective) at particular location is introduced.
    - starting with infecteds everywhere to achive the endemic level sooner.


    26.02.2014
    now we have an sis model and want to obtain the steady state prevalence

"""


chunck_size = 1000



g = {}

# intervention_time = 10000

inf_prob = 0.0 # prob. to infect an incoming patient in the community


# -------------------------------------------------------------------------
if __name__ == '__main__':

    import sys
    import string
    import re

    # reading the parameters or arguments

    av = args(sys.argv)

    random.seed(av.seed)

    # reading (initial) travel events into the queue
    Q = Queue.PriorityQueue()

    # ---------------------------------------------------------------------

    hosps = {}
    # real guys starting from the community
    # are all of them actually start from the community?

    a = np.loadtxt( av.hids_file ) # in hids_file there is just a list with available
    # HIDs
    for el in a:
        HID = int(el)
        hosps[HID] = hospital(0,set(),HID,g)

    hosps[0] = hospital(0,set(),0,g)


    for el in hosps:
        hosps[el].sus = set(hosps[el].patients) # otherwise would copy by reference
        hosps[el].inf = set()
        hosps[el].rec = set()



    now = datetime.datetime.now()
    f = open(av.file, "r")
    #outf = open(av.out_file + re.sub(" ","_",str(now)) + '.txt',"w")
    outf = open(av.out_file + '.txt', "a") # append modus, because we already wrote the
    # parameters into the out-file


    top_nodes = set()
    a = np.loadtxt( av.top_nodes_file ) # in hids_file there is just a list with available
    for el in a:
        top_nodes.add(int(el))


    # reading and sheduling travel events
    # Initializing the queue
    Q = InitQueue(Q, f, chunck_size, inf_prob)


    chunck_counter = 1

    # initialization

    # here to add initial infections in every hospital
    # to write

    #e = event(start, -1, init_hosp, init_hosp, 'I')

    # initial infection at time "start"!

    for h in hosps:
        e = event(av.start, -1, h, h, 'I') # <-- At time moment "start" there will
        # be an infection event in the hospital "h"
        Q.put( (e.t, e) )


    print "start", av.start


    # for h in sorted(hosps.keys()):
    #     print h


    # print sys.exit()

    # actual simulation!!!
    myeof = 0
    t = 0
    while t < av.max_t:


        if Q.empty() == False:
            e = Q.get()[1]
            Q.task_done()

        else:
            break

        # print "running", e.t, e.type

        fr = e.FROM
        to = e.TO

        if e.type != 'I' and e.type != 'NEXT_CHUNCK' and e.type != 'IN' and e.type != 'OUT': #i.e. type = "R" ot "T"
            try:
                l = g[e.PID].loc
            except KeyError:
                pass
                #print "KeyError", e.type, e.PID
                #sys.exit()


        f1 = 0  # indicates
        f2 = 0
        f3 = 0
        #print "here"

        if e.type  == 'NEXT_CHUNCK':
            chunck_counter += 1

            # print "NextChunck", chunck_counter

            # we combine sampling with reading the next chunck
            ProduceOutput(outf, hosps, t)

            # read the next chunck of travel events into the queue
            Q = ReadNextChunck(Q, f, chunck_size, inf_prob, myeof, av.NOI)


        else:

            if e.type == 'IN':

                if e.FROM == -3: # it was thought that one can infect only some
                                 # given hospitals
                    g[e.PID] = guys(1,0) # creating a new infected guy
                else:
                    g[e.PID] = guys(0,0) # creating a new susceptible guy

                try:
                    hosps[to].into(g, e.PID)
                except KeyError:
                    print to, "pid ", e.PID, (to in hosps)
                    sys.exit()


                # hosps[to].into(g, e.PID)
                l = g[e.PID].loc

            if e.type == 'OUT': # if we destroy

#                try:
                hosps[fr].out(g, e.PID)
#                except KeyError:
#                    print t, e.PID
#                if e.PID == 213024:
#                    print "213024 already del"
                del g[e.PID]


            if e.type == 'I':   # for initial infection, otherwise there are no such
                                # events anymore
                InitialInfection(Q, e, av.NOI, av.beta, av.start, hosps, g)

            else:
                if e.type == 'T':

                    # update means updating everything before e.t

                    # manipulate source                                 <=======
                    if fr != 0: # no dynamics in the community
                        f1 = hosps[fr].update(g, e, Q, av) # normally returns 1
                    else:
                        f1 = 1


                    # manipulate target                                 <=======
                    if to != 0: # no dynamics in the community
                        f2 = hosps[to].update(g, e, Q, av) # normally returns 1
                    else:
                        f2 = 1
                    if f1 != 1 or f2 != 1:  # if rt < e.t
                        Q.put( (e.t,e) ) # put it back?

                else:
                    # e.type == R; infections before recovery
                    if e.type != 'IN' and e.type != 'OUT' and e.type != 'NEXT_CHUNCK':

                        if l != 0:                      # <=======
                            f3 = hosps[l].update(g, e, Q, av) #?????
                        else:
                            f3 = 1
                        if f3 != 1:  # if rt < e.t
                            Q.put( (e.t,e) )
                        

                #  -----------    R E C O V E R Y   ----------------------

                # actual recovery (becoming healthy)
                # we need to take care of patients scheduled for recovery, but
                # who are now "screened", we cannot make them "susceptible" again
                
                if e.type == 'R' and f3 == 1: # rt > e.t
                
                    t = e.t
                    # checking of the patient is really still infected
                    # could be that the patient was already cured or isolated
                    # if 1:
                    # what
                    
                
                        
                    try:
                        zustand = g[e.PID].state

                        # HERE IS THE PART ADDED FOR SCREENING WITH FORGETTING

                        if zustand > 0: # <-----
                            # i.e. really   INFECTED
                            
#                            if e.PID == 213024:
#                                print "recovery 213024 happens"
#                            

                            g[e.PID].state = 0      # becoming healthy again!

                            try:
                                hosps[l].inf.remove(e.PID)
                            except ValueError:
                                sys.exit()
                            if e.PID == 1457:
                                print "here update sus", l
                            hosps[l].sus.add(e.PID)
                    
                        # Here the scenario A, always

                        # The lines below make sense only if we consider
                        # decolonization, otherwise the screened patients
                        # always remain isolated and thus does not play any
                        # role in the infectious dynamics anymore
                        
                        if zustand == -1: # <-----
                            if av.flag_scr_immer == 0:
                            
                                # i.e. SCREENED

                                # here we could take care of patients always remaining
                                # isolated (we do not allow them to become healthy again)
                                # or we decolonize them and do not consider them screened
                                # any more

                                g[e.PID].state = 0      # becoming healthy again!

                                try:
                                    hosps[l].rec.remove(e.PID)
                                except ValueError:
                                    print "Value Error", e.PID
                                    sys.exit()
                            
                                hosps[l].sus.add(e.PID)
                        
                            else:
                                pass
                        
                    except:
#                        print "exception"
                        pass # if we already took the patient out
                        # e.g. recovery was sheduled for the time after
                        # disappearance from the system
                        # do not forget output in the "infection" routine!

                # actual travel
                if e.type == 'T' and f1 == 1 and f2 == 1:  # travel
                
                    t = e.t # time of travel
                    
                    # 26.09.2016 the problem is in sus etc. counting
                    
#                    if e.PID == 213024:
#                        print "t, in, state, h, is rec, is inf, -----", t, (e.PID in hosps[fr].patients),\
#                        g[e.PID].state, fr, (e.PID in hosps[fr].rec), (e.PID in hosps[fr].inf)

#                    try:
                    hosps[fr].out(g, e.PID)
#                    except KeyError:
#                        print "time", t

                    
                    # try:
                    #     hosps[fr].out(g, e.PID)
                    # except KeyError:
                    #     print "here?", e.PID, "fr ", fr
                    #     try:
                    #         print (e.PID in hosps[fr].patients)
                    #     except:
                    #         pass
                    #     sys.exit()

                    """
                    # HERE THE NEW LINES FOR THE "FORGETTING"(Scenario B) WERE ADDED
                    # if g[e.PID] == -1: # independend from the source hospital
                    #     g[e.PID] == 1

                    # the case of returning home (more realistic)
                    if g[e.PID] == -1 and to == -1 : # independend from the source hospital
                        g[e.PID] == 1
                    """

                    ### --------------- S C R E E N I N G ---------------
#                    print "av.screening_fract", av.screening_fract
#                    Screening(t, e, to, g, av.intervention_time, av.screening_fract, top_nodes, av.decontamination, Q, hosps)

                    if av.screening_flag == 1:
                        
                        Screening(t, e, to, g, av, top_nodes, Q, hosps)
                        
                    try:
                        hosps[to].into(g, e.PID)
                    except KeyError:
                        print to, "pid ", e.PID
                        sys.exit()



    outf.close()
