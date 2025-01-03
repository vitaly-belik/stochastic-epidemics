# utilities for the sis.py

"""
    7.09.17 flag_scr_immer was introduced and sis.py and utils.py were rewritten
    
    22.04.2017 Now the parameter file is copied into results
    
    26.09.2016 The time of recovery is stored in "guy" object. The infected state
    corresponds to "guy.state" field > zero. This is needed to implement
    decontamination

    TOP_nodes are included, where TOP_nodes are read from a file. To think about
    implementing top_links.

    Parameter file now includes the variables names

"""



import re
import string
import random
import math

class guys(object):

    def __init__(self, state, loc):
        self.state = 0
        self.loc = loc


class event(object):

    def __init__(self, t, PID, FROM, TO, type):
        self.t = t
        self.PID = PID
        self.FROM = FROM
        self.TO = TO
        self.type = type
        # KIND = {"I","R","T"} type of event


    def __cmp__(self, other):
        return cmp(self.t, other.t)


def get_recovery_time(beta):
    # prototype
    xi = random.random()
    return -1.0 / (beta) * math.log(xi)


def get_recovery_time2(beta, d = 0.5):
    # prototype
    xi = random.random()
    if random.random() < d:
        return -1.0 / (beta) * math.log(xi)
    else:
        return -1.0 / (beta) * math.log(xi) / 2.0 # reducing the recovery time


def get_infection_time(alpha, I, S, N):
    xi = random.random()
    return -1.0 / (alpha / float(N) * S * I) * math.log(xi)


#alpha = 0.03



class hospital(object):

    __slots__ = ('last_event', 'patients', 'HID', 'sus', 'inf', 'rec')

    def __init__(self, last_event, patients, HID, g):
        # g is an dict with values of type "guy" describing patients

        self.patients = patients    # list of patients' IDs
        self.HID = HID
        self.sus = set()
        self.inf = set()
        self.rec = set()

        self.last_event = last_event

        for el in patients:
            if g[el].state > 0:
                self.inf.add(el)

            if g[el].state == 0:
                self.sus.add(el)

            if g[el].state == -1:
                self.rec.add(el)

    def out(self, g, PID):

        self.patients.remove(PID)

        if g[PID].state == 0:
            self.sus.remove(PID)
        if g[PID].state > 0 :
            self.inf.remove(PID)
        if g[PID].state == -1:
            self.rec.remove(PID)

    def into(self, g, PID):

        self.patients.add(PID)

        if g[PID].state == 0:
            self.sus.add(PID)

        if g[PID].state > 0:
            self.inf.add(PID)

        if g[PID].state == -1:
            self.rec.add(PID)

        g[PID].loc = self.HID # location is always set in "into"



    # in "update" infections happens and the corresponding recovery events are assigned
    def update(self, g, e, Q, av):

        if self.last_event < e.t:


            if len(self.inf) > 0 and len(self.sus) > 0:  # infection
                # N = len(self.inf) + len(self.sus) + len(self.rec)
                N = len(self.inf) + len(self.sus)

                it = self.last_event + get_infection_time(av.alpha, len(self.inf), len(self.sus), N)
            else:
                self.last_event = e.t
                return 1

            while it < e.t:  # we can infect
                self.last_event = it
                
                # infection
                # here we choose randomly a susceptible individual
                # one could do this according to some coupling matrix
                # based on some (vector of properties of agents).
                # To think using IDs of patients as vector labels?
                
                some = random.sample(self.sus,1)[0]

                # g[some].state = 1 # 26.09.2016


                # output ------------------------------------------------------
                """
                print it, len(self.inf), len(self.sus),
                print len(self.inf) + len(self.sus),
                print len(self.patients), self.HID
                """
                rt = it + get_recovery_time(av.beta)
                g[some].state = rt # < --------------- 26.09.2016

                self.sus.remove(some)
                self.inf.add(some)


                re = event(rt, some, -1, -1, 'R') # a new recovery event
                Q.put( (re.t, re) )  # putting the NEW event into the queue

                if rt < e.t:
                    return 0

                if len(self.inf) > 0 and len(self.sus) > 0:
                    # N = len(self.inf) + len(self.sus) + len(self.rec)
                    N = len(self.inf) + len(self.sus) # the recovered or screened
                    it = it + get_infection_time(
                        av.alpha,
                        len(self.inf),
                        len(self.sus),
                        N)
                else:
                    self.last_event = e.t
                    return 1
                    # break

            if it > e.t:
                self.last_event = e.t
                # it means there are no infection in before the next "event"
                # print "infection after event",self.HID
                return 1
        else:
            self.last_event = e.t
            return 1


########################
def ff(lines, i):
    return re.split('\s', string.strip(lines[i]))[-1]


class args(object):

    __slots__ = ('alpha',\
        # = float(sys.argv[1])  # infection rate
        'beta',\
        # recovery rate
        'max_t',\
        # maximal time of observation
        'file',\
        # file with links (+ artificial patients)
        'seed',\
        #
        'NOI',\
        # or percentage of patients, initially infected
        'start',\
        # time of the start of infection
        'scr_hosp',\
        # the screened hosp!!!!
        # to exclude screening just set scr_hosp to 0!
        'hids_file',\
        # file with the HIDs
        'pids_file',\
        # file with the PIDs
        'out_file',\
        'screening_fract',\
        'top_nodes_file',\
        'screening_flag',\
        'decontamination',\
        'intervention_time',\
        'flag_scr_immer')

    def __init__(self, sav):

############### ----------------- to get rid of line arguments!!!-----------

        if len(sav) == 13:

            pass
        else:
            # read parameter file:

            import string
            import re
            import shutil


            f = open(sav[1], 'r')
            
            lines = tuple(f)

            #alpha = float(re.split('\s', string.strip(lines[0]))[0])  # infection rate

            self.alpha = float(ff(lines, 0))  # infection rate

            self.beta = float(ff(lines, 1))  # recovery rate
            self.max_t = int(ff(lines, 2))  # maximal time of observation
            self.file = ff(lines, 3)  # file with links (+ artificial patients)
            self.seed = int(ff(lines, 4))
            self.NOI = float(ff(lines, 5))  # or percentage of patients
            self.start = int(ff(lines, 6))  # time of the start of infection
            #init_hosp = int(string.strip(lines[7]))
            self.scr_hosp = int(ff(lines, 7))
            self.hids_file = ff(lines, 8) # file with the HIDs
            self.pids_file = ff(lines, 9) # file with the PIDs

            self.top_nodes_file = ff(lines, 10)
            self.screening_flag = int(ff(lines, 11))
            self.decontamination = float(ff(lines, 12)) # gives us the reduction of time until
            self.intervention_time = int(ff(lines, 13))# recovery
            self.screening_fract = float(ff(lines, 14))
            self.flag_scr_immer = int(ff(lines, 15)) # flag for temporary (until recovery)
                                                        # or permanent isolation


            self.out_file = sav[2]
            f.close()
            
            
            f = open(sav[1], 'r')
            
            f1 = open(sav[2] + '.txt', "w")
            shutil.copyfileobj(f, f1) # <- here the content of the parameters file is copied
            # to the outputfile
#            f1.write("bla\n")

            f.close()
            f1.close()
            # self.screening_fract = float(sav[3])


#        if argv[1] == "--h" or argv[1] == "-help":
#            print "input parameters are:"
#            print "alpha, beta, max_t, file, seed, NOI, start, init_hosp, hids_file, pids_file"
#            sys.exit(0)
#



def InitQueue(Q, f, chunck_size, inf_prob):


    q = 0

    # reading chunck_size number of lines for initialization
    for i in range(chunck_size):

        s = f.readline()
        if len(s) > 0:

            a = re.split('\s', string.strip(s))

            source = int(a[0])
            target = int(a[1])
            event_time = float(a[2])  # extremely important
            PID = int(a[3])


            if source == -1: # first appearence of the agent in the system
                if random.random() < inf_prob: # or out of the community
                    e = event(event_time, PID, -3, target, 'IN')
                    Q.put( (e.t, e) )

                    """
                    #here we infect all new coming patients

                    rt = get_recovery_time(beta) + e.t
                    #rt = get_recovery_time2(beta) + start
                    e = event(rt, e.PID, -1, -1, 'R')
                    Q.put( (e.t, e) )
                    """

                else:
                    e = event(event_time, PID, -1, target, 'IN')
                    Q.put( (e.t, e) )


            if target == -1:
                e = event(event_time, PID, source, -1, 'OUT')

            if target != -1 and source != -1:
                e = event(event_time, PID, source, target, 'T')

            Q.put( (e.t, e) )

    print "INIT????"
    e = event(event_time + 0.01,-2,-2,-2,'NEXT_CHUNCK')

    print "INIT????"
    # 0.01 is very important to sort properly
    Q.put( (e.t, e) )

    return Q



def ReadNextChunck(Q, f, chunck_size, inf_prob, myeof, NOI):

    for i in range(chunck_size):
        try:
            a = re.split(' ', string.strip(f. readline()))
        except IOError:
            myeof = 1
            break

        if a[0] != '':
            source = int(a[0])

            target = int(a[1])
            event_time = float(a[2])  # "float" is extremely important
            PID = int(a[3])


            if source == -1:

                # infection_probability accounts for some prevalence in
                # the community

                if random.random() < inf_prob: # or out of the community
                    e = event(event_time, PID, -3, target, 'IN')
                    Q.put( (e.t, e) )

                    """
                    # infecting new incoming patients
                    rt = get_recovery_time(beta) + e.t
                    #rt = get_recovery_time2(beta) + start
                    e = event(rt, e.PID, -1, -1, 'R')
                    Q.put( (e.t, e) )
                    """

                else:
                    e = event(event_time, PID, -1, target, 'IN')
                    Q.put( (e.t, e) )




            if target == -1:
                e = event(event_time, PID, source, -1, 'OUT')
                Q.put( (e.t, e) )

            if target != -1 and source != -1:
                e = event(event_time, PID, source, target, 'T')
                Q.put( (e.t, e) )


    #                    print ";;;;"
    #                    e = event(event_time, PID, source, target, 'T')

        else:
            myeof = 1

    if myeof == 0:
        e = event(event_time + 0.01,-1,-1,-1,'NEXT_CHUNCK')
        Q.put( (e.t,e) )

    return Q


def ProduceOutput(outf, hosps, t):

    outf.write(str(t))
    for h in sorted(hosps.keys()):
        s = " " + str(len(hosps[h].inf)) + " " + str(len(hosps[h].sus)) + " " + str(len(hosps[h].rec))
        #s = " " + str(len(hosps[h].inf))
        outf.write(s)
    outf.write("\n")


def InitialInfection(Q, e, NOI, beta, start, hosps, g ):

    fr = e.FROM
    # how many to infect:
    J = NOI * len(hosps[fr].sus) # NOI is a fraction
    """
    print "J = ", J, "S = ", len(hosps[fr].sus), "t = ", e.t
    """

    # to assure the fraction NOI on average
    if J < 1 and J > 0:
        if random.random() < NOI:
            J = 1
        else:
            J = 0

    if J >= 1:

        #for el in random.sample(hosps[fr].sus, int(math.ceil(NOI * len(hosps[fr].sus)))):
        try:
            smpl = random.sample(hosps[fr].sus, int(math.floor(J))) # sample of susceptibles in a hospital

        except ValueError:
            print len(hosps[fr].sus), int(math.floor(J)), J
            sys.exit()
        except TypeError:
            print "hosps[fr].sus", hosps[fr].sus, "J", J
            sys.exit()

        for el in smpl:
            # infection
            # print "I"
             #?????????????????
            hosps[fr].sus.remove(el)
            hosps[fr].inf.add(el)

            # ASSUMPTION OF INTERVENTION AFTER THE INITIAL INFECTION
            rt = get_recovery_time(beta) + start
            #rt = get_recovery_time2(beta) + start

            g[el].state = rt

            # HERE WE ADD NEW EVENTS
            e = event(rt, el, -1, -1, 'R')

            Q.put( (e.t,e) )


#def Screening(t, e, to, g, intervention_time, screening_fract, top_nodes, decontamination, Q, hosps):
def Screening(t, e, to, g, av, top_nodes, Q, hosps):
    
    # what about args as ckass?


    ### --------------- S C R E E N I N G ---------------
    # 26.09.2016 Alongside with screening, decontamination is also implemented


    # modifying the state of the screened patients from I to R

    #print "Travel"
    #if g[e.PID].state == 1:
    #    print "infectious travel"

    # we infect here all the incoming guys
    # if g[e.PID].state == 0 and random.random() < 0.02:
    #    g[e.PID].state = 1
    # every time some one arrives into the hospital

    #if to != 0:
    #    print e.t, to, g[e.PID].state
    # here we check for the state of the incoming guys and
    # calculate their fraction. The question is over which period
    # of time is the average fraction defined?



    if t > av.intervention_time:
#        if e.PID == 213024:
#            print "intervention time", av.intervention_time

        if to != 0 and to in top_nodes: # here we could restrict screening to top_nodes
#        to = -1>
            if g[e.PID].state > 0: # if the incoming patient is infected

                if random.random() < av.screening_fract:
                
                    # here isolated either until recovery or forever
                    # this is decided in "sis.py" in the recovery part
                    
                    
                    
#                    hosps[to].inf.remove(e.PID)
#                    hosps[to].rec.add(e.PID)

                    if av.decontamination > 0.0: # scenario "C" (7.9.2017)
                    
                        if random.random() < av.decontamination: # reduce recovery time

#                            new_t = t + (g[e.PID].state - t)/3.0 # to think if this makes sense?
                            # could it happen that we already made "state = -1" ?
                            
                            new_t = t + (g[e.PID].state - t)/3.0 # to think if this makes sense?
            
                            new_event = event(new_t, e.PID, -1, -1, 'R')
                            Q.put( (new_event.t, new_event) )
                    else:
                        pass
                    
                    g[e.PID].state = -1 # Important, that this is after decontamination
                    
                    
                        
                else:
                    pass
