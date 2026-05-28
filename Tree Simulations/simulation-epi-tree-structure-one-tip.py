from random import random
from random import randrange
from random import expovariate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import bisect
import tqdm
import os


# seed(1234)
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# ##############################################################################################
#  Parameter class
# ##############################################################################################

class Para:
    # Keeps the parameters of the simulation
    def __init__(self):
        self.k       = 4                       # Number of downstream nodes
        self.fixed   = True                    # if True, degree is fixed; otherwise, it follows a Poisson distribution
        self.N       = 5000                    # Population size
        self.beta    = 1.5                     # Contact rate per individual
        self.mu      = 0.5                     # Spontaneous recovery
        self.sigma   = 0.5                     # Symptomatic/diagnosed recovery
        self.tagTime         = 0.6             # Time at which we try to tag an individual
        self.tagState        = (-1, self.k)    # Updated tagState
        self.timeHorizont    = 3.6             # Simulate until the last individual infected before that time did recover 
        self.tracingInterval = 1               # Tracing interval
        self.tracing         = False           # Tracing on/off
        
        # Specify different tracing possibilities
        self.forwardTracing     = False
        self.backwardTracing    = False
        self.clusterTracing     = False  # Not implemented

        self.rekursive          = True   # Recursive/one step
        self.p                  = 0.9    # Tracing probability

        # Derived parameters
        self.pobs = self.sigma / (self.sigma + self.mu)
        self.degree = 'fixed' if self.fixed else 'Poisson'


# ##############################################################################################
#  Event class
# ##############################################################################################


class Event:
    # store all parameters for an event
    def __init__(self, myType, idA, idB, eventTime):
        self.type       = myType     # type of event (infection or recovery)
        self.idA        = idA;       # upstream
        self.idB        = idB;       # downstream
        self.eventTime  = eventTime  # time of event
        self.deactivate = False      # in case of a recovery event; deactivated if individual is traced
    
    # def __lt__(self, other):
    #     """Ensure ordering by eventTime."""
    #     return self.eventTime < other.eventTime
    
    def printMe(self):
        print("type: ", self.type, "time = ", self.eventTime, " idA: ", self.idA, " idB :", self.idB);



# ###############################################################################################
#  Population class
# ###############################################################################################


class Population:

    def __init__(self, para_obj):
        
    ## init population object
        self.my_para          =  para_obj
        self.time             =  0   # initial time
        self.eventQueue       =  []  # no events in the beginning
        self.pop_size         =  0   # population count
        self.indiv_list       =  []  # individual stored in a list
        self.noInfecteds      =  0   # number of infecteds
        self.taggedInfected   = -1;  # the individual that has been tagged
        self.didNotFindIndiToTag = 0; # we put this var to 1 if we cannot tg an appropriate indi
        self.didDoTagIndi     = 0;  # we put this var to 1 if we did tag an appropriate indi
        
        self.TagInfectedTime = -1
        self.TagCreatedTime  = -1
        
        self.indiv_list.append(Individual(self.pop_size, -1, "no", self))  # create and store root individual in the individual class
        self.pop_size += 1          # increase population size by 1
        self.unobservedRecoveryEvent    =  []
        self.observedRecoveryEvent      =  []
        self.tmpTime                    =  0   # initial time until first observation in tagged clade is observed
        self.tmpPhyState                =  0   # type of observed individual
        
        self.indiv_list[0].HistoryTime.append(self.time)
               
        ## Lets infect the root at time t= 0 and generate downstreams
        self.indiv_list[0].infectMe(-1, False, self.time)

        ## produce a tag event
        tagEvent = Event("tagAnIndividual", para_obj.tagState[0], para_obj.tagState[1], para_obj.tagTime); # we want to tag some indi at time 0.3 why not?
        self.registerEvent(tagEvent);
        
        ## produce a time horizont event  
        thEvent = Event("timeHorizont", 0, 0, para_obj.timeHorizont); # we want to tag some indi at time 0.3 why not?
        self.registerEvent(thEvent);
        
    ## produce and return new individual
    def generateIndividual(self, infectorID, Tag, globalTime):
        newID = self.pop_size;
        myInd = Individual(newID, infectorID, Tag, self);
        myInd.HistoryTime.append(globalTime)
        self.indiv_list.append(myInd);
        self.pop_size += 1;
        return newID
             
    
    ## increase no of infected by 1 if we register a new infection
    def increaseInfecteds(self):
        self.noInfecteds += 1;
        
    ## decrease no of infected by 1 if we register a new infection
    def decreaseInfecteds(self):
        self.noInfecteds -= 1
    
    def registerEvent(self, eventIn):
        """Insert event into eventQueue while maintaining sorted order."""
        bisect.insort(self.eventQueue, eventIn, key=lambda x: x.eventTime)
  
    ## handle events
    def handleEvent(self):
        myEvent   = self.eventQueue.pop(0)
        self.time = myEvent.eventTime
        indA      = myEvent.idA
        indB      = myEvent.idB

        # if we have a tagged individualk: only hanlde this individual 
        # and its descendants
        if self.didDoTagIndi>0 and indA>0 and self.indiv_list[indA].IamTagged != True:
            return;
 

        if myEvent.deactivate:
            return
        
        ## contact event
        if myEvent.type == "contact":
            # self.infectionTimes.append(self.time)  #store infection time in infection times list
            if (self.indiv_list[indA].state == "I") and (self.indiv_list[indB].state == "S"):
                xTag  = self.indiv_list[indA].IamTagged; 
                self.indiv_list[indB].infectMe(indA, xTag, self.time);
                self.indiv_list[indA].registerContact(indB, self.time, self.time);
                return;
            # if (self.indiv_list[indB].state == "I") and (self.indiv_list[indA].state == "S"):
            #     xTag  = self.indiv_list[indB].IamTagged;
            #     self.indiv_list[indA].infectMe(indB, xTag, self.time);
            #     self.indiv_list[indB].registerContact(indA, self.time, self.time);
            #     return;
        
        # spontaneous recovery
        if myEvent.type == "SpntRcvr":
            self.indiv_list[indA].recoverUnobserved(self.time);
            self.indiv_list[indA].registerRecovery(self.time);
            self.unobservedRecoveryEvent.append(myEvent);
            return;
        
        # observed recovery
        if myEvent.type == "ObsRcvr":
            self.observed_type = self.indiv_list[indA].getType();
            self.indiv_list[indA].recoverObserved(self.time);
            self.indiv_list[indA].registerRecovery(self.time);
            if self.indiv_list[indA].IamTagged == True:
                self.tmpTime     =  self.time;
                self.tmpPhyState = self.indiv_list[indA].phyState;
            self.observedRecoveryEvent.append(myEvent);
            return;
        
        if myEvent.type == "timeHorizont":
            # time horizont reached
            self.time = self.time+0.0001; # increase a bit to ensure that the prog stops
            return;
        
        if myEvent.type == "tagAnIndividual":
            # find an ifected individual of type (4,4) and tag that
            myList=[];
            for i in range(self.pop_size):
                if self.indiv_list[i].state == "I":
                    myList.append(i);
            if len(myList) == 0:
                # did not find ANY indi of the approprate type
                self.didNotFindIndiToTag = 1;  # sign that no indi could be tagged
                return;
                
                self.tagState 
                
            # now select randomly an indi out of the myList
            ii = randrange(len(myList));
            idx  = myList[ii];
            self.my_para.tagState = self.indiv_list[idx].phyState;
            
            self.indiv_list[idx].IamTagged = True;
            self.taggedInfected = self.indiv_list[idx];
            self.didDoTagIndi   = 1;
            
            self.TagInfectedTime = self.indiv_list[idx].infectionEventTime
            self.TagCreatedTime  = self.indiv_list[idx].HistoryTime[0]
            #self.my_para.k = self.my_para.kPrime;

    def printEvents(self):
        for i in range(0, len(self.eventQueue)):
            self.eventQueue[i].printMe()
        print()

    def printMe(self):
        for i in range(0, self.pop_size):
            self.indiv_list[i].printMe()

    def printState(self):
        print(self.time, "\t", self.noInfecteds, end="")
        
    # def cladeStatistics(self):
    #     cladeInfInd = [i for i in range(self.pop_size) if self.indiv_list[i].IamTagged and self.indiv_list[i].state == 'I']
    #     cladeSusInd = [i for i in range(self.pop_size) if self.indiv_list[i].IamTagged and self.indiv_list[i].state == 'S']
    #     cladeRecInd = [i for i in range(self.pop_size) if self.indiv_list[i].IamTagged and self.indiv_list[i].state == 'R']
        
    #     print(len(cladeInfInd), len(cladeSusInd), len(cladeRecInd))


# #####################################################################################################
#  Individual class
# #####################################################################################################


class Individual:
    
    ## initialization
    def __init__(self, my_id_in, infectorID, IamTagged, pop_obj_in):
        self.my_id       = my_id_in            # own id
        self.infectorID  = infectorID          # infector ID
        self.pop_obj     = pop_obj_in          # my population object
        self.paraObj     = pop_obj_in.my_para  # my population para object
        self.state       = "S"                 # disease state of an individual
        self.contactees      = []              # infected downstreams
        self.contactTimes    = []              # infected times of infected downstreams
        self.HistoryTime     = []
        self.IamPrimaryIndex = False  # initial primary index case is set to false
        self.myRecoveryEvent = 0      # empty event
        self.myDownstreams   = []     # all downstreams
        self.myDegree        = self.paraObj.k if self.paraObj.fixed else np.random.poisson(self.paraObj.k);
        self.j               = 0           # number of taken infected downstreams. 0 when infected
        self.phyState        = [self.j,self.myDegree]  # state of an individual 
        self.IamTagged       = IamTagged   # identity of an individual
        self.myRecoverType   = "x"         # 1 for observed recovery, 0 for unobserved recovery
        self.tmpInfectTime   = 0           # time to be infectious
        self.infectionEventTime = -1

    ## infect an individual
    def infectMe(self, infectorID, xTag, globalTime):
        
        # generate downstream individuals
        for ii in range(0, self.myDegree):
            #print("degree: ", self.myDegree)
            ind       = self.pop_obj.generateIndividual(self.my_id, xTag, globalTime);
            contactee = ind;
            self.myDownstreams.append(contactee);
            eventTime = globalTime+expovariate(self.paraObj.beta);    ## contact per edge
            myEvent   = Event("contact", self.my_id, contactee, eventTime);
            self.pop_obj.registerEvent(myEvent);
        
        # set my own state - Infection, update IamTagged.
        self.state     = "I"
        self.IamTagged = xTag;
        self.infectionEventTime = globalTime
 
        self.pop_obj.increaseInfecteds() 
        self.HistoryTime.append(globalTime)
        # self.j        = 0;
        # self.phyState = [0,self.paraObj.k];
        
        # define recovery event
        recoverTime = globalTime+expovariate(self.paraObj.mu + self.paraObj.sigma)
        self.tmpInfectTime  = recoverTime - globalTime
 
        if random() < self.paraObj.pobs:
            myEvent = Event("SpntRcvr", self.my_id, -1, recoverTime);
        else:
            myEvent = Event("ObsRcvr", self.my_id, -1, recoverTime);
        self.pop_obj.registerEvent(myEvent);
        self.myRecoveryEvent = myEvent;
        return;   

    ## register contact: I did infect a downstream individual
    def registerContact(self, indi2, contactTime, globalTime):
        self.contactees.append(indi2);
        self.contactTimes.append(contactTime);
        self.j        = self.j + 1;
        self.phyState = [self.j, self.myDegree];

    ## register recovery   
    def registerRecovery(self, globalTime):
        self.recoveryTime = globalTime;

    ## register unobserved event
    def recoverUnobserved(self, globalTime):
        self.myRecoverType = 0
        self.state         = "R"
        self.pop_obj.decreaseInfecteds()
        self.HistoryTime.append(globalTime)

    ## register observed event
    def recoverObserved(self, globalTime):
        self.myRecoverType = 1
        self.state = "R"
        self.pop_obj.decreaseInfecteds() 
        self.HistoryTime.append(globalTime)
    
    def getType(self):
        # return (i,k), (# no infected downstreams, degree)
        return(self.phyState);
        
    def printMe(self):
        print("my ID:  ", self.my_id, "state ",self.state)
#%%
# Press the green button to run the script.
if __name__ == '__main__':
    
    j_obs_list        = []
    k_obs_list        = []
    cladeLifeSpanList = []
    i_tag_list        = []
    k_tag_list        = []
    
    # taggedIndiInfectTimeList  = []      # infection time of tagged individuals
    # taggedIndiCreatedTimeList = []      # creation time of tagged individuals
    taggedTimeList            = []

    N = 500000
    for ii in tqdm.tqdm(range(0, N)):
        
    
        para = Para()
        popu = Population(para);

        ttt = 0;
        OK = True
        nextTime = 1.0;
        span  = -1.234; # mke sure that we are initialized
        while OK: 
            popu.handleEvent();
            
            # we need to decide if the simul is done.
            simulIsDone = False;
            if len(popu.eventQueue)==0:
                #print("empty queue!");
                simulIsDone = True;
            if popu.noInfecteds ==0:
                simulIsDone = True;
            if popu.time > para.timeHorizont:
                simulIsDone = True;
            if popu.didNotFindIndiToTag > 0:
                simulIsDone = True;
            if popu.tmpTime > 0:
                simulIsDone = True;

            span = -999;
            # if we are done, get the result:
            if simulIsDone==True:         # Invalid run    
                if popu.didDoTagIndi == 0:
                    span         = -100
                    j_obs,k_obs  = -5,-5
                    i_tag, k_tag = -5,-5
                    
                    taggedIndiInfect      = -200
                    taggedIndiCreated     = -400
                    taggedTime            = -500
                    break

                elif popu.tmpTime > 0: # Observed
                    span               = popu.tmpTime
                    j_obs,k_obs        = popu.tmpPhyState[0], popu.tmpPhyState[1]
                    taggedIndiInfect   = popu.TagInfectedTime
                    taggedIndiCreated  = popu.TagCreatedTime   
                    taggedTime         = para.tagTime
                    i_tag, k_tag       = para.tagState[0], para.tagState[1]
                    break

                else: # Extinct
                    span  = float('inf')
                    j_obs,k_obs = None,None
                    taggedIndiInfect      = 202
                    taggedIndiCreated     = 302  
                    taggedTime            = para.tagTime
                    i_tag, k_tag          = para.tagState[0], para.tagState[1]
                    break

            if span>-200:
               break;                 
        cladeLifeSpanList.append(span)
        j_obs_list.append(j_obs)
        k_obs_list.append(k_obs)
        i_tag_list.append(i_tag)
        k_tag_list.append(k_tag)
        degree = para.degree
        
        # taggedIndiInfectTimeList.append(taggedIndiInfect)
        # taggedIndiCreatedTimeList.append(taggedIndiCreated)
        
        taggedTimeList.append(taggedTime)

        # if span >-99:
        #     print(span," ", (j_obs, k_obs));
        #     # print(span);
#%%
df = pd.DataFrame({
    'i_foc': i_tag_list,
    'k_foc': k_tag_list,
    't_foc': taggedTimeList,
    'j_obs': j_obs_list,
    'k_obs': k_obs_list,
    't_obs': cladeLifeSpanList

})

df = df[df['t_obs'] != -100]
df.to_csv(f'phylo-epi-sim-data-{degree}.csv', index = False)