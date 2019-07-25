
# coding: utf-8

# In[ ]:


#simulate the model of epidemics with different incubation times

    #a) agent-based model

        #agents are in four compartments: SEIR
        #S are susceptible, E are infected but in the latency state, I are infected and infectious and R recovered.
        #S+E+R move freely on a spatial network (2d lattice)
        #I do not move
        
        #S+I->I+E with prob beta
        #E->I after a certain latent time T_l (deterministic or from a distribution)
        #I->R after a certain infectious time T_i (deterministic or from a distribution)
        
    #b) network minimal model for comparison of basic mechanism
        
        #metapopulation-like model
        #interaction network is a sum of powers of the mobility network in a) up to power T_l
        #Time is in generations, where a genera1tion is takes a time of T_l
        
import networkx as nx
import numpy as np
from numpy.random import binomial,seed
import random
import matplotlib.pyplot as plt
import os
from scipy.stats import erlang
import csv
import pandas as pd

random.seed(12345)
seed(seed=12345)
#a) ABM

#initialize

## space and demographics

l=15 # side of the rectangular lattice
l2=10 # side of rectangular lattice
n_agloc=300 # number of agents in each location initially

N_loc=l*l2 # number of locations
N_agents=n_agloc*l*l2# total number of agents
G=nx.grid_2d_graph(l,l2) # network for the spatial substrate

## epidemic model parameters

m=1 # number of steps per day

alpha=0.1
beta=0.3
mu_recover=0.2
delta=0.001 #flight prob


exp=2 #exponent for the gravity law for mobility



T_l=1 # incubation period
alpha_step=1-np.power((1-alpha),1.0/float(m)) # probability of moving
beta_step=beta/(m) # probability of disease transmission
mu_step=1-np.power((1-mu_recover),1.0/float(m)) # probability of recovery
gamma=0 # increase/decrease on moving probability for infectious individuals

# other parameters

end_time=500

#choose how the incubation periods are chosen from an Erlang distribution
# parameters of the distribution: k=shape, mu=scale; 
                                #T_l=k*mu (mean), sigma^2=k*mu^2 (variance);
                                #k=(T_l/sigma)^2, mu=sigma^2/T_l
    #1- constant variance sigma_0^2 across different T_l's
    #2- constant k (not 1, because with k=1 we have an exponential distribution),
    #sigma increases with T_l 

erlang_type=2
sigma_0=2 # will use only in case erlang_type=1 (constant variance)
k=10 # will use only in case erlang_type=2 (constant k)
if erlang_type==1:
    a=T_l/sigma_0
    mu=sigma_0/a
    k=a*a
elif erlang_type==2:
    mu=T_l/k
    
rv=erlang(k,scale=mu)

## give locations to all agents

Nruns=10


peak_height=list()
peak_time=list()
start_time=list()
inf_per_day=list()
day_50=list()
day_100=list()
day_150=list()
day_200=list()


for irun in range(Nruns):
    print (Nruns-irun,'start')

    loc=dict() # location of each agent
    agents=dict() # agents in each location
    S=dict() # S agents in each location
    E=dict() # E agents in each location
    I=dict() # I agents in each nrunslocation
    R=dict() # R agents in each location
    day_count=dict()



    peak_height.append(np.zeros((l,l2)))
    peak_time.append(np.zeros((l,l2)))
    start_time.append(np.zeros((l,l2)))
    day_50.append(np.zeros((l,l2)))
    day_100.append(np.zeros((l,l2)))
    day_150.append(np.zeros((l,l2)))
    day_200.append(np.zeros((l,l2)))
    time_series=np.zeros((end_time+2,l*l2))
    active_agents=np.zeros((end_time+2,l*l2))


    k=0
    for i in range(l):
        for j in range(l2):
            iloc=(i,j)
            agents[iloc]=set()
            S[iloc]=set()
            E[iloc]=set()
            I[iloc]=set()
            R[iloc]=set()
            day_count[iloc]=set()
            for j in range(n_agloc):
                loc[k]=iloc
                agents[iloc].add(k)
                S[iloc].add(k)
                k+=1
    
    pop_size=np.zeros((l,l2))
    for x in range(l):
        for y in range(l2):
            pop_size[x][y]=n_agloc
    
    mob_net = pd.read_csv('SL_mobnet_equalpop.csv', sep=',',header=0)
    mob_net=mob_net.values

            
    state=np.zeros(N_agents) # epidemic state of each agent
    inf_time=np.zeros(N_agents) # time at which agent became E
    ## initial condition for the epidemics

    iloc=(int(l/2-0.5),int(l2/2)) # the place where we start having a percentage perc of E individuals 
    perc=0.05 # percentage of E agents initially in location iloc; the rest of agents are S


    for agent in random.sample(agents[iloc],int(perc*(n_agloc))):
        state[agent]=1
        tt=rv.rvs(1)
        #tt=np.random.uniform(0,10)
        #print (tt)
        inf_time[agent]=tt
        #inf_time[agent]=T_l # change here the time T_l for one from a distribution to get variation on it
        S[iloc].remove(agent)
        E[iloc].add(agent)
        
    ###SAVE INITIAL CONDITION/PRINT


    #data structure to save all epidemic curves

    epi_curves=dict()
    #for ix in range(l):
    #    epi_curves[ix]=dict()
     #   for iy in range(l):
      #      epi_curves[ix][iy]=list()

    i_plot=np.zeros((l,l2))
    #for ix in range(l):
     #   for iy in range(l):
      #      kk=float(len(I[(ix,iy)]))/float(len(agents[(ix,iy)]))
       #     i_plot[ix][iy]=kk
        #    epi_curves[ix][iy].append(kk)
    #fig=plt.figure()
    ##plt.subplot(221,title='$P$ only space')
    #plt.imshow(i_plot,vmin=0, vmax=1, cmap='jet')
    #fig.savefig('a_%.3i.png' % 0,bbox_inches='tight')
    ##plt.show()
    #plt.close()



    I_tot=0
    #print(0,I_tot)
    ## start the epidemics

    locations=G.nodes()
    for time in range(end_time):
        for k in range(m):
            # first move
            for iloc in locations:
                #print(time,iloc,len(agents[iloc]))
                fro=iloc[0]+l*iloc[1] 
                #print(iloc,fro)
                N_move=binomial(len(agents[iloc])-len(I[iloc]),alpha) # number of agents moving
                agents_move=random.sample(agents[iloc]-I[iloc],N_move) # set of agents moving
                #N_flight=binomial(len(agents_move),delta)
                #agents_flight=random.sample(agents_move,N_flight)
                #agents_flight=set(agents_flight)
                agents_move=set(agents_move)
                #agents_diff=set()
                #agents_diff=agents_move-agents_flight
                # places where they will move
                #p=mob_net[fro][:]
                #print(p)
                dest=np.random.choice(np.arange(l*l2),size=N_move,replace=True,p=mob_net[fro][:])
                i=0
                for agent in agents_move:
                    to=dest[i]
                    xdest=to-l*int(to/l)
                    ydest=int(to/l)
                    #jloc=random.sample(G.neighbors(iloc),1)[0]
                    jloc=(xdest,ydest)
                    loc[agent]=jloc
                    agents[iloc].remove(agent)
                    agents[jloc].add(agent)
                    if state[agent]==0:
                        S[iloc].remove(agent)
                        S[jloc].add(agent)
                    elif state[agent]==1:
                        E[iloc].remove(agent)
                        E[jloc].add(agent)
                    else:
                        R[iloc].remove(agent)
                        R[jloc].add(agent)
                    i+=1
                #for agent in agents_flight:
                    #jloc=random.choice(locations)
                    #loc[agent]=jloc
                    #agents[iloc].remove(agent)
                    #agents[jloc].add(agent)
                    #if state[agent]==0:
                        #S[iloc].remove(agent)
                        #S[jloc].add(agent)
                    #elif state[agent]==1:
                        #E[iloc].remove(agent)
                        #E[jloc].add(agent)
                    #else:
                        #R[iloc].remove(agent)
                        #R[jloc].add(agent)
                ###HAVE TO INCLUDE THE TELEPORTATIONS FOR I IN CASE GAMMA!=0
                N_move=binomial(len(I[iloc]),alpha*gamma)
                dest=np.random.choice(np.arange(l*l2),size=N_move,replace=True,p=mob_net[fro][:])
                i=0
                for agent in random.sample(I[iloc],N_move):
                    #jloc=random.sample(G.neighbors(iloc),1)[0]
                    to=dest[i]
                    xdest=to-l*int(to/l)
                    ydest=int(to/l)
                    jloc=(xdest,ydest)
                    loc[agent]=jloc
                    agents[iloc].remove(agent)
                    agents[jloc].add(agent)
                    I[iloc].remove(agent)
                    I[jloc].add(agent)
            # second infection dynamics
            col=0
            for iloc in locations:
                beta_step=beta/len(agents[iloc])
                #print(iloc)
                N_inf=binomial(len(S[iloc]),1.0-np.power(1.0-beta_step,len(I[iloc])))
                N_rem=binomial(len(I[iloc]),mu_step)
                for agent in random.sample(S[iloc],N_inf):
                    state[agent]=1
                    tt=rv.rvs(1)
                    #tt=np.random.uniform(0,10)
                    inf_time[agent]=time+tt
                    #inf_time[agent]=time+T_l # change here the time T_l for one from a distribution to get variation on it
                    S[iloc].remove(agent)
                    E[iloc].add(agent)
                    day_count[iloc].add(agent)
                for agent in random.sample(I[iloc],N_rem):
                    state[agent]=3
                    I[iloc].remove(agent)
                    R[iloc].add(agent)
                a=set(E[iloc])
                New_inf=0
                for agent in a:
                    if inf_time[agent]<=time:
                        state[agent]=2
                        E[iloc].remove(agent)
                        I[iloc].add(agent)
                        New_inf+=1
                time_series[0][col]=iloc[0]
                time_series[1][col]=iloc[1]
                time_series[time+2][col]=New_inf
                active_agents[0][col]=iloc[0]
                active_agents[1][col]=iloc[1]
                active_agents[time+2][col]=len(I[iloc])
                col+=1
        I_tot=0
        for iloc in locations:
            I_tot+=len(I[iloc])
        inf_per_day.append(I_tot)
        for ix in range(l):
            for iy in range(l2):
                if float(len(agents[(ix,iy)]))>0:
                    kk=float(len(I[(ix,iy)]))/float(len(agents[(ix,iy)]))
                else:
                    kk=0 
                i_plot[ix][iy]=kk
                #epi_curves[ix][iy].append(kk)
                if kk > peak_height[irun][ix][iy]:
                    peak_height[irun][ix][iy]=kk
                    peak_time[irun][ix][iy]=time+1
                if start_time[irun][ix][iy]==0:
                    if len(I[(ix,iy)])>0:
                        start_time[irun][ix][iy]=time+1
                if time==50:
                    day_50[irun][ix][iy]=len(day_count[ix,iy])
                if time==100:
                    day_100[irun][ix][iy]=len(day_count[ix,iy])
                if time==150:
                    day_150[irun][ix][iy]=len(day_count[ix,iy])
                if time==200:
                    day_200[irun][ix][iy]=len(day_count[ix,iy])

            fig=plt.figure()
        #plt.subplot(221,title='$P$ only space')
        plt.imshow(i_plot,vmin=0, vmax=0.5, cmap='jet')
        plt.colorbar()
        fig.savefig('a_%.3i.png' % (time+1),bbox_inches='tight')
        #plt.show()
        plt.close()
    np.savetxt('inf_per_day'+str(T_l)+'_'+str(irun)+'.csv', inf_per_day)
    np.savetxt('time_series'+str(T_l)+'_'+str(irun)+'.csv', time_series,delimiter=',')
    np.savetxt('active_agents'+str(T_l)+'_'+str(irun)+'.csv', active_agents,delimiter=',')


    os.system('mencoder mf://a_*.png -mf w=800:h=600:fps=10:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o latent_times_l_'+str(l)+'_Tl_'+str(T_l)+'_'+str(m)+'_'+str(irun)+'.avi')

    os.system('rm a_*.png')
    print(Nruns-irun,'end')


    

av_peak_height=np.zeros((l,l2))

for ix in range(l):
    for iy in range(l2):
        for irun in range(Nruns):
            av_peak_height[ix][iy]+=peak_height[irun][ix][iy]/float(Nruns)

fig=plt.figure()
#plt.subplot(221,title='$P$ only space')
plt.imshow(av_peak_height,vmin=0, vmax=1, cmap='jet')
plt.colorbar()
fig.savefig('av_peak_height_l_'+str(l)+'_Tl_'+str(T_l)+'_'+str(m)+'.png',bbox_inches='tight')
np.savetxt('av_peak_height'+str(T_l)+'.csv', av_peak_height, delimiter=",")
#plt.show()
plt.close()

#figure av_peak_time 2d


av_peak_time=np.zeros((l,l2))
av_start_time=np.zeros((l,l2))
av_day_50=np.zeros((l,l2))
av_day_100=np.zeros((l,l2))
av_day_150=np.zeros((l,l2))
av_day_200=np.zeros((l,l2))


for ix in range(l):
    for iy in range(l2):
        for irun in range(Nruns):
            av_peak_time[ix][iy]+=peak_time[irun][ix][iy]/float(Nruns)
            av_start_time[ix][iy]+=start_time[irun][ix][iy]/float(Nruns)
            av_day_50[ix][iy]+=day_50[irun][ix][iy]/float(Nruns)
            av_day_100[ix][iy]+=day_100[irun][ix][iy]/float(Nruns)
            av_day_150[ix][iy]+=day_150[irun][ix][iy]/float(Nruns)
            av_day_200[ix][iy]+=day_200[irun][ix][iy]/float(Nruns)

fig=plt.figure()
#plt.subplot(221,title='$P$ only space')
plt.imshow(av_peak_time,vmin=0,vmax=end_time,cmap='jet')
plt.colorbar()
fig.savefig('av_peak_time_l_'+str(l)+'_Tl_'+str(T_l)+'_'+str(m)+'.png',bbox_inches='tight')
np.savetxt('av_peak_time'+str(T_l)+'.csv', av_peak_time, delimiter=",")
np.savetxt('av_start_time'+str(T_l)+'.csv', av_start_time, delimiter=",") 
#np.savetxt('av_day_50_'+str(T_l)+'.csv', av_day_50, delimiter=",") 
#np.savetxt('av_day_100_'+str(T_l)+'.csv', av_day_100, delimiter=",") 
#np.savetxt('av_day_150_'+str(T_l)+'.csv', av_day_150, delimiter=",") 
#np.savetxt('av_day_200_'+str(T_l)+'.csv', av_day_200, delimiter=",") 
#plt.show()
plt.close()

#figure av_peak_time as a function of distance

maxdist=int(3.0*l/4.0)

av_peak_time=np.zeros((maxdist)) #!!!!  BE CAREFUL
av2_peak_time=np.zeros((maxdist))
norm_av_peak_time=np.zeros((maxdist))
a=int(l/2)
for ix in range(l):
    for iy in range(l2):
        d=abs(ix-a)+abs(iy-a)
        if d<maxdist:
            #print(d)
            for irun in range(Nruns):
                kk=float(peak_time[irun][ix][iy])
                av_peak_time[d]+=kk
                av2_peak_time[d]+=kk*kk
                norm_av_peak_time[d]+=1.0

for i in range(maxdist):
    #print(av_peak_time[i],av2_peak_time[i])
    av_peak_time[i]=av_peak_time[i]/norm_av_peak_time[i]
    av2_peak_time[i]=av2_peak_time[i]/norm_av_peak_time[i]
    av2_peak_time[i]=np.sqrt(abs(av2_peak_time[i]-av_peak_time[i]*av_peak_time[i]))

fig=plt.figure()
#plt.subplot(221,title='$P$ only space')
plt.ylim(0,400)
plt.errorbar(np.arange(maxdist),av_peak_time,yerr=av2_peak_time)
fig.savefig('av_peak_time_d_l_'+str(l)+'_Tl_'+str(T_l)+'_'+str(m)+'.png',bbox_inches='tight')
#plt.show()
plt.close()






