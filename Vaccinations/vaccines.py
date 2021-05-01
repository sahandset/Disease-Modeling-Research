import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pylab as plt
import random as rnd
rnd.seed()
import copy

class SIR():
    
    def drawGz(self,G,z,t):
        # DO NOT MODIFY THIS FUNCTION
        # This function draws G with node labels from partition z
        #
        # input  : G is a networkx graph
        #        : z is a dictionary of group labels for G's nodes
        # output : none
        # 
        # WARNING: function is optimistic: assumes inputs are properly formatted

        colors = ['#d61111','#11d646','#11c6d6','#d67711','#1b11d6','#d611cc'] # map node labels to colors (for the visualization)

        node_colors = []
        for i in G.nodes():
            if z[i]=='S':
                cid = 0
            if z[i]=='I':
                cid = 1
            if z[i]=='R':
                cid = 2
            if z[i]=='V':
                cid = 3
            node_colors.append(colors[int(cid)])
        nsize  = 600
        flabel = True

        if G.order() > 50:
            nsize  = 100
            flabel = False

        nx.draw_kamada_kawai(G,with_labels=flabel,node_size=nsize/2,width=2,node_color=node_colors) # draw it prettier
        #nx.draw_networkx(G,with_labels=flabel,node_size=nsize,width=2,node_color=node_colors) # draw it pretty
        limits=plt.axis('off')                                      # turn off axes
        plt.show() 

        return

    def SIR_Simulation(self,G,beta,gamma, group):
        # This function simulates a SIR model based on network G with beta and alpha
        #
        # input  : G is a networkx graph
        #        : beta, infection rate per S-I contact
        #        : gamma, I->R recovery rate
        # output : none
        # 
        # WARNING: function is optimistic: assumes inputs are properly formatted

        n  = G.order()
        zt = dict.fromkeys(range(n), 'S') # all nodes S, initially
        #vz is a dictionary containing three labels, 0,1,2 representing no vaccines, 1-dose, 2-dose
        vz = dict.fromkeys(range(n), 0) #all nodes 0 vaccine dosage, initially

        St = [] # S(t), time series of number of S nodes per time step t
        It = [] # I(t), time series of number of I nodes per time step t
        Rt = [] # R(t), time series of number of R nodes per time step t
        Vt = [] # V(t), time series of number of V nodes per time step t

        seed     = int(rnd.randint(0,n-1)) # pick a random node is patient 0
        zt[seed] = 'I'
        t        = 1 #t is in the scale of day
        v_dict = {1 : 0.9, 2 : 0.95} #possibility of S -> V

        print(f'time step {t}')
        self.drawGz(G,zt, t)

        Sc,Ic,Rc, Vc = n-1,1,0,0 # S,I,R,V node counts, initial
        St.append(Sc)
        It.append(Ic)
        Rt.append(Rc)
        Vt.append(Vc)
     
        while any(xi < 2 for xi in vz.values()):
            zu = copy.deepcopy(zt) # nodes states for next time step (synchronous updates)
            
            # do S/I/R -> V transitions
            if t >= 7 and t%7 == 0:
                # number of nodes that will get 1 dose of vaccine:
                vac = int(group * rnd.uniform(1, 2))
                # after a number of time, more vaccines will be provided
                if t >= 20:
                    vac = int(group * rnd.uniform(2, 4))
                    
                # which node will get the vaccine
                v_key = [k for k,h in vz.items() if h < 2]
                
                for i in range(vac):
                    v_ch = rnd.choice(v_key) 
                    vz[v_ch] += 1 #node's vaccine dosage status 0,1,2
                    v_key.remove(v_ch) #one person can only get one vaccine per month
                    v_pr = vz[v_ch] # get the current possibily of vaccination from v_dict
                    if rnd.random() < v_dict[v_pr]:
                        zu[v_ch] = 'V'           # S -> V
                        Vc = Vc+1     # update counts
                        
            # do S -> I transitions
            for e in G.edges():
                i,j = e[0],e[1]           # this edge (i,j)
                if zt[i]=='I' and zt[j]=='S' and zu[j]!='I':
                    if rnd.random() < beta:
                        zu[j] = 'I'       # i infects j for next round
                        Sc,Ic = Sc-1,Ic+1 # update counts

                if zt[i]=='S' and zt[j]=='I' and zu[i]!='I':
                    if rnd.random() < beta:
                        zu[i] = 'I'       # j infects i for next round
                        Sc,Ic = Sc-1,Ic+1 # update counts

            # do I -> R transitions
            for i in G.nodes():
                if zt[i] == 'I' and rnd.random() < gamma:
                    zu[i] = 'R'           # i recovers (R)
                    Ic,Rc = Ic-1,Rc+1     # update counts
                                        
                        
            # update all states synchronously, update clock
            zt = copy.deepcopy(zu)
            t  = t+1
            print(f'time step {t}')
            self.drawGz(G,zt, t)

            St.append(Sc)
            It.append(Ic)
            Rt.append(Rc) 
            Vt.append(Vc)# append these counts to the time series
            
            if t>= 30:
                break
            

        print(f'number of steps in epidemic: {t-1}')
        print(f'final number of S: {Sc}')
        print(f'final number of R: {Rc}')
        print(f'final number of V: {Vc}')
        # plot the S(t),I(t),R(t) time series nicely
        fig = plt.figure()
        ax1 = fig.add_subplot(111) # put multiple 
        plt.plot(range(t), St, 'bo-', alpha=0.5,label='S(t)')  # plot the log-likelihood trajectory
        plt.plot(range(t), It, 'rv-', alpha=0.5,label='I(t)')  # plot the log-likelihood trajectory
        plt.plot(range(t), Rt, 'gs-', alpha=0.5,label='R(t)')  # plot the log-likelihood trajectory
        plt.plot(range(t), Vt, 'cd-', alpha=0.5,label='V(t)')  # plot the log-likelihood trajectory
        plt.ylabel('number of nodes')
        plt.xlabel('time, t')
        plt.legend(loc='upper right');
        plt.show()

        # Changing Vaccine Release Intervals - Connected graph
        # plt.savefig('standard_interval_85_90_connected')
        # plt.savefig('increased_interval_85_90_connected')
        # plt.savefig('decreased_interval_85_90_connected')

        # Vaccine Efficiencies - Connected graph
        # plt.savefig('standard_interval_75_80_connected')
        # plt.savefig('standard_interval_90_95_connected')

        # Changing Vaccine Release Intervals - Relaxed graph
        # plt.savefig('standard_interval_85_90_relaxed')
        # plt.savefig('increased_interval_85_90_relaxed')
        # plt.savefig('decreased_interval_85_90_relaxed')

        # Vaccine Efficiencies - Relaxed graph
        # plt.savefig('standard_interval_75_80_relaxed')
        # plt.savefig('standard_interval_90_95_relaxed')


def connected_caveman_random_partition_graph():
    G = nx.connected_caveman_graph(6,6)
    
    n = 30
    s = 5 #mean cluster size, 2 kids + 2 parents + 2 grandparents
    num_clust = n/s
    v = 4 #variance s/v
    p_in = 0.1 #intra edges probability within clusters, within the family, all members are connected
    p_out = 0.1
    
    for i in range(6):
        g = nx.gaussian_random_partition_graph(n=30, s=5, v=4, p_in=0.1 ,p_out=0.1, directed=False)
        F = nx.disjoint_union(G,g)
        F.add_edge(6*i,len(G))
        G = F
    return G

def relaxed_caveman_random_partition_graph():
    G = nx.relaxed_caveman_graph(6,10,0.1)
    for i in range(6):
        g = nx.gaussian_random_partition_graph(n=30, s=5, v=4, p_in=0.1 ,p_out=0.1, directed=False)
        F = nx.disjoint_union(G,g)
        F.add_edge(6*i,len(G))
        G = F
    return G


class Main():

    # G = connected_caveman_random_partition_graph()
    # SIR = SIR()
    # SIR.SIR_Simulation(G, beta = 0.6, gamma = 0.3, group = 6)
    
    G = relaxed_caveman_random_partition_graph()
    SIR = SIR()
    SIR.SIR_Simulation(G,beta = 0.6,gamma = 0.3, group = 6)

Main()
