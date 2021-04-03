import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pylab as plt
#%matplotlib inline
import random as rnd
rnd.seed()
import copy

class SIR():
    def drawGz(self,G,z):
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
            node_colors.append(colors[int(cid)])
        nsize  = 600
        flabel = True

        if G.order() > 50:
            nsize  = 100
            flabel = False

        nx.draw_kamada_kawai(G,with_labels=flabel,node_size=nsize,width=2,node_color=node_colors) # draw it prettier
        #nx.draw_networkx(G,with_labels=flabel,node_size=nsize,width=2,node_color=node_colors) # draw it pretty
        limits=plt.axis('off')                                      # turn off axes
        plt.show() 

        return

    def SIR_Simulation(self,G,beta,gamma):
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

        St = [] # S(t), time series of number of S nodes per time step t
        It = [] # I(t), time series of number of I nodes per time step t
        Rt = [] # R(t), time series of number of R nodes per time step t

        seed     = int(rnd.randint(0,n-1)) # pick a random node is patient 0
        zt[seed] = 'I'
        t        = 1

        print(f'time step {t}')
        self.drawGz(G,zt)

        Sc,Ic,Rc = n-1,1,0 # S,I,R node counts, initial
        St.append(Sc)
        It.append(Ic)
        Rt.append(Rc)
        while any(xi == 'I' for xi in zt.values()):
            zu = copy.deepcopy(zt) # nodes states for next time step (synchronous updates)

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
            self.drawGz(G,zt)

            St.append(Sc)
            It.append(Ic)
            Rt.append(Rc) # append these counts to the time series

        # report how it went
        print(f'number of steps in epidemic: {t-1}')
        print(f'final number of S: {Sc}')
        print(f'final number of R: {Rc}')

        # plot the S(t),I(t),R(t) time series nicely
        fig = plt.figure()
        ax1 = fig.add_subplot(111) # put multiple 
        plt.plot(range(t), St, 'bo-', alpha=0.5,label='S(t)')  # plot the log-likelihood trajectory
        plt.plot(range(t), It, 'rv-', alpha=0.5,label='I(t)')  # plot the log-likelihood trajectory
        plt.plot(range(t), Rt, 'gs-', alpha=0.5,label='R(t)')  # plot the log-likelihood trajectory
        plt.ylabel('number of nodes')
        plt.xlabel('time, t')
        plt.legend(loc='upper right');
        plt.show()