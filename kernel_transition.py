import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
label_size = 25
plt.rcParams['xtick.labelsize'] = label_size 
plt.rcParams['ytick.labelsize'] = label_size 
mpl.rc('text', usetex = True)

def kernel(A, nodes, Nsq):
    '''
    compute spatial kernel for a mobility matrix A
    '''
    K = np.zeros(Nsq/2)
    norm = np.zeros(Nsq/2)
    x = [i + 0.5 for i in range(Nsq/2)]
    N = len(nodes)
    for i in range(N):
        for j in range(i, N):
            d = dist(nodes[i], nodes[j])
            d = int(d)
            if d <= Nsq/2 - 1:
                norm[d] += 1.0
                K[d] += A[i][j]
    for d in range(Nsq/2):
        K[d] /= norm[d]
    return x, K

def dist(node1, node2):
    '''
    distance in the 2d lattice
    '''
    d = np.sqrt((node1[0] - node2[0])**2 + (node1[1] - node2[1])**2)
    return d

#generate a square network
Nsq = 50
G = nx.grid_2d_graph(Nsq, Nsq, periodic = False)
nodes = G.nodes()
N = len(nodes)

#get the adjacency matrix
T = np.asarray(nx.to_numpy_matrix(G))
#add staying at the same cell
for i in range(N):
    T[i,i] = 4.0
#normalize transition matrix
for i in range(N):
    kk = np.sum(A[i,:])
    T[i,:] /= kk

x, K = kernel(A, nodes, Nsq)

plt.plot(x, K, label = '$\\tau_i = 1$')

for i in range(8):
    print(i)
    T = T.dot(T)
    x, K = kernel(A, nodes, Nsq)
    plt.plot(x, K, label = '$\\tau_i = %i$' % (i+2))

plt.legend(fontsize = 10)
plt.xlabel('$d$', fontsize = 30)
plt.ylabel('$K(d)$', fontsize = 30)

#no log scales, looks awful
#plt.savefig('kernel_transition_incubation_nolog.png',bbox_inches='tight', dpi=300)

plt.yscale('log')

#log in y axis
plt.savefig('kernel_transition_incubation_logy.png',bbox_inches='tight', dpi=300)


#loglog plot
#plt.xscale('log')
#plt.savefig('kernel_transition_incubation_logxy.png',bbox_inches='tight', dpi=300)

plt.show()
