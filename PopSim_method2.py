from __future__ import division
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import time
from matplotlib.pyplot import pause
from math import *
from math import log


degfile = open('degree.csv', 'w')
matfile = open('matrix.csv', 'w')



# calculate C(n,k)
def comb(n, k):
    return factorial(n) / factorial(k) / factorial(n - k)


# set initial fully connected network as 4 nodes
m0 = 4
init_network_degree = 2 * comb(m0, 2)
pause_time =0.1

#FLAGS - {show degree, plot powerlaw/distribution, plot network respectively}
show_degr = 1
draw_powerlaw = True
plot_network = 0


# node position map for plotting
node_pos = {}

# draw the initial network G with random positions
def plot_initial_network(G):
    for i in G.nodes():
       npos = nx.spring_layout(G)

    #Initial plot with only start nodes
    fig = plt.figure('PoP-Sim Network')
    fig.text(0, 0.97, 'starting nodes: green ',style='italic',fontsize=13)
    fig.text(0, 0.94, 'New node: blue ',style='italic',fontsize=13)
    fig.text(0, 0.91, 'Previously added nodes: red', style='italic',fontsize=13)
    nx.draw_networkx(G,npos,node_color = 'green')
    plt.draw()


#plot network edges

def plot_new_edges(G, new_edges,i):
    plt.clf()
    # Create a color map for nodes to differenctiate between: stating nodes (green), new node (blue) and already added nodes (red)
    color_map = []
    for j in G.nodes():
        if j < m0:
            color_map.append('green')
        elif j == m0 + i :
            color_map.append('blue')
        else: color_map.append('red')
    # Define new node's position and draw the graph
    node_pos= nx.spring_layout(G)
    nx.draw_networkx(G, node_pos, node_color=color_map)
    nx.draw_networkx_edges(G, node_pos,new_edges, width = 2.0 , edge_color = 'b' )
    fig = plt.figure('PoP-Sim Network')
    fig.text(0, 0.97, 'starting nodes: green                 Iteration: '+ str(i+1),style='italic',fontsize=14)
    fig.text(0, 0.94, 'New node: blue ['+str(m0 + i) + ']',style='italic',fontsize=14)
    fig.text(0, 0.91, 'Previously added nodes: red', style='italic',fontsize=14)
    plt.draw()
    pause(pause_time)


def Pop_Sim2( m, gamma):
  #  n=read_nodes()
    sim_matrix = gen_matrix()
    # initialize graph
    G = nx.Graph()
    # Connect all the initial m0 nodes to the graph
    # need to use the nodelist sorted by birth date
    for i in range(1,m0+1):
        G.add_node(i)
        for j in G:
            if (i != j):
                G.add_edge(i, j)
    # could be removed if plotting is not needed
    if plot_network:
        plot_initial_network(G)

    # Adding new node
    N = len(sim_matrix)
    for i in range(1,N - m0+1):
        new_edges = []
        # Compute duration of calculations for consistent timing
        loop_start_time = time.time()
        # select neighbors the new node will connect to according to the  hyperbolic distance
        node_index = choose_neighbour (G,gamma, sim_matrix,m,i+m0)
        G.add_node( i+m0)
        # m must not be less than m0
        for node in node_index:
            G.add_edge(m0+i,node)   #connect new node to m nodes with smallest distance
            new_edges.append((i+m0,node))
            print(new_edges)

        degfile.write(str(new_edges) + "\n")
        if plot_network:
            plot_new_edges(G, new_edges, i)  #plot edges
    plt.close()

    loop_duration = time.time() - loop_start_time
    # Pause for the needed time, taking the calculation time into account
    if pause_time - loop_duration > 0:
        pause(pause_time - loop_duration)

    if draw_powerlaw:
        gen_degree(G)
    # nx.draw(G)
    if show_degr:
        print('Press any key to continue')
        input()

        degr(G)
    else:
        print('Press any key to exit')
        input()



#choose m closest node as returned by hyperbolic function
def choose_neighbour(G,gamma,sim,m,t):
    neigbours = []
    dist_index=new_hyperbolic_dist(gamma,sim,t) #estimate distance and return theri indices
    for i in range(m):     #select m closest(smallest) distances by index
        neigbours.append(dist_index[i])
    nodes=list(G.nodes())
    #use distance indices to select the nodes i.e.
    #if the shortest distance has index 3, connect nodes[3] and so on until m nodes connected
    T=[nodes[i] for i in neigbours]  
    return T


#calculate hyperbolic distance using sim as theta
def new_hyperbolic_dist(gamma,sim,t):
    distance=[]
    beta = 1 / (gamma - 1)
    for s in range (1,t) :  #time starts from 1 not 0(node 0 arrives at t=1)
    # Move all nodes to their new radial coordinates to simulate popularity fading
    # s are existing nodes at time t
        rs = (beta * 2*log(s)) + (1 - beta) *2* (log(t))
        # New node is added to the network and acquires polar coordinates
        rt=2*log(t)
        theta =((sim[t-1][s-1])* (pi)) # sim[t][s]*pi calculated for theta
        d=cosh (rs)* cosh(rt)-sinh(rs) * sinh(rt)*cos(theta)
    # Due to precision problems, numbers d < 1 should be set to 1 to get the right final hyperbolic dis
        if (d<1):
            d=1
            d=acosh(d)
            distance.append(d)
        else:
            d = acosh(d)
            distance.append(d)
    #normalise distance list
    sum_dist = sum(distance)
    distance = [i / sum_dist for i in distance]
    print(distance)
    #sort distance and return distance index
    dist_index = np.argsort(distance)
    return dist_index


#read and process RWR matrix
def gen_matrix():
    rwr_matrix=np.genfromtxt('real_rwr.csv', delimiter=';')
    np.fill_diagonal(rwr_matrix, 0)  # zero diagonal of rwr matrix
    rwr_matrix=normalize_matrix(rwr_matrix)   # normalize rwr matrix
    uniform_matrix=(1/(len(rwr_matrix)-1))+np.zeros((len(rwr_matrix),len(rwr_matrix))) #generate unifortm matrix
    np.fill_diagonal(uniform_matrix, 0)#set uniform matrix diagonal to zero
    sim_matrix = rwr_matrix + uniform_matrix#  add two matrix
    sim_matrix = (sim_matrix / 2) # divide matrix by 2
    sim_matrix = one_minus(sim_matrix)#one_minus trick to make smaller theta (in term of sim) attract
    matrix=normalize_matrix(sim_matrix)#normalize again after 1- trick
    return(matrix)


#read sorted graph nodes ****this has been replaced with time***
def read_nodes():
    nodes = []
    node_file = open('test_serv.txt', 'r')
    for node in node_file:
        temp = node.split(';')
        nodes.append(temp[1])
    return (nodes)


#normalize utility
def normalize_matrix(matrix):
    row_sums = matrix.sum(axis=1)
    new_matrix = matrix / row_sums[:, np.newaxis]
    return new_matrix

#one_minus trick utility
def one_minus(matrix):
    one_matrix=(len(matrix),len(matrix))#generate matrix of ones with size N=sim_matrix
    one_mat=np.ones(one_matrix)
    np.fill_diagonal(one_mat,0) #set one matrix diagonal to zero
    new_mat=np.subtract(one_mat,matrix)# subtract one_matrix from sim_matrix
    return new_mat


# save degree in a text file for further analysis
def gen_degree(G):
    plt.close()
    degfile = open('degreelist.csv', 'w')
    deg_list = []
    for n in G.nodes():
        deg_list.append(G.degree(n))
    degrees = np.array(deg_list)
    for d in deg_list:
        degfile.write(str(d) + '\n')
    return (degrees)



#plot degree distribution of nodes.
def degr(G, scale='log',alpha=.8,low=1,high=1,ec=1):
    plt.close()
    num_nodes=G.number_of_nodes()
    max_degree=0
    for n in G.nodes():
        if G.degree(n)>max_degree:
            max_degree=G.degree(n)
    x=[]
    y_tmp=[]
    for i in range (max_degree+1):
        x.append(i)
        y_tmp.append(0)
        for n in G.nodes():
            if G.degree(n)==i:
                y_tmp[i] += 1
    y=[i/num_nodes for i in y_tmp]
    deg, =plt.plot(x,y,label='Degree Distribution',linewidth=0 , marker='o', markersize=7, color='blue',alpha=alpha)
    # Check for the lin / log parameter and set axes scale
    if scale == 'log':
        plt.xscale('log')
        plt.yscale('log')
        plt.title('Degree distribution (log-log scale)')
        # add theoretical distribution line k^-3
        w = [a for a in range(low,high)]
        z = []
        for i in w:
            x = (i ** -3) * ec  # set line's length and fit intercept
            z.append(x)
        plt.plot(w, z, 'k-', color='#7f7f7f')
    else:
        plt.title('Degree distribution (linear scale)')

    plt.ylabel('Pk')
    plt.xlabel('k')
    plt.title( 'Degree Distribution (Closest Node Method )')
    plt.show()



if __name__ == "__main__":

    Pop_Sim2(4,3) # inputs are m and gamma respectively3
