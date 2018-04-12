from __future__ import division
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import time
from matplotlib.pyplot import pause
from math import *



degfile=open('degree.csv','w')
matfile=open('matrix.npy','w')


# calculate C(n,k)

def comb(n,k):
    return factorial(n) / factorial(k) / factorial(n - k)

        
# set initial fully connected network as 4 nodes
m0 = 4
init_network_degree = 2*comb(m0,2)
pause_time = 0.1

# Flags - {show degree, plot powerlaw, plot network respectively}
show_degr = 1
draw_distr = 1
plot_network=0

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
    nx.draw_networkx(G,npos, node_color = 'green')
    plt.draw()

def plot_new_edges(G, new_edges, i):
    plt.clf()
    # Create a color map for nodes to differenctiate between: stating nodes (green), new node (blue) and already added nodes (red)
    color_map = []
    for j in G.nodes():
        if int(j) < m0:
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

def Pop_Sim(m,gamma):
    sim_matrix=gen_matrix()
    #n = read_nodes() # read API nodes== t
    # initialize graph
    G = nx.Graph()
    # Connect all the initial m0 nodes to the graph
    # need to use the nodelist sorted by birth date

    for i in range(m0):
        G.add_node(i)
        for j in G:
            if (i != j):
                G.add_edge(i, j)

     # could be removed if plotting is not needed
    if plot_network:
        plot_initial_network(G)

    # Adding new nodes
    N=len(sim_matrix)
    for i in range(1,((N-m0)+1)):
        # Compute duration of calculations for consistent timing
        loop_start_time = time.time()
       # select neighbors the new node will connect to according to the  hyperbolic distance
        neighbors = choose_neighbour(G,gamma,m,sim_matrix,(i+m0))
        # A Check to make sure the correct number of neighbors are chosen
        if (len(neighbors) != m):
            print ("Error, number of neighbors is not as expected")
            return
        # Add the new node to the graph
        G.add_node(m0 + i)
        # Save new edges in a list for drawing purposed
        new_edges = []

        for nb in neighbors:
            G.add_edge(m0 + i, nb)
            new_edges.append((m0 + i,nb))
            print(new_edges)
        degfile.write(str(new_edges)+"\n")

        if plot_network:
            plot_new_edges(G, new_edges, i)

    plt.close()

    loop_duration = time.time() - loop_start_time
    # Pause for the needed time, taking the calculation time into account
    if pause_time - loop_duration > 0:
        pause(pause_time - loop_duration)

    if draw_distr:
        gen_degree(G)

    if show_degr:
        print ('Press any key to continue')
        input()
        degr(G)
    else:
        print ('Press any key to exit')
        input()



#Randomly select m nodes without replacement based on hyperpolic distance
def choose_neighbour(G, gamma, m, sim, t):
    dist = new_hyperbolic_dist(gamma, sim, t)
   # print(dist)
    connect_to = np.random.choice(list(G.nodes()), m, replace=False, p=dist)
    return (connect_to)


# calculate hyperbolic distance using sim as theta
def new_hyperbolic_dist(gamma, sim, t):
    distance = []
    beta = 1 / (gamma - 1)
    for s in range(1, t):
        # Move all nodes to their new radial coordinates to simulate popularity fading
        # s are existing nodes at time t
        rs = (beta * 2 * log(s)) + (1 - beta) * 2 * (log(t))
        # New node is added to the network and acquires polar coordinates
        rt = 2 * log(t)
        theta = ((sim[t - 1][s - 1]) * (pi))  # sim[t][s]*pi calculated for theta
        d = cosh(rs) * cosh(rt) - sinh(rs) * sinh(rt) * cos(theta)
        # Due to precision problems, numbers d < 1 should be set to 1 to get the right final hyperbolic dis
        if (d < 1):
            d = 1
            d = acosh(d)
            distance.append(d)
        else:
            d = acosh(d)
            distance.append(d)
    # normalise distance list
    sum_dist = sum(distance)
    distance = [i / sum_dist for i in distance]
    return distance

#read sorted graph nodes*****ds can be replaced with time****
def read_nodes():
    nodes = []
    node_file = open('test_serv.txt', 'r')
    for node in node_file:
        temp = node.split(';')
        nodes.append(temp[1])
    return (nodes)


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

#normalize utility
def normalize_matrix(matrix):
    row_sums = matrix.sum(axis=1) 
    new_matrix = matrix / row_sums[:, np.newaxis]# row normalization
    return new_matrix

#one_minus trick utility
def one_minus(matrix):
    one_matrix=(len(matrix),len(matrix))#generate one_martrix of size input matrix
    mat_one=np.ones(one_matrix)  #generate matrix of ones
    np.fill_diagonal(mat_one,0) #set matrix of one to zero
    new_mat=np.subtract(mat_one,matrix)  #subtract matrix of 1 from the rwr_matrix
    return new_mat


#save degree in a text file for further analysis 
def gen_degree(G):
    plt.close() 
    degfile = open('degreelist.csv','w')
    deg_list=[]
    for n in G.nodes():
        deg_list.append(G.degree(n))
       # deg_list.append((n, G.degree(n)))
    degrees=np.array(deg_list)
    for d in deg_list:
        degfile.write(str(d)+'\n')
    return(degrees)



#plot degree distibution
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
    deg, =plt.plot(x,y,label='Degree Distribution',linewidth=0 , marker='o', markersize=5, color='r',alpha=alpha)
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

    plt.ylabel('P(k)')
    plt.xlabel('k')
    plt.title( 'Degree Distribution (Choice Method)')
    plt.show()


if __name__ == "__main__":
    #gen_matrix()
    Pop_Sim(1,2.1)


    
