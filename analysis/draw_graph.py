'''
Created on 10. apr. 2014

@author: H, inspired by udacity example
'''

import networkx as nx
import matplotlib.pyplot as plt

COLORS = ['blue', 'red', 'green', 'black', 'orange',
          'magenta', 'yellow', 'cyan', 'pink', 'brown']


def draw_graph(nodesets, edges, nodepos):
    
    print "\n DRAW GRAPH"
    G=nx.Graph()
    
#     G.add_node(0, pos=(1.5, 1.5), color='m')
#     G.add_node(1, pos=(1.5, 2.5))
#     G.add_node(2, pos=(0.5, 2.5))
    
    # add all nodes to graph
    for node in nodepos:
#         print node
        p = (float(node[1]), float(node[2]))
#         print p
        G.add_node(int(node[0]),pos=p )
#             G.node[node[0]]['pos'] = p
#     G.add_edge(0,1)
    
    # add all edges to the graph
    for edge in edges:
        G.add_edge(int(edge[0]),int(edge[1]) )


    # get position map
    pos=nx.get_node_attributes(G,'pos')
#     col = nx.get_node_attributes(G,'color')
    
#     print "positions?", pos
#     nx.draw_networkx_nodes(G, pos=pos, nodelist, node_size, node_color, node_shape, alpha, cmap, vmin, vmax, ax, linewidths, label)
        
    nx.draw(G, pos)
    
    # draw nodes with different colors
    l = len(nodesets)
    for nodeset, c in zip(nodesets, COLORS[:l]):
        l = [int(el) for el in nodeset]
        print "adding nodes color:", c, l, sorted(nodeset)
        nx.draw_networkx_nodes(G, pos, node_color=c, nodelist=l)
        
    nx.draw_networkx_edges(G, pos)

    plt.show()
    
def draw_graph4(nodepos, edges, nodecolors):
    """ draw graph using:
    list of positions for all nodes, [nodenr, x, y]
    list of edge pairs in the graph, [a, b]
    list of colors corresponding to nodepos
    """
    # using the map node--> color
    print "\n DRAW GRAPH"
    G=nx.Graph()
    
    # add all nodes to graph
    for node in nodepos:
        p = (float(node[1]), float(node[2]))
        G.add_node(int(node[0]),pos=p )

    # add all edges to the graph
    for edge in edges:
        G.add_edge(int(edge[0]),int(edge[1]) )

    # get position map
    pos=nx.get_node_attributes(G,'pos')
    nx.draw(G, pos) # ugly way to get labels:p
    nx.draw_networkx_nodes(G, pos, node_color=nodecolors)
    nx.draw_networkx_edges(G, pos)
    
    plt.show()

def draw_graph3(graph):

    # extract nodes from graph
    nodes = set([n1 for n1, n2 in graph] + [n2 for n1, n2 in graph])
    
    for n in nodes:
        print n

    # create networkx graph
    G=nx.Graph()

    # add nodes
    for node in nodes:
        G.add_node(node)

    # add edges
    for edge in graph:
        G.add_edge(edge[0], edge[1])

    # draw graph
#     pos = nx.shell_layout(G)
    pos = nx.graphviz_layout(G)
    nx.draw(G, pos)

    # show graph
    plt.show()

# draw example
# graph = [(20, 21),(21, 22),(22, 23), (23, 24),(24, 25), (25, 20)]
# draw_graph(graph)



def draw_graph2(graph, labels=None, graph_layout='shell',
               node_size=1600, node_color='blue', node_alpha=0.3,
               node_text_size=12,
               edge_color='blue', edge_alpha=0.3, edge_tickness=1,
               edge_text_pos=0.3,
               text_font='sans-serif'):

    # create networkx graph
    G=nx.Graph()

    # add edges
    for edge in graph:
        G.add_edge(edge[0], edge[1])

    # these are different layouts for the network you may try
    # shell seems to work best
    if graph_layout == 'spring':
        graph_pos=nx.spring_layout(G)
    elif graph_layout == 'spectral':
        graph_pos=nx.spectral_layout(G)
    elif graph_layout == 'random':
        graph_pos=nx.random_layout(G)
    else:
        graph_pos=nx.shell_layout(G)

    # draw graph
    nx.draw_networkx_nodes(G,graph_pos,node_size=node_size, 
                           alpha=node_alpha, node_color=node_color)
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
                           alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,
                            font_family=text_font)

    if labels is None:
        labels = range(len(graph))

    edge_labels = dict(zip(graph, labels))
    nx.draw_networkx_edge_labels(G, graph_pos, edge_labels=edge_labels, 
                                 label_pos=edge_text_pos)

    # show graph
    plt.show()
    
    


# graph = [(0, 1), (1, 5), (1, 7), (4, 5), (4, 8), (1, 6), (3, 7), (5, 9),
#          (2, 4), (0, 4), (2, 5), (3, 6), (8, 9)]
# 
# you may name your edge labels
# labels = map(chr, range(65, 65+len(graph)))
#draw_graph(graph, labels)

# if edge labels is not specified, numeric labels (0, 1, 2...) will be used
# draw_graph2(graph)