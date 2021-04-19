def connected_caveman_random_partition_graph():
    G = nx.connected_caveman_graph(6,6)
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
