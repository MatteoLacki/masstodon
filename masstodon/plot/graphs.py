try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass
import networkx as nx


def plot_numbered_graph(G, show=True):
    colors= ['red' if N > 0 else 'blue' for N in G]
    pos   = nx.spring_layout(G)
    nodes = nx.draw_networkx_nodes(G,
                                   pos = pos,
                                   node_color = colors)
    edges = nx.draw_networkx_edges(G,
                                   pos = pos)
    if show:
        plt.show()
