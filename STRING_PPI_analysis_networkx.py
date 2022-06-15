import pandas as pd
import networkx as nx
import numpy as np
import collections
import matplotlib.pyplot as plt

cutoff = 400

def plot_degree(g,filename):
    degree_sequence = sorted([d for n, d in g.degree()], reverse=True)
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())

    fig, ax = plt.subplots()
    plt.bar(deg, cnt, width=0.80, color='b')

    plt.title("Degree Histogram")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    ax.set_xticks([d + 0.4 for d in deg])
    ax.set_xticklabels(deg)

    plt.axes([0.4, 0.4, 0.5, 0.5])
    pos = nx.spring_layout(g)
    plt.axis('off')
    nx.draw_networkx_nodes(g, pos, node_size=20)
    nx.draw_networkx_edges(g, pos, alpha=0.4)

    plt.show()
    filename = filename+".png"
    fig.savefig(filename)

column_names = ["Organism","Nodes","Edges", "Avg CC", "Avg BC", "Max Degree", "Avg Degree", 
                "Avg SP", "Connected","Diameter"]
networks = pd.DataFrame(columns = column_names)

k=0

organisms = {
    "Homo sapiens": 9606,
    "Mus musculus": 10090,
    "Arabidopsis thaliana": 3702,
    "Saccharomyces cerevisiae": 4932,
    "Escherichia coli": 511145,
    "Caenorhabditis elegans": 6239,
    "Rattus norvegicus": 10116,
    "Drosophila melanogaster": 7227,
    "Bacillus subtilis": 224308,
    "Pseudomonas aeruginosa PAO1": 287
}

for org in organisms:

    id_string = organisms[org]
    name = org

    ppi_link = "https://stringdb-static.org/download/protein.links.v11.5/"+str(id_string)+".protein.links.v11.5.txt.gz"
    df = pd.read_csv(ppi_link, compression='gzip', header=0, sep='\s', engine="python")

    final_df = df.loc[df["combined_score"] >= cutoff]
    #print(name)

    g = nx.from_pandas_edgelist(final_df,'protein1','protein2',edge_attr=True)

    #metrics
    nodes = g.number_of_nodes()
    edges = g.number_of_edges()
    cc = nx.closeness_centrality(g)
    avg_cc = np.mean(list(cc.values()))
    bc = nx.betweenness_centrality(g)
    avg_bc = np.mean(list(bc.values()))
    degree_sequence = sorted((d for n, d in g.degree()), reverse=True)
    dmax = max(degree_sequence)
    davg = np.mean(degree_sequence)
    connect = nx.is_connected(g)
    diam = nx.diameter(g) if connect else "None"
    avg_sp = nx.average_shortest_path_length(g) if connect else "None"

    networks.loc[k] = [name,nodes,edges,avg_cc,avg_bc,dmax,davg,avg_sp,connect,diam]
    k=k+1

    plot_degree(g,name)

networks.to_csv("all_metrics.csv")
