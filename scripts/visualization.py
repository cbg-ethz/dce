# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import os

import numpy as np
import pandas as pd
import networkx as nx

import seaborn as sns
import matplotlib.pyplot as plt
# -

sns.set_context('talk')

# # Load data

# +
df_normal = pd.read_csv(snakemake.input['normal_file'])
df_tumor = pd.read_csv(snakemake.input['tumor_file'])

df_normal['type'] = 'normal'
df_tumor['type'] = 'tumor'

df_graph = pd.concat([df_normal, df_tumor])
df_graph.head()

# +
graph_normal = nx.from_pandas_edgelist(
    df_graph[df_graph['type'] == 'normal'],
    source='source', target='sink', edge_attr=True,
    create_using=nx.MultiDiGraph)

graph_tumor = nx.from_pandas_edgelist(
    df_graph[df_graph['type'] == 'tumor'],
    source='source', target='sink', edge_attr=True,
    create_using=nx.MultiDiGraph)
# -

# # Basic statistics

# +
plt.figure(figsize=(8,6))

sns.distplot(df_graph.loc[df_graph['type'] == 'normal', 'causal.effect'], kde=False, label='normal')
sns.distplot(df_graph.loc[df_graph['type'] == 'tumor', 'causal.effect'], kde=False, label='tumor')

plt.legend(loc='best')


# -

# # Plot graph

@np.vectorize
def symlog(x):
    if x > 0:
        return np.log(x)
    elif x < 0:
        return -np.log(abs(x))
    else:
        return 0


def plot(graph, ax):
    # layout
    pos = nx.drawing.nx_agraph.graphviz_layout(graph, prog='neato', args='-Goverlap=false')

    # visualization
    nx.draw_networkx_nodes(graph, pos, ax=ax)
    nx.draw_networkx_labels(
        graph, pos, ax=ax,
        labels={n: n for n in graph.nodes()})

    nx.draw_networkx_edges(
        graph, pos, ax=ax,
        edge_color=symlog([e[-1]['causal.effect'] for e in graph.edges(data=True)]),
        edge_cmap=plt.get_cmap('bwr'))
    nx.draw_networkx_edge_labels(
        graph, pos, ax=ax,
        edge_labels={e[:2]: '{:.2}'.format(e[-1]['causal.effect'])
                     for e in graph.edges(data=True)})


# +
plt.figure(figsize=(24,12))
plt.gcf().set_facecolor('w')

ax = plt.subplot(121)
ax.set_axis_off()
plot(graph_normal, ax)
plt.title('Normal')

ax = plt.subplot(122)
ax.set_axis_off()
plot(graph_tumor, ax)
plt.title('Tumor')

plt.tight_layout()
plt.savefig(snakemake.output['img_file'])
