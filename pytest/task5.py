################################################
### TASK 5
################################################


from itertools import combinations

from graphviz import Digraph
from graphviz import Graph
from scipy.stats.stats import pearsonr

import numpy as np
import pandas as pd


########################################################
##### Reading data and setting static variables
########################################################
df = pd.read_excel("Data_Cortex_Nuclear.xls")
class_list = ['c-CS-m','t-CS-m','c-CS-s','t-CS-s']
col=["#0000ff","#391326","orange","green"]
shp = ['circle','triangle','box','invtriangle']



################################################
### question a)
################################################

## reduction of data
index = []
for i in range(0,len(df)):
    if df['class'].iloc[i] in class_list:
        index.append(i)
        
df = df.ix[index,:]
class_labels = df['class']

## Interpolation to fill missing values
a = df.ix[:,1:78]
a = a.interpolate()
a = a.interpolate(axis=1)

mids = {}                   ## object to contain mouse IDs
i=0
for m in df.MouseID:
    key = m.split("_")[0]
    if key not in mids:
        mids[key] = [i]
    else:
        mids[key].append(i)
    i+=1
# Mean expression values for protein in each mouse have
# been stored in the m_ave_values dictionary

m_ave_values ={}
for key in mids:
    m_ave_values[key] = a.iloc[mids[key],:].mean()

# class for each mouse is in m_class
m_class = {}
for key in mids:
    m_class[key] = df['class'].iloc[mids[key][0]]
    

################################################
### question b)
################################################

## Sample Graph

dot = Digraph(name='pet-shop', node_attr={'shape': 'plaintext'})
dot.node('parrot')
dot.node('dead')
dot.edge('parrot', 'dead')
dot.graph_attr['rankdir'] = 'LR'
dot.edge_attr.update(arrowhead='vee', arrowsize='2')
dot.render(view=True)


################################################
### question c) and d)
################################################

## Shape corresponds to those used in previous tasks
g = Graph(name = 'mouse_cluster')
key_list = mids.keys()                      ## list of mouse IDs
for key in key_list:
    g.attr('node', color=col[class_list.index(m_class[key])], shape = shp[class_list.index(m_class[key])])      ## setting node properties based of mouse class
    g.node(key)                     ## Initialising node
for x in combinations(key_list,2):
    if pearsonr(m_ave_values[x[0]],m_ave_values[x[1]])[0] > 0.98:           ## Check for correlation score
        g.edge(x[0], x[1], penwidth = str(0.03/float(1-pearsonr(m_ave_values[x[0]],m_ave_values[x[1]])[0])))    ## setting up edge, if the condition satisfies, to express the correlation score 
g.view()


################################################
### question e)
################################################


## initialising graph
g = Graph(name = 'Alternate Mouse Cluster')

##Initialising subgraphs based on classed
c0 = Graph('cluster_0')
c0.node_attr.update(color=col[0],shape =shp[0])
for key in key_list:
    if m_class[key] == class_list[0]:
        c0.node(key)

c1 = Graph('cluster_'+class_list[1])
c1.node_attr.update(color=col[1],shape =shp[1])
for key in key_list:
    if m_class[key] == class_list[1]:
        c1.node(key)


c2 = Graph('cluster_'+class_list[2])
c2.node_attr.update(color=col[2],shape =shp[2])
for key in key_list:
    if m_class[key] == class_list[2]:
        c2.node(key)


c3 = Graph('cluster_'+class_list[3])
c3.node_attr.update(color=col[3],shape =shp[3])
for key in key_list:
    if m_class[key] == class_list[3]:
        c3.node(key)
        
## adding subgraphs to main graph
g.subgraph(c0)
g.subgraph(c3)
g.subgraph(c1)
g.subgraph(c2)

## adding edges based on the given condition
for x in combinations(key_list,2):
    if pearsonr(m_ave_values[x[0]],m_ave_values[x[1]])[0] > 0.98:
        g.edge(x[0], x[1], penwidth = str(0.03/float(1-pearsonr(m_ave_values[x[0]],m_ave_values[x[1]])[0])) )
g.view()