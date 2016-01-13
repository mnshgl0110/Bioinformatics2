################################################
### TASK 3
################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats import gaussian_kde
from scipy.spatial.distance import euclidean as de
from itertools import combinations as ic




################################################

## Function to plot density plots
def plot_density(data,xs,clr):
    density = gaussian_kde(data)
    density.covariance_factor = lambda : .25
    density._compute_covariance()
    plt.plot(xs,density(xs),color=clr)
    
    
## Function to calculate distance consistency
def distance_consistency(pro1, pro2):
    centers={}
    for cl in class_list:
        index = np.where(df['class'] == cl)[0]
        centers[cl] = [np.mean(pro1[index]),np.mean(pro2[index])]
    count = 0
    for cl in class_list:
        index = np.where(df['class'] == cl)[0]
        for i in index:
            if de(centers[cl], [pro1[i],pro2[i]]) < np.max([de(centers[c],[pro1[i],pro2[i]]) for c in class_list]):
                count+=1
    return 1 - (count/len(pro1))


################################################
##### Reading data and setting static variables
################################################
df = pd.read_csv("reduced_df.csv")
class_list = ['c-CS-m','t-CS-m','c-CS-s','t-CS-s']
col=["#0000ff","#391326","orange","green"]
marker_style = ['o','^','s','v']



################################################
### question a)
################################################

plt.figure(figsize=(2000/96, 2000/96), dpi=96)      ##opening an image object and setting image size

for i in range(0,5):                                ## for top 5 proteins of the original data set
    plt.subplot(5,5,5*i+i+1)                        ## the density function to be plotted on the diagonals
    plt.title(df.columns[i+1], fontsize=10)
    min_value = float(np.min(df.ix[:,i+1]))
    max_value = float(np.max(df.ix[:,i+1]))
    xs = np.linspace(min_value,max_value,500)       ## setting linespace for possible values of the protein
    k=0
    for cl in class_list:
        index = np.where(df['class'] == cl)[0]      ##selecting samples consisting of individuals classes
        plot_density(df.ix[index,i+1],xs,col[k])    
        k+=1
    k=0
################################################
### question b)
################################################
    
    for j in range(0,5):
        if i!=j:
            plt.subplot(5,5,5*i+j+1)                ## non-diagonal subplot for scatter plots
            plt.title(df.columns[i+1] + " vs " + df.columns[j+1], fontsize=10)
            plt.xlabel(df.columns[i+1], fontsize=10)
            plt.ylabel(df.columns[j+1], fontsize=10)
            for k in range(0,4):
                index = np.where(df['class'] == class_list[k])[0]
                plt.scatter(df.ix[index,1+i],df.ix[index,1+j],s=10,color = col[k],marker = marker_style[k], facecolors = 'none')


## setting the legend 
line=[]
for i in range(0,4):
    line.append(mlines.Line2D([], [], color=col[i], marker=marker_style[i],markersize=10))
plt.tight_layout()
plt.figlegend(handles = line, labels = class_list, loc = "upper right")

plt.savefig("task3_figure1.png")
plt.show() #refer to the saved image for proper visualization



################################################
### question c)
################################################


#Working on the assumption that sample class would be the dominating clustering factor,
#we can calculate distance consistency from the data frame itself.
#Also, we need not to calculate the distance consistency values twice for the (A,B)/(B,A) protein
#pairs because the results calculation would be order independent.
          
proteinPairList = [list(x) for x in ic(range(1,6),2)]       ## making all possible pairs of index for 5 proteins
max_d = 0
for x in proteinPairList:
    d = distance_consistency(df.ix[:,x[0]],df.ix[:,x[1]])   ## calculation of distance consistency score
    print("distance consistency for",df.columns[x[0]],"and",df.columns[x[1]],"=",d)
    if d > max_d:
        max_d = d
        proteinPair = x
print("maximum distance consistency is for",df.columns[proteinPair[0]],"and",df.columns[proteinPair[1]],"=",max_d)


################################################
### question d)
################################################
k=0
plt.title('t-CS-m')
plt.xlabel(df.columns[1], fontsize=10)
plt.ylabel(df.columns[5], fontsize=10)
index = np.where(df['class'] == 't-CS-m')[0]
plt.scatter(df.ix[index,1],df.ix[index,5],s=10,color = col[k],marker = marker_style[k], facecolors = 'none')
plt.show()

plt.title('c-CS-s')
plt.xlabel(df.columns[1], fontsize=10)
plt.ylabel(df.columns[5], fontsize=10)
index = np.where(df['class'] == 'c-CS-s')[0]
plt.scatter(df.ix[index,1],df.ix[index,5],s=10,color = col[k],marker = marker_style[k], facecolors = 'none')
plt.show()