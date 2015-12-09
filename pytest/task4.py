'''
Created on Dec 9, 2015

@author: Manish
'''
################################################
### TASK 4
################################################


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA


################################################


def plot_density(data,xs,col):
    density = gaussian_kde(data)
    density.covariance_factor = lambda : .25
    density._compute_covariance()
    plt.plot(xs,density(xs),color = col)
    
    
    
def scatter_plot_matrix(data,class_labels,filename):
    plt.figure(figsize=(2000/96, 2000/96), dpi=96)
    for i in range(0,5):
        plt.subplot(5,5,5*i+i+1)
        plt.title(str(data.columns[i]), fontsize=10)
        min_value = float(np.min(data.ix[:,i]))
        max_value = float(np.max(data.ix[:,i]))
        xs = np.linspace(min_value,max_value,500) 
        k=0
        for cl in class_list:
            ind = np.where(class_labels == cl)[0]
            plot_density(data.ix[ind,i],xs,col[k])
            k+=1
        k=0
        for j in range(0,5):
            if i!=j:
                plt.subplot(5,5,5*i+j+1)
                plt.title(str(data.columns[i]) + " vs " + str(data.columns[j]), fontsize=10)
                plt.xlabel(data.columns[i], fontsize=10)
                plt.ylabel(data.columns[j], fontsize=10)
                for k in range(0,4):
                    ind = np.where(class_labels == class_list[k])[0]
                    plt.scatter(data.ix[ind,i],data.ix[ind,j],s=10,color = col[k],marker = marker_style[k], facecolors = 'none')
    line=[]
    for i in range(0,4):
        line.append(mlines.Line2D([], [], color=col[i], marker=marker_style[i],markersize=10))
    plt.tight_layout()
    plt.figlegend(handles = line, labels = class_list, loc = "upper right")
    plt.savefig(filename)
    plt.show() #refer to the saved image for proper visualization



    
def fValue(protein, classLabel1, classLabel2, data):
    group1 = subgroup(data,classLabel1)
    group2 = subgroup(data,classLabel2)
    meanX = (np.sum(group1[protein])+np.sum(group2[protein]))/(len(group1[protein])+len(group2[protein]))
    mean1 = np.mean(group1[protein])
    mean2 = np.mean(group2[protein])
    var1 = np.var(group1[protein])
    var2 = np.var(group2[protein])
    f = (np.square(mean1-meanX) + np.square(mean2-meanX))/(var1 + var2)
    return(f)



def subgroup(data, cl):
    subDataFrame = pd.DataFrame()
    ind = np.where(class_labels == cl)[0]
    return data.iloc[ind,:]


########################################################
##### Reading data and setting static variables
########################################################

df = pd.read_excel("Data_Cortex_Nuclear.xls")
class_list = ['c-CS-m','t-CS-m','c-CS-s','t-CS-s']
col=["#0000ff","#391326","orange","green"]
marker_style = ['o','^','s','v']



################################################
### question a)
################################################

index = []
for i in range(0,len(df)):
    if df['class'].iloc[i] in class_list:
        index.append(i)

df = df.ix[index,:]
class_labels = df['class']

a = df.ix[:,1:78]
a = a.interpolate()
a = a.interpolate(axis=1)

pca = PCA()
pca.fit(a)
plt.plot(np.cumsum(pca.explained_variance_ratio_))
print("Number of components to cover >=95% variance: ", np.where(np.cumsum(pca.explained_variance_ratio_) >= 0.95)[0][0] + 1)
plt.title("Task4_figureA_Variance_vs_V")
plt.savefig("Task4_figureA_Variance_vs_V.png",type="png")
plt.show()



################################################
### question b)
################################################

pca = PCA(n_components=5)
pca_df = pd.DataFrame(pca.fit_transform(a))
scatter_plot_matrix(pca_df,class_labels,"task4_figureB_PCA_top5.png")


################################################
### question c)
################################################

# From the visual inspection we can see that the outliers are predominant in the 3 pca component
# The outliers are have a 3 component value of more than 2

outlier_in_pca = []
[outlier_in_pca.append(x) for x in np.where(pca_df.ix[:,2] > 2)[0]]
outlier_index = df['MouseID'][[index[x] for x in outlier_in_pca]]
print("\n outlier from mode 3\n")
print(outlier_index)




################################################
### question d)
################################################

a_new = a.drop(a.index[outlier_in_pca])
pca = PCA(n_components=5)
pca_df = pd.DataFrame(pca.fit_transform(a_new))
class_labels_new = pd.Series([cl for i, cl in enumerate(class_labels) if i not in outlier_in_pca])
scatter_plot_matrix(pca_df,class_labels_new,"task4_figureD_PCA_after_removing_outlier.png")

# Individual mice cluster analysis

# from the scatter plot matrix, first possible cluster has component1<-2.5 and component2<0.5
outlier_in_pca = []
l1 = np.where(pca_df.ix[:,1] < -2.5)[0]
l2 = np.where(pca_df.ix[:,2] < 0.5)[0]
outlier_in_pca = list(set(l1).intersection(l2))
outlier_index = df['MouseID'][[index[x] for x in outlier_in_pca]]
print("\n\n First possible cluster mouseIds\n")
print(outlier_index)
# Shows that the cluster is composed of 1 mice samples only; MouseID ==> 361

# another cluster has component2>0.8 and component4<-0.5
outlier_in_pca = []
l1 = np.where(pca_df.ix[:,2] > 0.8)[0]
l2 = np.where(pca_df.ix[:,4] < -0.5)[0]
outlier_in_pca = list(set(l1).intersection(l2))
outlier_index = df['MouseID'][[index[x] for x in outlier_in_pca]]
print("\n\n Second possible cluster mouseIds\n")
print(outlier_index)
# Shows that the cluster is composed of 1 mice samples only; MouseID ==> 322


################################################
### question e)
################################################


c_cs_ind = np.where(class_labels == "c-CS-s")[0]

# To calculate deviation from baseline I am using ratio as the measure, 
# alternatively absolute change can also be used
# but the statistical power of the two approaches would vary and 
# calculation of that is beyond the scope of this study

a_norm = a.div(a.iloc[c_cs_ind,:].mean())
pca = PCA()
pca.fit(a_norm)
print("Number of components to cover >=95% variance: ", np.where(np.cumsum(pca.explained_variance_ratio_) >= 0.95)[0][0] + 1)




################################################
### question f)
################################################

## computing Fvalue
proteinFvalue = []
for i in range(0,len(a_norm.columns)):
    proteinFvalue.append(fValue(a_norm.columns[i],"t-CS-s","c-CS-s",a_norm))
## Multiply by FValue
a_norm_F = a_norm.mul(proteinFvalue)

# Scatter plot without re-weighing
plt.subplot(1,2,1)
plt.title("without re-weighing")
pca = PCA(n_components=2)
pca_df = pd.DataFrame(pca.fit_transform(a_norm))

for k in range(0,4):
    ind = np.where(class_labels == class_list[k])[0]
    plt.scatter(pca_df.ix[ind,0],pca_df.ix[ind,1],s=10,color = col[k],marker = marker_style[k], facecolors = 'none')

# Scatter plot with re-weighing
plt.subplot(1,2,2)
plt.title("with re-weighing")
pca = PCA(n_components=2)
pca_df = pd.DataFrame(pca.fit_transform(a_norm_F))
for k in range(0,4):
    ind = np.where(class_labels == class_list[k])[0]
    plt.scatter(pca_df.ix[ind,0],pca_df.ix[ind,1],s=10,color = col[k],marker = marker_style[k], facecolors = 'none')

line=[]
for i in range(0,4):
    line.append(mlines.Line2D([], [], color=col[i], marker=marker_style[i],markersize=10))
plt.tight_layout()
plt.figlegend(handles = line, labels = class_list, loc = "upper right")
plt.savefig("task4_figureF_PCA_after_removing_outlier.png")
plt.show()