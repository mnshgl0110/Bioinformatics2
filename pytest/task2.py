################################################
### TASK 2
################################################

import pandas as pd
import numpy as np



################################################

###### Function to extract subgroups
def subgroup(df, classLabel):
    subDataFrame = pd.DataFrame()
    for i in range(0,len(df)):
        if(df["class"].iloc[i] == classLabel):
            subDataFrame[i] = df.iloc[i]
    return subDataFrame.T


###### Function to calculate F-score
def fValue(protein, classLabel1, classLabel2, df):
    group1 = subgroup(df,classLabel1)
    group2 = subgroup(df,classLabel2)
    meanX = (np.sum(group1[protein])+np.sum(group2[protein]))/(len(group1[protein])+len(group2[protein]))
    mean1 = np.mean(group1[protein])
    mean2 = np.mean(group2[protein])
    var1 = np.var(group1[protein])
    var2 = np.var(group2[protein])
    f = (np.square(mean1-meanX) + np.square(mean2-meanX))/(var1 + var2)
    return(f)


################################################
### question a)
################################################

# Read file
df = pd.read_excel("Data_Cortex_Nuclear.xls")
print("Number of Instances:",len(df))  #Number of instances
print("Number of Columns:",len(df.columns)) #Number of columns
print("Colnames:", list(df.columns)) #Colnames



################################################
### question b)
################################################

## Interpolation strategy requires that First matrix element (1,1) be present and not NaN
#Interpolating to fill NaN values in each column
df = df.interpolate()

#Interpolation to fill NaN values which are in the beginnning of column (case 'BCL_2')
df.ix[:,1:78] = df.ix[:,1:78].interpolate(axis=1)


################################################
### question c)
################################################
            
# Extracting subgroups
group_tCSsal = subgroup(df,"t-CS-s")
group_cCSsal = subgroup(df,"c-CS-s")
group_tCSmem = subgroup(df,"t-CS-m")
group_cCSmem = subgroup(df,"c-CS-m")

# Printing number of instances
print("Number of instances of class t-CS-s", len(group_tCSsal))
print("Number of instances of class c-CS-s",len(group_cCSsal))
print("Number of instances of class t-CS-m",len(group_tCSmem))
print("Number of instances of class c-CS-m",len(group_cCSmem))



################################################
### question d)
################################################
proteinFvalue = {}
for i in range(1,len(df.columns)-4):
    proteinFvalue[df.columns[i]] = fValue(df.columns[i],"t-CS-s","c-CS-s",df)
# Top 5 proteins are 
top5 = sorted(proteinFvalue.items(),key = lambda x: x[1], reverse = True)[:5]



################################################
### question e)
################################################

proteinIDs = []
for values in top5:
    proteinIDs = proteinIDs + [values[0]]
reduceDF = pd.concat([group_cCSmem,group_cCSsal,group_tCSmem,group_tCSsal])
reduceDF = reduceDF.loc[:,["MouseID"]+proteinIDs+["Genotype","Treatment","Behavior","class"]]
print("top 5 proteins:", list(reduceDF.columns[1:6]))
reduceDF.to_csv("reduced_df.csv",index=False)