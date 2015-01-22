"""Usage: metagenome_PCA.py (--infile FILE --dim=<dimensions>) [--kmeans=<clust_no>]
       metagenome_PCA.py -h | --help
    Arguments:
    FILE  the infile.
    dimensions  the dimensions for PCA.

Options:
    -h --help  This is a help message.
    -i=FILE --infile=FILE  The input file from IMG/M.
    --kmeans=<clust_no>, -k=<clust_no>  Number of clusters to generate with kmeans [default: 1].
    -d=<dimensions>, --dim=<dimensions>  The dimensions to reduce to by PCA [default: 3]. 

"""
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from docopt import docopt
import numpy as np
import pandas as pd
import plotly.plotly as py
from plotly.graph_objs import *
import random
py.sign_in("wesfield", "8u7gqn6uul")

def k_means(data,clusters=8):
    model = KMeans(n_clusters=clusters, max_iter=300, init='k-means++')
    labels=model.fit(data).predict(data)
    to_return= [[] for x in range(clusters)]
    for  i,entry in enumerate(data):
        to_return[labels[i]].append(entry)
    return to_return
    
def make_trace(scatter_data,color='rgba(217, 217, 217, 0.14)'):
    X=[]
    Y=[]
    Z=[]
    for entry in scatter_data:
        X.append(entry[0])
        Y.append(entry[1])
        Z.append(entry[2])
        
    trace=Scatter3d(
        x=X,y=Y,z=Z, mode='markers',
        marker=Marker(
            size=12,
            line=Line(
                color=color,
                width=0.5
                ),
                opacity=0.8
            )
        )
    return trace
        

if __name__=='__main__':
    
    #Clean data and discard all non-bacterial info
    arguments = docopt(__doc__)
    df= pd.DataFrame.from_csv(arguments['--infile'][0], sep='\t', header=0,index_col=None)
    df=df[df['Domain']=='Bacteria']
    df=df.drop('Genome Count',1)
    df=df.filter(regex='^(?!.*Count).*$').transpose()
    df.columns=df.iloc[1]
    df.reindex(df.index.drop(['Domain','Phylum']))
    newdf=pd.DataFrame(data=df[2:-1],columns= df.columns)
    newdf.to_csv(path_or_buf='testing.csv', sep='\t')
    
    #Perform PCA
    metaPCA=PCA(n_components=int(arguments['--dim'][0]))
    metaPCA.fit(newdf)
    scatter_data= metaPCA.transform(newdf)
    

    clusters= k_means(scatter_data,clusters=int(arguments['--kmeans'][0]))
    #make plotly scatter
    traces=[]
    for cluster in clusters:
        r = lambda: random.randint(0,255)
        rcolor = 'rgba('+str(r())+','+str(r())+','+str(r())+',.89)'
        traces.append(make_trace(cluster,color=rcolor))
    
    data = Data(traces)
    
    layout = Layout(
        margin=Margin(
            l=0,
            r=0,
            b=0,
            t=0
            )
        )
    fig=Figure(data=data,layout=layout)
    plot_url=py.plot(fig,filename='Metagenome_3dScatter')
    
    