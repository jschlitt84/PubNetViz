import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sb
from datetime import datetime
from Bio import Entrez
from itertools import chain

from networkx.algorithms.community import greedy_modularity_communities



def search(query,results='99999'):
    """Runs Entrez search and returns article handes"""
    Entrez.email = email
    handle = Entrez.esearch(db='pubmed', 
                            sort='date', 
                            retmax=str(results),
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    """Pulls article contents by article ID list"""
    ids = unicode(','.join(id_list))
    Entrez.email = email
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id = ids)
    results = Entrez.read(handle)
    return results


def pullResults(query,results):
    """Returns papers json from pubmed query"""
    articleHits = search(query,results)
    articleIDs = articleHits['IdList']
    print query, len(articleIDs)
    if True:
        papers = fetch_details(articleIDs)
    else:
        papers = []
    return papers
    
def queryToDf(query,results,show=False,dumpRaw=False):
    """Creates df of select pubmed data fields from query"""
    papers = pullResults(query,results)
    data = pd.DataFrame(columns=['year','month','day','title',
                                 'abstract','qualifiers','grantAgency',
                                 'grantID'])
    years = []
    titles = []
    abstracts = []
    noLoads = []
    toIterate = papers['PubmedArticle'] + papers['PubmedBookArticle']
    if dumpRaw:
        return toIterate,0
    for paper in toIterate:
        title = abstract = year = 'null'
        try:
            title = paper['MedlineCitation']['Article']['ArticleTitle']
        except:
            noLoads.append(['bad title',paper])
        try:
            year = int(paper['MedlineCitation']['DateCreated']['Year'])
            month = int(paper['MedlineCitation']['DateCreated']['Month'])
            day = int(paper['MedlineCitation']['DateCreated']['Day'])
            citeFormat = 'Medline'
        except:
            try:
                year = int(paper['PubmedData']['History'][0]['Year'])
                month = int(paper['PubmedData']['History'][0]['Month'])
                day = int(paper['PubmedData']['History'][0]['Day'])
                citeFormat = 'PubmedData'
            except:
                noLoads.append(['no year',paper])
                citeFormat = 'NotFound'
        try:
            abstract = unicode(paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0])
        except:
            noLoads.append(['no abstract',paper])
            
        try:
            meshHeadings = paper['MedlineCitation']['MeshHeadingList']
            qualifiers = set(chain(*[[str(j) for j in i['QualifierName']] for i in meshHeadings]))
        except:
            qualifiers = {}
        
        try:
            grantList = paper['MedlineCitation']['Article']['GrantList']
            grantAgencies = set([i['Agency'] for i in grantList])
            grantIDs = set([i['GrantID'] for i in grantList])
        except:
            grantAgencies = {}
            grantIDs = {}
        
        data.loc[len(data)] = [year,month,day,title,abstract,qualifiers,grantAgencies,grantIDs]

    return data,noLoads


def getDates(x):
    """Returns series of datetimes from string date information"""
    x = x[(x.year != 'null') & (x.month != 'null') & (x.day != 'null')]
    return pd.to_datetime(x.year*10000+x.month*100+x.day,format='%Y%m%d')


def prepMethodDates(subset,query,limit=120):
    """Returns subset of df based on recency"""
    for methodName,methodData in subset.iteritems():
        t0 = datetime.now()
        methodData['elapsed'] = (getDates(methodData) - t0).dt.days
        methodData['method'] = methodName
        methodData = methodData[methodData.elapsed.abs() < limit]
        subset[methodName] = methodData
    return subset


def getCategories(df,catCol):
    """Get categorical labels within a given column"""
    try:
        categories = set(df[catCol])
    except:
        categories = set()
        for i in df[catCol]:
            categories.update(i)
    return sorted(list(categories))


def makeSet(i):
    """Converts series entry to set"""
    if i == {}:
        i = set()
    elif type(i) is not set:
        i = set([i])
    return i


def cleanGrantID(grantIDs):
    """Removes punctuation and spaces from grantIDs"""
    IDsOut = set()
    for grantID in grantIDs:
        grantID = grantID.replace(' ','-').replace('/','-')
        grantID = grantID.replace('-','')
        IDsOut.add(grantID)
    return IDsOut
                                               


cleaners = {'qualifiers':lambda x:x,
            'grantAgency':lambda x:x,
            'grantID':cleanGrantID,
            'method':lambda x: x}


def getFlatMethodDf(query,catCol,maxEntries:
    """Makes a flattened df for a given query and category column"""
    year = datetime.now().year
    subset = {query:queryToDf('"%s" AND %s:%s[dp]' % (query,year-1,year),results=maxEntries)[0]}
    subset = prepMethodDates(subset,query)
    flattened = pd.concat(subset.values())
    flattened.loc[:,catCol] = [makeSet(i) for i in flattened[catCol]]
    flattened.loc[:,catCol] = [cleaners[catCol](i) for i in flattened[catCol]]
    categories = getCategories(flattened,catCol)

    return flattened, categories
 

def getEdgeWeightDf(query,
                    catCol='qualifiers',
                    freqMin = 1,
                    limit=1000):
    """Creates an edge weight dataframe, a nodeWeight dictionary, and a categories list from a query"""
    
    flattened,categories = getFlatMethodDf(query,
                                           catCol,
                                           limit)

    nodeWeights = {i:0 for i in categories}
    edgeWeights = pd.DataFrame(0,index=categories,columns=categories)

    for index, row in flattened.iterrows():
        rowCats = row[catCol]
        for i in rowCats:
            nodeWeights[i] += 1
            for j in rowCats:
                if i != j:
                    edgeWeights.loc[i,j] += 1
    
    edgeWeights = edgeWeights.loc[(edgeWeights.sum(axis=1) >= freqMin), (edgeWeights.sum(axis=0) >= freqMin)]
    categories = edgeWeights.columns.tolist()
    return edgeWeights, nodeWeights, categories
    
    
def edgeWeightDfToNet(edgeWeights,nodeWeights,categories):
    """Creates a weighted undirected network from an edgeweight DF"""
    pubGraph = nx.Graph()
    numCats = len(categories)
    for i in range(numCats):
        catI = categories[i]
        if edgeWeights[catI].sum() != 0:
            pubGraph.add_node(catI,size=nodeWeights[catI])
        for j in range(i+1,numCats):
            catJ = categories[j]
            weight = edgeWeights.loc[catI,catJ]
            if weight != 0:
                pubGraph.add_edge(catI,catJ,weight=edgeWeights.loc[catI,catJ])
            
    return pubGraph  
    

def labelCommunityColors(graph,nodeWeights,minSize=2):
    """Detects communities in weighted, undirected graph and returns a color key for drawn nodes"""
    communities = list(greedy_modularity_communities(graph,weight='weight'))
    membership = {}
    solo = 0
    for i, community in enumerate(communities):
        for j in list(community):
            membership[j] = i
    communitySizes = {i:sum([nodeWeights[node] for node in community]) for i,community in enumerate(communities)}
    sortedCommunities = [key for key,value in sorted(communitySizes.iteritems(), key=lambda (k,v): (v,k),reverse=True)]
    rank = {community:i for i,community in enumerate(sortedCommunities)}
    nRanks = len(rank)
    palette = sb.color_palette("Paired",nRanks)
    
    def pickColor(node,solo):
        if len(communities[membership[node]]) <= 1:
            solo += 1
            return [.75,.75,.75,1.]
        else:
            return palette[membership[node]]

    nodeColors = {node:pickColor(node,solo) for node in graph.nodes()}
    

    
    return nodeColors,len(communities)-solo


def plotCategoryGraph(catGraph,
                      nodeWeights,
                      title,
                      bg=(.5,.5,.5),
                      textScale=.5):
    """Renders a viz for a category graph"""
                
    def darken(color,by=2):
        return [i/float(by) for i in color] 
                    
    pos = nx.kamada_kawai_layout(catGraph,weight='weight',scale=30)
    
    sizes=nx.get_node_attributes(catGraph,'size')
    maxSize = max(sizes.values())
    weights=nx.get_edge_attributes(catGraph,'weight')
    maxWeight = max(weights.values())
    fig, ax = plt.subplots(figsize=(20,20))
    
    nodeColors, nCommunities = labelCommunityColors(catGraph,nodeWeights)
    ax.set_title(title+"with %s Distinct Communities" % nCommunities,size=20)
    ax.set_facecolor(bg)


    posHigher = {i:[x,y+.01] for i,[x,y] in pos.iteritems()}

    for edge in catGraph.edges:
        edgeColorV = [(i+j)/2. for (i,j) in zip(nodeColors[edge[0]],nodeColors[edge[1]])]
        nx.draw_networkx_edges(catGraph,
                               pos,
                               edgelist=[edge],
                               width= (50*weights[edge]/maxWeight)**.5+1,
                               edge_color=[edgeColorV],
                               zorder = 1,
                               alpha = .3+.4*(weights[edge]/maxWeight),
                               ax=ax)


    for node in catGraph.nodes:
        drawn = nx.draw_networkx_nodes(catGraph,
                               pos,
                               node_size=(500000*sizes[node]/maxSize)**.5,
                               nodelist=[node],
                               ax=ax,
                               node_color=nodeColors[node],
                               alpha = 1,
                               zorder = 2)
        drawn.set_edgecolor(darken(nodeColors[node],1.5))
        
        nodeName = node.split()
        nodeName = ''.join(['\n'*(len(i)>3)+' '*(len(i)<4)+i for i in nodeName])
        
        if not interactive:
          nx.draw_networkx_labels(catGraph,
                                  posHigher,
                                  {node:nodeName},
                                  font_size= textScale* ((2000*sizes[node]/maxSize)**.35+10),
                                  ax=ax,
                                  font_color = darken(nodeColors[node],3),
                                  zorder = 3)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.axis('off')
    if not interactive:
      plt.show()
    else:
      mpld3.display(fig)
      
                    
def prepEdgeData(query,
                catCol='grantAgency',
                freqMin = 1,
                limit=10000,
                showNetwork = True,
                showCommunities = True):
    """One in all function generates a network viz from a pubmed query"""
                    
    edgeWeights,nodeWeights,categories = getEdgeWeightDf(query = query,
                                                        catCol=catCol,
                                                        freqMin=freqMin,
                                                        limit=limit)
        
    categoryGraph = edgeWeightDfToNet(edgeWeights,nodeWeights,categories)
    enumSubGraphs = enumerate(nx.connected_component_subgraphs(categoryGraph))

    for i,subGraph in enumSubGraphs:
      if len(subGraph) >= minSize:
        titleStr = ("''%s'' Categories for Query ''%s''\nAnd Component %s " % 
                    (catCol,query,i))
        print len(subGraph)
        ngTitleStr = "Community Network of " + titleStr
        plotCategoryGraph(subGraph,nodeWeights,ngTitleStr)
                    
    #return edgeWeights,categoryGraph,communities 
