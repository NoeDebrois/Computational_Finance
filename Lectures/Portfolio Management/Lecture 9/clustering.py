# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 14:43:44 2024

@author: ginevra.angelini
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

def get_performance_metrics(ce):
    # Compute Performance Metrics
    aret = ann_ret(ce)
    avol = ann_vol(ce.pct_change().dropna())
    sr = sharpe_ratio(ce)
    mxdd = max_drawdown(ce)
    cr = calmar_ratio(ce)
    df_perf = pd.DataFrame([aret,avol, sr, mxdd, cr], index = ['AnnRet', 'AnnVol', 'Sharpe', 'MaxDD', 'Calmar'], columns = ['Metrics'] )
    return df_perf
    
## Ex1 : Create a strategy with clustering algorithms: K-means, Hierarchical Clustering
##       Universe: sectorial and factorial index
##       Idea: grouping index with similar characteristics
# Read data
path_ = your_path+'\\prices.xlsx'
data = pd.read_excel(path_)
data.index = data.iloc[:,0]
data = data.iloc[:,1:]
data = data.loc[data.index >= pd.to_datetime('2010-01-01')]
data.plot()
returns = data.pct_change().dropna()

# Step 1: Compute Correlation Matrix
correlation_matrix = returns.corr()
sns.heatmap(correlation_matrix, annot=False, cmap="coolwarm", square=True)
plt.title("Matrice di correlazione")
plt.show()


# Step 2: Features Computation

features = pd.DataFrame({"mean_return": returns.mean(),  # Exp Return
                         "volatility" : returns.std(),   # Volatility
                         "skewness"   : returns.skew(),  # Asymmetry
                         "kurtosis"   : returns.kurt()   # Kurtosis
                         })

# Step 3: Scaling Features
scaler = StandardScaler()
normalized_features = scaler.fit_transform(features)

# Step 4: Clustering computation
# Algorithm 1: K-Means Clustering
kmeans = KMeans(n_clusters = 4, random_state=42)
clusters_kmeans = kmeans.fit_predict(normalized_features)
features["Cluster_KMeans"] = clusters_kmeans

# Algorithm 2: Hierarchical Clustering
linkage_matrix = linkage(normalized_features, method="ward")
dendrogram(linkage_matrix, labels=features.index, leaf_rotation=90)
plt.title("Dendrogram (Hierarchical Clustering)")
plt.show()

clusters_hierarchical = fcluster(linkage_matrix, t=4, criterion="maxclust")
features["Cluster_Hierarchical"] = clusters_hierarchical

# Step 5: Building Strategies -> Equal weight for each asset in cluster
weights = pd.DataFrame(index=returns.columns)  
for cluster in features["Cluster_KMeans"].unique():
    cluster_assets = features[features["Cluster_KMeans"] == cluster].index
    weights[cluster] = [1 / len(cluster_assets) if asset in cluster_assets else 0 for asset in returns.columns]

weights = weights.fillna(0)

# Step 6: Portfolio/ Equity computation and evaluation
portfolio_returns = returns.iloc[:, :] @ weights  
cumulative_returns = (1 + portfolio_returns).cumprod()

# Benchmark
ce_ew = (1 + returns.mean(1)).cumprod()
#plt.figure(figsize=(10, 6))
cumulative_returns.plot()
ce_ew.plot(c = 'k')
plt.title("Cluster Portfolio Performance")
plt.xlabel("Date")
plt.ylabel("Cumulative Returns")
plt.legend(title="Cluster")
plt.grid(True)
plt.show()

# Performance metrics
list_tbl = []
for col in cumulative_returns.columns:
    perf_metrics = get_performance_metrics(cumulative_returns.loc[:, col])
    perf_metrics.columns = [col]
    list_tbl.append(perf_metrics)

perf_table = pd.concat(list_tbl,1)
perf_ew = get_performance_metrics(ce_ew)

### In the first there are some imprecisions and in general we can make some improvements
## Ex 2: Enhancement of Ex1

def max_drawdown(series):
    cumulative = (1 + series).cumprod()
    peak = cumulative.cummax()
    drawdown = (cumulative - peak) / peak
    return drawdown.min()

from sklearn.model_selection import train_test_split
returns = data.pct_change().dropna()

# Imprecision: There is not train- test splitting, this can lead to overfitting
# 1) Train-Test splitting
train_returns, test_returns = train_test_split(returns, test_size=0.2, shuffle=False)

# 2) New Features: Enhance the feature set (computed on the TRAIN SET)
benchmark = train_returns.mean(1)
features_train = pd.DataFrame({"mean_return" : train_returns.mean(), 
                               "volatility"  : train_returns.std(),
                               "mean_corr"   : train_returns.corr().mean(), # Mean corr between index
                               "max_drawdown": train_returns.apply(max_drawdown), # Max Drawdown
                               "sharpe_ratio": train_returns.mean()/train_returns.std(),
                               "stability"   : 1/train_returns.rolling(window=30).std().mean(), # return stability
                               "beta"        : train_returns.apply(lambda x: np.cov(x, benchmark)[0, 1] / np.var(benchmark)) # sensibility of an index respect to a benchmark
                               })

# 3): Features Scaling -> SCALER IS FITTED ON TRAIN SET AND THEN APPLIED TO THE TEST SET
scaler = StandardScaler()
normalized_train_features = scaler.fit_transform(features_train)
normalized_features_df = pd.DataFrame(normalized_train_features, columns=features_train.columns, index=features_train.index)

#normalized_features_df.plot() #check stazionarità

## ------------------------------------  Feature Analysis ------------------------------------ ##
## Pairwise analysis
sns.pairplot(pd.DataFrame(normalized_train_features, columns=features_train.columns), diag_kind='kde', corner=True)
plt.suptitle("Scatter Matrix of Normalized Features", y = 1.02)
plt.show()

## Analysis of intra-feature variance
# The variance of a normalized feature tell us if the feature has a wide distribution (Feature with low variance (< 0.1) could be Uninformative
feature_variance = pd.DataFrame({"Feature": features_train.columns,
                                 "Variance": normalized_train_features.var(axis=0)
                                 })
print("Intra-feature Variance:")
print(feature_variance)

## Analysis of feature correlation
# If 2 features are higly correlated (absolute value > 0.9) could be redoundant
correlation_matrix = pd.DataFrame(normalized_train_features, columns=features_train.columns).corr()
plt.figure(figsize=(10, 8))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f", vmin=-1, vmax=1)
plt.title("Features correlation matrix")
plt.show()
## ------------------------------------------------------------------------------------------- ##

## 4) Which is the optimal number of clusters? 
# Elbow 
inertia = []
k_values = range(1, 16)
for k in k_values:
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(normalized_train_features)
    inertia.append(kmeans.inertia_)

plt.figure(figsize=(8, 5))
plt.plot(k_values, inertia, marker='o')
plt.title("Elbow Method: searching for optimal number of clusters")
plt.xlabel("Number of clusters (k)")
plt.ylabel("Inertia")
plt.grid(True)
plt.show()

# 5) Clustering with the optimal number of clusters (result from the elbow method)
optimal_k = 4  
kmeans = KMeans(n_clusters=optimal_k, random_state=42)
features_train["Cluster"] = kmeans.fit_predict(normalized_train_features)


# Let's see the pairplot for different clusters
sns.pairplot(features_train,
             vars = features_train.columns[:-1],  # Exclusion of columns 'cluster' 
             hue = 'Cluster',                     # Coloring dots based on cluster
             palette = 'Set2',                    # color palette
             diag_kind = 'kde',                   # Distribution on diagonal
             corner = True                        # Showing only the lower part of the matrix
             )
plt.suptitle("Scatter Matrix of normalized features", y = 1.02)
plt.show()


# 6) Application on the TEST SET
benchmark_test = test_returns.mean(1)
normalized_test_features = scaler.transform(pd.DataFrame({"mean_return" : test_returns.mean(),
                                                          "volatility"  : test_returns.std(),
                                                          "mean_corr"   : test_returns.corr().mean(),
                                                          "max_drawdown": test_returns.apply(max_drawdown),
                                                          "sharpe_ratio": test_returns.mean()/test_returns.std(),
                                                          "stability"   : 1/test_returns.rolling(window=30).std().mean(),
                                                          "beta"        : test_returns.apply(lambda x: np.cov(x, benchmark_test)[0, 1] / np.var(benchmark_test))
                                                          }))
test_clusters = kmeans.predict(normalized_test_features)

# 7) Building Portfolios: equal weights in each cluster
weights_train = pd.DataFrame(index=train_returns.columns)
for cluster in features_train["Cluster"].unique():
    cluster_assets = features_train[features_train["Cluster"] == cluster].index
    weights_train[cluster] = [1 / len(cluster_assets) if asset in cluster_assets else 0 for asset in train_returns.columns]

weights_train = weights_train.fillna(0)

# 8) Portfolio evaluation
portfolio_returns_test = test_returns @ weights_train
cumulative_returns_test = (1 + portfolio_returns_test).cumprod()
cumulative_returns_test = 100*cumulative_returns_test/cumulative_returns_test.iloc[0,:]

# Benchmark
ce_ew = (1 + returns.mean(1)).cumprod()
ce_ew = ce_ew.loc[ce_ew.index >= cumulative_returns_test.index[0]]
ce_ew = 100*ce_ew/ce_ew.iloc[0]

#plt.figure(figsize=(10, 6))
cumulative_returns_test.plot()
ce_ew.plot(c = 'k')
plt.title("Clustering Portfolio Performances (Test Set)")
plt.xlabel("Date")
plt.ylabel("Cumulative Returns")
plt.legend(title="Cluster")
plt.grid(True)
plt.show()

# Performance Metrics
list_tbl = []
for col in cumulative_returns_test.columns:
    perf_metrics = get_performance_metrics(cumulative_returns_test.loc[:, col])
    perf_metrics.columns = [col]
    list_tbl.append(perf_metrics)

perf_table = pd.concat(list_tbl,1)
perf_ew = get_performance_metrics(ce_ew)


## EX3: Identify Market regimes using Kmeans algorithm on spx drawdown data and find the best portfolio in each regime
##      Portfolios will be built with sectorial and factorial indexes
def get_top_sectors_for_cluster(sector_performance, cluster, N=3):
    """
    Select the top N sectors for a given cluster based on their average performance.
    Inputs:
        - sector_performance: A DataFrame containing the average performance of sectors for each cluster.
        - cluster: The cluster for which to select the top sectors (e.g., 0, 1, 2, etc.)
        - N: Number of sectors to select.
    Output:
        A list of the top sectors for the specified cluster
    """
    # Sort the sectors based on the performance of the selected cluster (in descending order)
    sorted_sectors = sector_performance.loc[cluster].sort_values(ascending=False)
    # Selecting first N sectors
    top_sectors = sorted_sectors.head(N).index.tolist()
    return top_sectors

# Step 1: collect data
# Read data
path_spx = your_path+'\\spx_price.p'
df = pd.read_pickle(path_spx)
df.index.name = 'Date'

# Calculate the daily returns of the S&P 500
df['SPX_returns'] = df['SPX Index'].pct_change()
df = df.loc[df.index >= pd.to_datetime('2005-01-01')]


# Compute rolling drawdown of SPX Index
window = 252
df['SPX_max'] = df['SPX Index'].rolling(window=window, min_periods=1).max()
df['drawdown'] = df['SPX Index'] / df['SPX_max'] - 1  # Calcola il drawdown

# Plot drawdowns
plt.figure(figsize=(10, 6))
plt.plot(df.index, df['drawdown'], label='Drawdown S&P 500')
plt.title("SPX Drawdown")
plt.xlabel("Data")
plt.ylabel("Drawdown")
plt.legend()
plt.show()

# Step 2: drawdown Clustering
drawdown_data = df[['drawdown']].dropna()  # Remove NaN

# Scaling Data
scaler = StandardScaler()
drawdown_scaled = scaler.fit_transform(drawdown_data)

## Elbow method
inertia = []
k_values = range(1, 10)
for k in k_values:
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(drawdown_data)
    inertia.append(kmeans.inertia_)

plt.figure(figsize=(8, 5))
plt.plot(k_values, inertia, marker='o')
plt.title("Elbow Method: searching for optimal number of clusters")
plt.xlabel("Number of clusters (k)")
plt.ylabel("Inertia")
plt.grid(True)
plt.show()
##

# Kmeans clustering
kmeans = KMeans(n_clusters = 3, random_state = 42)  
df['drawdown_cluster'] = kmeans.fit_predict(drawdown_scaled)

# Plot clustering
plt.figure(figsize=(10, 6))
sns.lineplot(data = df, x = 'Date', y = 'drawdown', hue = 'drawdown_cluster', palette = 'Set2', linewidth = 2)
plt.title("Drawdown Clustering of SPX")
plt.xlabel("Data")
plt.ylabel("Drawdown")
plt.legend(title='Cluster')
plt.show()

# Step 3: Analyze the performance of sector indices for each market regime
path_ = your_path+'\\prices.xlsx'
data = pd.read_excel(path_)
data.index = data.iloc[:,0]
data = data.iloc[:,1:]
data = data.loc[data.index >= df.index[0]]
sector_columns = data.columns  
sector_returns = data[sector_columns].pct_change() # compute index returns

# Adding cluster
sector_returns['drawdown_cluster'] = df['drawdown_cluster']

# Calculate the average performance of sectors for each cluster
sector_performance = sector_returns.groupby('drawdown_cluster').mean()

# Plot index performance for each cluster
plt.figure(figsize=(10, 6))
sector_performance.T.plot(kind='bar', figsize=(10, 6))
plt.title("Index Performance by Market Regime (Drawdown Clusters)")
plt.xlabel("Index")
plt.ylabel("Mean Return")
plt.legend(title='Cluster')
plt.show()

# Step 4: Portfolio Construction Based on Clusters: we choose top N indexes for each cluster
top_N = 3 # selecting top 3 sectors
portfolio_allocation = {}
for cl in [0,1,2]:
    cluster_to_analyze = cl
    best_sectors_for_cluster = get_top_sectors_for_cluster(sector_performance, cluster_to_analyze, N=top_N)
    portfolio_allocation[cl] = best_sectors_for_cluster


# Assign portfolio weights based on the current cluster
def get_portfolio_weights(cluster, dt):
    sectors_to_invest = portfolio_allocation.get(cluster, [])
    weights = pd.DataFrame(columns = data.columns, index = [dt])
    for sector in sectors_to_invest:
        weights[sector] = 1 / len(sectors_to_invest)  # Pesare equamente tra i settori scelti
    return weights.fillna(0)

# Assign portfolio weights based on the current cluster, but if we are in the cluster of "market crisis" we put exposition = 0
def get_portfolio_weights_null(cluster, dt):
    if cluster != 2:
        sectors_to_invest = portfolio_allocation.get(cluster, [])
        weights = pd.DataFrame(columns = data.columns, index = [dt])
        for sector in sectors_to_invest:
            weights[sector] = 1 / len(sectors_to_invest)  # Pesare equamente tra i settori scelti
    else:
        weights = pd.DataFrame(np.zeros((1,len(data.columns))), columns = data.columns, index = [dt])
    return weights.fillna(0)


l_weights = []
for idx in data.index:
    print(idx)
    sample_cluster = int(df.loc[idx]['drawdown_cluster']) 
    weights = get_portfolio_weights_null(sample_cluster, idx)
    #weights = get_portfolio_weights(sample_cluster, idx)
    l_weights.append(weights)

weights_df = pd.concat(l_weights)
returns = data.pct_change()
returns = returns.loc[:, weights_df.columns]

# Compute Equity
ret_strat = pd.DataFrame(weights_df.iloc[0:-1].values*returns.iloc[1:,:].values, index = returns.iloc[1:,:].index).sum(1)
eq = (ret_strat+1).cumprod()
eq = 100*eq/eq.iloc[0]

spx_price = df.loc[:, 'SPX Index']
spx_price =spx_price.reindex(eq.index)
spx_price = 100*spx_price/spx_price.iloc[0]

# PLot
eq.plot(c = 'limegreen',label = 'Portfolio', legend = True)
spx_price.plot(c = 'navy', label = 'SPX', legend = True)
plt.grid()

# Compute Performance Metrics
perf_eq = get_performance_metrics(eq)
perf_spx = get_performance_metrics(spx_price)