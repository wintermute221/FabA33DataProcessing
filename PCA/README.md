# PCA

 For the aim to use the structural clustering to then identify residues that account for the structural clusters

 ## Rational

It seems that it makes more sense to cluster residues one by one across conditions. In such way different residues might end up in different clusters, but this exactly reveals how residues respond to external conditions differently. So I would probably only focus on six replicas for one condition once a time. 


**Also I need to be aware of this step actually constitutes two sub-steps. The first step is PCA, which reduces the dimensionality of the data and provide clustering signpost to feed in the next step. The second step is doing classification, based on the number of groups decided by the PCA step, by k-means.** 
