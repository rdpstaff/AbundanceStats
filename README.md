## RDP AbundanceStats

### Intro
This package computes Shannon and Chao diversity indices, jaccard or sorensen abundance stats and rarefaction curve. See more details on http://rdp.cme.msu.edu/tutorials/cluster/RDPtutorial_CLUST-RESULTS.html.

### Setup


### Usage

* Compute Shannon and Chao1 indices

		java -cp /path/to/AbundanceStats.jar edu/msu/cme/rdp/abundstats/ShannonChao all_seq_complete.clust shannchao

* Compute jaccard or sorensen abundance stats

		java -jar /path/toAbundanceStats.jar
		usage: Main [options] <cluster file>
 		-j,--jaccard              Compute jaccard abundance stats
 		-l,--lower-cutoff <arg>   Lowest cutoff in the cluster file to compute stats for
 		-r,--result-dir <arg>     Directory to put the result files in (default=.)
 		-R,--R-location <arg>     Triggers the R plotter subsystem, provide the path to the R command
 		-s,--sorensen             Compute sorensen abundance stats
 		-t,--otu-table            input file is an otu table, not rdp cluster file
 		-u,--upper-cutoff <arg>   Highest cutoff in the cluster file to compute stats for
 		
 	An example command to calculate jaccard distance matrix, perform UPGMA clustering and plot dendrogram in R environment. 
 	
 		java -jar /path/to/AbundanceStats.jar -jaccard -R /usr/bin/R all_seq_complete.clust

* Calculate rarefaction curve

		java -cp /path/to/AbundanceStats.jar edu/msu/cme/rdp/rarefaction/Rarefaction all_seq_complete.clust rarefaction plot
