#include "dbscan.h"
#include "distance.h"
#include <iostream>

int main(){

	using namespace Metrics;

	Clustering::Points ps;

	// random init points dataset (dims, points)
	Clustering::randomInit(ps, 10, 100);   

	// init: sim threshold, minPts
	Clustering::DBSCAN clusters(ps, 0.8, 2); 

	// uniform distribution dataset
	//	clusters.uniformPartition();          

	// build similarity  matrix
	Distance<Cosine<Clustering::Point> > d;
	clusters.computeSimilarity(d);     

	// run clustering
	clusters.run_cluster();

	std::cout << clusters;

	int test;
	std::cin >> test;

}