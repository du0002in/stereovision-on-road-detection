#include <vector>
#include <cmath>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/foreach.hpp>
#include "distance.h"

namespace Clustering{

	// a single point is made up of vector of doubles
	typedef boost::numeric::ublas::vector<double> Point;
	typedef std::vector<Point> Points;

	typedef unsigned int ClusterId;
	typedef unsigned int PointId;	

	// a cluster is a vector of pointid
	typedef std::vector<PointId> Cluster;
	// a set of Neighbors is a vector of pointid
	typedef std::vector<PointId> Neighbors;


	void randomInit	(Points & ps, unsigned int dims = 5, 
						unsigned int num_points = 10);

	class Clusters
	{
	public:
		Clusters (Points & ps) : _ps(ps) 
		{
			_pointId_to_clusterId.resize(_ps.size(), 0);
		};

		// assign each point to a new cluster
		void uniformPartition();

		// compute similarity
		template <typename Distance_type>
		void computeSimilarity(Distance_type & d)
		{
			unsigned int size = _ps.size();
			_sim.resize(size, size, false);
			for (unsigned int i=0; i < size; i++)
			{
				for (unsigned int j=i+1; j < size; j++)
				{
					_sim(j, i) = _sim(i, j) = d.similarity(_ps[i], _ps[j]);
	//				std::cout << "(" << i << ", " << j << ")=" << _sim(i, j) << " ";
				}
				std::cout << std::endl;
			}
		};
	
		//
		// findNeighbors(PointId pid, double threshold)
		//
		// this can be implemented with reduced complexity by using R+trees
		//
		Neighbors findNeighbors(PointId pid, double threshold);
	
	protected:
		// the collection of points we are working on
		Points& _ps;
		
		// mapping point_id -> clusterId
		std::vector<ClusterId> _pointId_to_clusterId;

		// the collection of clusters
		std::vector<Cluster> _clusters;

		// simarity_matrix
		boost::numeric::ublas::matrix<double> _sim;

		friend	
			std::ostream& operator << (std::ostream& o, const Clusters& c);
	};

}