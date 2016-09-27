#include "clusters.h"
#include <boost/foreach.hpp>

namespace Clustering
{
	void randomInit	(Points & ps, unsigned int dims, 
					unsigned int num_points)
	{
		for (unsigned int j = 0; j < num_points; j++)
		{
			Point p(dims);
			for (unsigned int i = 0; i < dims; i++)	
			{
				p(i) = (-1.0 + rand() * (2.0) / RAND_MAX);
	//			std::cout << p(i) << ' ';
			}
			ps.push_back(p);
	//		std::cout << std::endl;
		}
	}

	// assign each point to a new cluster
	void Clusters::uniformPartition()
	{
		PointId pid = 0;
		ClusterId cid = 0;
		BOOST_FOREACH(Point p, _ps)
		{
			// create a new cluster for this current point
			Cluster c;
			c.push_back(pid++);
			_clusters.push_back(c);
			_pointId_to_clusterId.push_back(cid++);			
		}
	}

	// findNeighbors(PointId pid, double threshold) 	
	Neighbors Clusters::findNeighbors(PointId pid, double threshold)
	{
		Neighbors ne;

		for (unsigned int j=0; j < _sim.size1(); j++)
		{
			if 	((pid != j ) && (_sim(pid, j)) > threshold)
			{
				ne.push_back(j);
				//std::cout << "sim(" << pid  << "," << j << ")" << _sim(pid, j) << ">" << threshold << std::endl;
			}
		}
		return ne;
	};


	// single point output
	std::ostream& operator<<(std::ostream& o, const Point& p)
	{
		o << "{ ";
		BOOST_FOREACH(Point::value_type x, p)
		{
			o << " " << x;
		}
		o << " }, ";

		return o;
	}

	// single cluster output
	std::ostream& operator<<(std::ostream& o, const Cluster& c)
	{
		o << "[ ";
		BOOST_FOREACH(PointId pid, c)
		{
			o << " " << pid;
		}
		o << " ]";

		return o;
	}

	// clusters output
	std::ostream& operator<<(std::ostream& o, const Clusters& cs)
	{
		ClusterId cid = 1;
		BOOST_FOREACH(Cluster c, cs._clusters)
		{
			o << "c(" << cid++ << ")=";
			
			BOOST_FOREACH(PointId pid, c)
			{
				o << cs._ps[pid];
			}
			o << std::endl;
		}
		return o;
	}


};
