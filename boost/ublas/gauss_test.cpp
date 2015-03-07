#include "gauss_elimination.hpp"
#include <boost/numeric/ublas/io.hpp>

int main(){
	using namespace boost::numeric::ublas;
	int m, n, temp;
	std::cin>>m>> n;
	matrix<double> mat(m, n);
	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			std::cin>>mat(i, j);

	gaussian_elimination(mat);
	std::cout<<mat<<std::endl;
	return 0;
}
