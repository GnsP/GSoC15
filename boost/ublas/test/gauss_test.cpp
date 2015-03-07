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

	permutation_matrix<> pm(mat.size1());
	gaussian_elimination(mat, pm);
	swap_rows(pm, mat);
	row(mat,1).swap(row(mat,2));
	std::cout<<mat<<std::endl;
	std::cout<<pm<<std::endl;
	vector<double> v(n);
	lu_substitute(mat, pm, v);
	std::cout<<v<<std::endl;
	return 0;
}
