#include "../gauss_elimination.hpp"
#include <boost/numeric/ublas/io.hpp>

int main(){
	using namespace boost::numeric::ublas;
	int m, n, temp;
	std::cin>>m>> n;
	matrix<double> mat(m, n);
	vector<double> v(m);
	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			std::cin>>mat(i, j);
	for(int i=0; i<m; i++) std::cin>>v(i);

	std::cout<<mat<<std::endl<<v<<std::endl;

//	permutation_matrix<> pm(mat.size1());
	gaussian_elimination(mat, v);
	//swap_rows(pm, mat);
//	row(mat,1).swap(row(mat,2));
//	std::swap(v(1),v(2));
	std::cout<<mat<<std::endl;
	std::cout<<v<<std::endl;
	//std::cout<<pm<<std::endl;
//	vector<double> v(n);
	gauss_substitute(mat, v);
	std::cout<<v<<std::endl;
	return 0;
}
