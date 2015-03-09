#include "../gauss_elimination.hpp"
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

static const double eps(1.0e-5);
bool equal(double a, double b){
	return std::abs(a-b) < eps;
}

template<typename T>
bool verify(matrix<T> &a, matrix<T> &b, int m){
	for(int i=0; i<m; i++)
		for(int j=0; j<m; j++)
			if(!equal(a(i,j),b(i,j)))
				return false;
	return true;
}


template<typename T>
bool verify(vector<T> &a, vector<T> &b, int m){
	for(int i=0; i<m; i++)
		if(!equal(a(i),b(i)))
			return false;
	return true;
}

enum testType{ NO_PIVOT=0, NO_PERM, WITH_PERM };

struct testGaussianElimination{
	int m,n;
	matrix<double> a;
	vector<double> b;
	matrix<double> u;
	vector<double> x;

	explicit testGaussianElimination(){
		std::cin>>m>>n;
		a=matrix<double>(m,n);
		b=vector<double>(m);
		u=matrix<double>(m,n);
		x=vector<double>(m);
		for(int i=0; i<m; i++)
			for(int j=0; j<n; j++)
				std::cin>>a(i, j);
		for(int i=0; i<m; i++) std::cin>>b(i);
		
		for(int i=0; i<m; i++)
			for(int j=0; j<n; j++)
				std::cin>>u(i, j);
		for(int i=0; i<m; i++) std::cin>>x(i);
	}

	void test(testType t){
		switch(t){
			case NO_PIVOT:
				gaussian_elimination_no_pivot(a,b);
				if(verify(a,u,m))
					std::cout<<"Passed gaussian_elimination_no_pivot"<<std::endl;
				else std::cout<<"Failed gaussian_elimination_no_pivot"<<std::endl;
							  //<<a<<std::endl<<u<<std::endl;
				try{
					gauss_substitute(a,b);
					if(verify(b,x,m))
						std::cout<<"Passed solver_without_pivot"<<std::endl;
					else std::cout<<"Failed solver_without_pivot"<<std::endl
								  <<b<<std::endl<<x<<std::endl;
					 
				}
				catch(...){
					std::cout<<"Exception in solver_without_pivot (instability)"
							 <<std::endl;
				}
				break;
			case NO_PERM:
				gaussian_elimination(a,b);
				if(verify(a,u,m))
					std::cout<<"Passed gaussian_elimination"<<std::endl;
				else std::cout<<"Failed gaussian_elimination"<<std::endl
							  <<a<<std::endl<<u<<std::endl;
				gauss_substitute(a,b);
				if(verify(b,x,m))
					std::cout<<"Passed solver"<<std::endl;
				else std::cout<<"Failed solver"<<std::endl
							  <<b<<std::endl<<x<<std::endl;
				break;
			default:
				permutation_matrix<> pm(m);
				gaussian_elimination(a,pm,b);
				if(verify(a,u,m))
					std::cout<<"Passed gaussian_elimination_with_permutation_matrix"<<std::endl;
				else std::cout<<"Failed gaussian_elimination_with_permutation_matrix"<<std::endl
							  <<a<<std::endl<<u<<std::endl;
				gauss_substitute(a,b);
				if(verify(b,x,m))
					std::cout<<"Passed solve_with_permutation_matrix"<<std::endl;
				else std::cout<<"Failed solver_with_permutation_matrix"<<std::endl
							  <<b<<std::endl<<x<<std::endl;
				break;
		}
	}

};


int main(){
	for(int i=0; i<200; i++){
		for(int j=0; j<3; j++){
			std::cout<<"TEST "<<i<<"."<<j<<std::endl<<"=================="<<std::endl;
			testGaussianElimination t = testGaussianElimination();
			t.test(j);
			std::cout<<std::endl<<std::endl;
		}
	}
	return 0;
}




/*int main(){
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
}*/
