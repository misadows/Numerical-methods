#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>

using namespace std;

double fun(int i, int j){
	if(i==0) return 1;
	else return 1.0/(i+1+j+1-1);
}

double rand_fun(int i, int j){
	return rand()%6+1;
}

double fun_x(int i){
	if(rand()%2) return 1;
	else return -1;
}

int iden(int i){
	return i;
}

void print_vector(vector<double> v){
	for(int i=0; i<v.size(); i++){
		cout<<v[i]<<" ";
	}
	cout<<endl;
}

void print_matrix(vector<vector<double> > v){
	for(int i=0; i<v.size(); i++){
		print_vector(v[i]);
	}
	cout<<endl;
}

vector<double> multiply_matrices(vector<vector<double> > m1, vector<double> m2){
	vector<double> result(m1.size());
	for(int i=0; i<m1.size(); i++){
		if(m1[i].size()!=m2.size()) throw "Cannot multiply matrices";
		double sum=0;
		for(int j=0; j<m2.size(); j++){
			sum+=m1[i][j]*m2[j];
		}
		result[i]=sum;
	}

	return result;
}

vector<vector<double> > multiply_matrices(vector<vector<double> > m1, vector<vector<double> > m2){
	vector<vector<double> > result(m1.size());
	if(m2.size()<1){
		return result;
	}
	int width=m2[0].size();
	for(int i=0; i<m2.size(); i++){
		if(m2[i].size()!=width) throw "Non square matrix";
	}
	for(int i=0; i<m1.size(); i++){//row of m1
		result[i]=vector<double>(width);
		for(int j=0; j<width; j++){//column of m2
			if(m1[i].size()!=m2.size()) throw "Cannot multiply matrices";
			double sum=0;
			for(int k=0; k<m2.size(); k++){//elements in i row of m1 and j column of m2
				sum+=m1[i][k]*m2[k][j];
			}
			result[i][j]=sum;
		}
	}
	return result;
}

vector<vector<double> > create_matrix(int n, double (*fun)(int,int)){
	vector<vector<double> > matrix;
	matrix.resize(n, vector<double>(n, 0));

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			matrix[i][j]=fun(i,j);
		}
	}
	return matrix;
}

vector<vector<double> > create_matrix(int i, int j, double (*fun)(int,int)){
	vector<vector<double> > matrix;
	matrix.resize(i);

	for(int it=0; it<i; it++){
		matrix[it]=vector<double>(j);
		for(int jt=0; jt<j; jt++){
			matrix[it][jt]=fun(it,jt);
		}
	}
	return matrix;
}

vector<double> create_one_dim_matrix(int n, double (*fun)(int)){
	vector<double> matrix(n);
	for(int i=0; i<n; i++){
		matrix[i]=fun(i);
	}
	return matrix;
}

vector<int> create_permutation(vector<vector<double> > matrix){
	//matrix must be square
	vector<int> per(matrix.size());
	vector<double> maxes(matrix.size());
	for(int i=0; i<matrix.size(); i++){
		per[i]=i;
		maxes[i]=max(abs(*max_element(matrix[i].begin(), matrix[i].end())), abs(*min_element(matrix[i].begin(), matrix[i].end())));
	}
	for(int k=0; k<matrix.size(); k++){
		int max_pos=k;
		for(int w=k; w<matrix.size(); w++){
			if(abs(matrix[per[max_pos]][k]/maxes[per[max_pos]]) < abs(matrix[per[w]][k]/maxes[per[w]])) max_pos=w;
		}
		swap(per[max_pos],per[k]);
	}

	return per;
}

vector<int> create_permutation_iden(vector<vector<double> > matrix){
	//matrix must be square
	vector<int> per(matrix.size());
	for(int i=0; i<matrix.size(); i++){
		per[i]=i;
	}
	return per;
}

void gauss_elimination(vector<vector<double> > &matrix, vector<double> &y ,vector<int> p){
	for(int k=0; k<matrix.size()-1; k++){
		for(int i=k+1; i<matrix.size(); i++){
			double z=matrix[p[i]][k]/matrix[p[k]][k];
			matrix[p[i]][k]=0;
			for(int j=k+1; j<matrix.size(); j++){
				matrix[p[i]][j] -= z*matrix[p[k]][j];
			}
			y[p[i]] -= z*y[p[k]];
		}
	}
}

vector<double> calculate_x(vector<vector<double> > matrix, vector<double> y, vector<int> p){
	vector<double> x(matrix.size(), 0);
	for(int i=matrix.size()-1; i>=0; i--){
		double sum=0;
		for(int j=i+1; j<matrix.size(); j++){
			sum+=matrix[p[i]][j]*x[p[j]];
		}
		x[p[i]]=(y[p[i]]-sum)/matrix[p[i]][i];
	}
	return x;
}

int main(){
	// srand(time(NULL));
	// int n=15;
	// vector<vector<double> > matrix = create_matrix(n, fun);
	// vector<double> x = create_one_dim_matrix(n, fun_x);
	// vector<double> y = multiply_matrices(matrix, x);
	// //vector<int> p = create_permutation(matrix);
	// vector<int> p = create_permutation_iden(matrix);

	// cout<<endl;
	// for(int i=0; i<matrix.size(); i++) {
	//  	print_vector(matrix[i]);
	// }
	// cout<<endl;
	// gauss_elimination(matrix, y, p);
	// print_vector(x);
	// vector<double> res_x = calculate_x(matrix, y, p);
	// for(int i=0; i<matrix.size(); i++){
	// 	cout<<res_x[p[i]]<<" ";
	// }
	// cout<<"\n\n";

	//  for(int i=0; i<matrix.size(); i++) {
	//  	print_vector(matrix[p[i]]);
	// }
	try{
		vector<vector<double> > v1=create_matrix(3,2, rand_fun);
		vector<vector<double> > v2=create_matrix(2,4, rand_fun);
		vector<vector<double> > res=multiply_matrices(v1,v2);
		print_matrix(v1);
		print_matrix(v2);
		print_matrix(res);
	}
	catch(char const* error){
		cout<<error<<endl;
	}
	



	return 0;
}