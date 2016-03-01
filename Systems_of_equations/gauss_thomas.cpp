#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <ctime>
#define FLOAT_TYPE double

using namespace std;

FLOAT_TYPE fun1(int i, int j){
	if(i==0) return 1;
	else return 1.0/(i+1+j+1-1);
}

FLOAT_TYPE fun2(int i, int j){
	if(j>=i) return 2.0*(i+1)/(FLOAT_TYPE)(j+1);
	else return 2.0*(j+1)/(FLOAT_TYPE)(i+1);
}

FLOAT_TYPE fun3(int i, int j){
	int m=3;
	int k=6;
	if(i==j) return k;
	if(j==(i+1)) return 1.0/(i+1+m);
	if(i>0 && j==(i-1)) return (FLOAT_TYPE)k/(i+1+m+1);
	return 0;
}

FLOAT_TYPE fun4(int i, int j){
	int m=3;
	int k=6;
	if(i==j) return (FLOAT_TYPE)m*(i+1.0)*(-1.0)-k;
	if(j==(i+1)) return (FLOAT_TYPE)(i+1);
	if(i>0 && j==(i-1)) return (FLOAT_TYPE)m/(i+1.0);
	return 0;
}

FLOAT_TYPE fun_x(int i){
	if(rand()%2) return 1;
	else return -1;
}


vector<FLOAT_TYPE> multiply_matrices(vector<vector<FLOAT_TYPE> > m1, vector<FLOAT_TYPE> m2){
	vector<FLOAT_TYPE> result(m1.size());
	for(int i=0; i<m1.size(); i++){
		if(m1[i].size()!=m2.size()) throw "Cannot multiply matrices";
		FLOAT_TYPE sum=0;
		for(int j=0; j<m2.size(); j++){
			sum+=m1[i][j]*m2[j];
		}
		result[i]=sum;
	}

	return result;
}

vector<vector<FLOAT_TYPE> > create_matrix(int n, FLOAT_TYPE (*fun)(int,int)){
	vector<vector<FLOAT_TYPE> > matrix;
	matrix.resize(n, vector<FLOAT_TYPE>(n, 0));

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			matrix[i][j]=fun(i,j);
		}
	}
	return matrix;
}

vector<FLOAT_TYPE> create_one_dim_matrix(int n, FLOAT_TYPE (*fun)(int)){
	vector<FLOAT_TYPE> matrix(n);
	for(int i=0; i<n; i++){
		matrix[i]=fun(i);
	}
	return matrix;
}



vector<int> create_permutation_iden(vector<vector<FLOAT_TYPE> > matrix){
	//matrix must be square
	vector<int> per(matrix.size());
	for(int i=0; i<matrix.size(); i++){
		per[i]=i;
	}
	return per;
}

void gauss_elimination(vector<vector<FLOAT_TYPE> > &matrix, vector<FLOAT_TYPE> &y ,vector<int> p){
	for(int k=0; k<matrix.size()-1; k++){
		for(int i=k+1; i<matrix.size(); i++){
			FLOAT_TYPE z=matrix[p[i]][k]/matrix[p[k]][k];
			matrix[p[i]][k]=0;
			for(int j=k+1; j<matrix.size(); j++){
				matrix[p[i]][j] -= z*matrix[p[k]][j];
			}
			y[p[i]] -= z*y[p[k]];
		}
	}
}

vector<FLOAT_TYPE> calculate_x(vector<vector<FLOAT_TYPE> > matrix, vector<FLOAT_TYPE> y, vector<int> p){
	vector<FLOAT_TYPE> x(matrix.size(), 0);
	for(int i=matrix.size()-1; i>=0; i--){
		FLOAT_TYPE sum=0;
		for(int j=i+1; j<matrix.size(); j++){
			sum+=matrix[p[i]][j]*x[p[j]];
		}
		x[p[i]]=(y[p[i]]-sum)/matrix[p[i]][i];
	}
	return x;
}

void print_vector(vector<FLOAT_TYPE> v){
	for(int i=0; i<v.size(); i++){
		cout<<setw(6)<<v[i]<<" ";
	}
	cout<<endl;
}

void print_matrix(vector<vector<FLOAT_TYPE> > v){
	for(int i=0; i<v.size(); i++){
		print_vector(v[i]);
	}
	cout<<endl;
}

FLOAT_TYPE calculate_euclidan_distance(vector<FLOAT_TYPE> v1, vector<FLOAT_TYPE> v2){
	FLOAT_TYPE sum=0;
	if(v1.size()!=v2.size()) throw "Cannot calculate euclidan distance";
	for(int i=0; i<v1.size(); i++){
		sum+=pow(v1[i]-v2[i],2);
	}
	return sqrt(sum);
}

void test_one_matrix(int n, FLOAT_TYPE (*fun)(int,int)){
	vector<vector<FLOAT_TYPE> > matrix = create_matrix(n, fun);
	vector<FLOAT_TYPE> x = create_one_dim_matrix(n, fun_x);
	vector<FLOAT_TYPE> y = multiply_matrices(matrix, x);
	vector<int> p = create_permutation_iden(matrix); // 1-1, 2-2, 3-3... change p, to change order of rows

	cout << std::setprecision(2);
	cout<<"Matrix A\n";
	for(int i=0; i<matrix.size(); i++) {
	 	print_vector(matrix[i]);
	}
	cout<<endl;

	cout<<"Matrix x\n";
	print_vector(x);
	cout<<endl;

	gauss_elimination(matrix, y, p);
	vector<FLOAT_TYPE> res_x = calculate_x(matrix, y, p);

	cout<<"\nRow echelon form\n";

	for(int i=0; i<matrix.size(); i++) {
	 	print_vector(matrix[p[i]]);
	}

	cout<<"\n\nCalculated x\n";
	for(int i=0; i<matrix.size(); i++){
		cout<<std::setprecision(6)<<res_x[p[i]]<<" ";
	}
	cout<<"\n\n";
	cout<<"Euclidan distance between x: " << calculate_euclidan_distance(x, res_x)<<endl<<endl;
}


void thomas(vector<FLOAT_TYPE> a, vector<FLOAT_TYPE> b, vector<FLOAT_TYPE> c, vector<FLOAT_TYPE> &d) {
	int n=a.size();
    /*
    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|
	Result is in d
    */
    n=n-1;
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
	c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i=n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }
}

vector<FLOAT_TYPE> test_thomas(vector<vector<FLOAT_TYPE> > matrix, vector<FLOAT_TYPE>y){
	int n = matrix.size();

	vector<FLOAT_TYPE> a(n-1, 0);
	vector<FLOAT_TYPE> b(n, 0);
	vector<FLOAT_TYPE> c(n-1, 0);

	for(int i=0; i<n; i++){
		b[i]=matrix[i][i];
		if(i+1<n){
			a[i]=matrix[i+1][i];
			c[i]=matrix[i][i+1];
		}
	}

	thomas(a,b,c,y);

	return y;
}

void test_gauss(int n, FLOAT_TYPE (*fun)(int,int)){
	vector<vector<FLOAT_TYPE> > matrix = create_matrix(n, fun);
	vector<FLOAT_TYPE> x = create_one_dim_matrix(n, fun_x);
	vector<FLOAT_TYPE> y = multiply_matrices(matrix, x);
	vector<int> p = create_permutation_iden(matrix); // 1-1, 2-2, 3-3... change p, to change order of rows


	gauss_elimination(matrix, y, p);
	vector<FLOAT_TYPE> res_x1 = calculate_x(matrix, y, p);
	FLOAT_TYPE eucl1 =calculate_euclidan_distance(x, res_x1);


	printf("%d,%f\n",n, eucl1);
}

void gauss_statistics(int n, FLOAT_TYPE (*fun)(int,int)){
	for(int i=10; i<n; i+=10){
		test_gauss(i, fun);
	}
}

void compare_thomas_gauss(int n, FLOAT_TYPE (*fun)(int,int)){
	vector<vector<FLOAT_TYPE> > matrix = create_matrix(n, fun);
	vector<int> p = create_permutation_iden(matrix);
	vector<FLOAT_TYPE> x = create_one_dim_matrix(n, fun_x);
	vector<FLOAT_TYPE> y = multiply_matrices(matrix, x);

	vector<vector<FLOAT_TYPE> > matrix_g = matrix;
	vector<FLOAT_TYPE> y_g = y;

	clock_t begin_time_gauss = clock();
	gauss_elimination(matrix, y_g, p);
	vector<FLOAT_TYPE> res_x = calculate_x(matrix_g, y_g, p);
	clock_t end_time_gauss = clock();


	FLOAT_TYPE eucl_dist_gauss = calculate_euclidan_distance(x, res_x);


	n = matrix.size();

	vector<FLOAT_TYPE> a(n-1, 0);
	vector<FLOAT_TYPE> b(n, 0);
	vector<FLOAT_TYPE> c(n-1, 0);

	for(int i=0; i<n; i++){
		b[i]=matrix[i][i];
		if(i+1<n){
			a[i]=matrix[i+1][i];
			c[i]=matrix[i][i+1];
		}
	}
	clock_t begin_time_thomas = clock();
	thomas(a,b,c,y);
	clock_t end_time_thomas = clock();

	vector<FLOAT_TYPE> res_x_th = y;
	FLOAT_TYPE eucl_dist_thomas = calculate_euclidan_distance(x, res_x_th);

	float gauss_time = (float)(end_time_gauss-begin_time_gauss)/CLOCKS_PER_SEC;
	float thomas_time = (float)(end_time_thomas-begin_time_thomas)/CLOCKS_PER_SEC;

	printf("%d,%f,%f,%f,%f\n",n, eucl_dist_gauss, eucl_dist_thomas, gauss_time, thomas_time);
}

void thomas_gauss_statistics(int n, FLOAT_TYPE (*fun)(int,int)){
	for(int i=5; i<n; i+=5){
		compare_thomas_gauss(i, fun);
	}
}


int main(){
	srand(time(NULL));
	try{
		cout<<"Zadanie 1"<<endl;
		test_one_matrix(10, fun1);

		cout<<"Zadanie 2"<<endl;
		test_one_matrix(20, fun2);

		//thomas_gauss_statistics(500, fun4);
		//gauss_statistics(500, fun2);
		
	}
	catch(char const* error){
		cout<<error<<endl;
	}
	
	return 0;
}