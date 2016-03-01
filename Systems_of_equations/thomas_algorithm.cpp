#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
using namespace std;
#define TYPE double


class Matrix {
private:
    vector<vector<TYPE> > elements;


    vector<TYPE> calculatePermutations() {

    }

public:
    Matrix(vector<vector<TYPE> > a) {
        this->elements = a;
    }

    Matrix(unsigned size) {
        elements.resize(size, vector<TYPE>(size));

        for(unsigned i=0; i<size; i++) {
            for(unsigned j=0; j<size; j++) {
                elements[i][j] = 0;
            }
        }
    }

    Matrix(unsigned size, TYPE (*fun)(unsigned, unsigned)) {
        elements.resize(size, vector<TYPE>(size));

        for(unsigned i=0; i<size; i++) {
            for(unsigned j=0; j<size; j++) {
                elements[i][j] = fun(i, j);
            }
        }
    }

    unsigned getSize() {
        return elements.size();
    }

    void print() {
        unsigned size = this->elements.size();
        for(unsigned i=0; i<size; i++) {
            for(unsigned j=0; j<size; j++) {
                cout << this->elements[i][j] << " ";
            }
            cout << "\n";
        }
    }

    void gaussianElimination(vector<TYPE> &b) {
        unsigned size = getSize();
        TYPE multiplier;

        for(unsigned k=0; k<size-1; k++) {
            for(unsigned i=k+1; i<size; i++) {
                multiplier = elements[i][k] / elements[k][k];
                elements[i][k] = 0;
                for(unsigned j=k+1; j<size; j++) {
                    elements[i][j] -= multiplier*elements[k][j];
                }
                b[i] -= multiplier*b[k];
            }
        }
    }

    vector<TYPE> operator[](unsigned i) {
        return elements[i];
    }

    Matrix operator*(Matrix &a) {
        unsigned size = getSize();

        if(size != a.getSize()) {
            throw "Matrices cannot be multiplied";
        }

        vector<vector<TYPE> > result(size, vector<TYPE>(size, 0));

        for(unsigned i=0; i<size; i++) {
            for(unsigned j=0; j<size; j++) {
                for(unsigned k=0; k<size; k++) {
                    result[i][j] += elements[i][k]*a[k][j];
                }
            }
        }

        return Matrix(result);      
    }

    vector<TYPE> operator*(vector<TYPE> &a) {
        unsigned size = getSize();

        if(size != a.size()) {
            throw "Matrice and vector cannot be multiplied";
        }

        vector<TYPE> result(size, 0);

        for(unsigned i=0; i<size; i++) {
            for(unsigned j=0; j<size; j++) {
                result[i] += elements[i][j]*a[j];
            }
        }

        return result;
    }

};

vector<TYPE> create_vector(unsigned size, TYPE (*fun)(unsigned i)) {
    vector<TYPE> result(size);

    for(unsigned i=0; i<size; i++) {
        result[i] = fun(i);
    }

    return result;
}

void print_vector(vector<TYPE> &a) {
    for(unsigned i=0; i<a.size(); i++) {
        cout << a[i] << " ";
    }
    cout << "\n";
}

TYPE function_a(unsigned i, unsigned j) {
    int m = 7, k = 2;
    i = i + 1;
    j = j + 1;

    if(i == j) {
        return -m*((TYPE)i) - k;
    } else if(i == j-1) {
        return (TYPE)i;
    } else if(j == i-1 && i > 1) {
        return m/(TYPE)i;
    } else return 0;
}

TYPE random(unsigned i) {
    if(rand() % 2 == 0) {
        return 1;
    }
    return -1;
}

vector<TYPE> solve_linear_equation(Matrix &matrix, vector<TYPE> b) {
    vector<TYPE> result(matrix.getSize(), 0);
    for(int i=matrix.getSize()-1; i>=0; i--){
        TYPE sum=0;
        for(int j=i+1; j<matrix.getSize(); j++){
            sum+=matrix[i][j]*result[j];
        }
        result[i]=(b[i]-sum)/matrix[i][i];
    }
    return result;
}


TYPE calculate_vectors_norm(vector<TYPE> &a, vector<TYPE> &b) {
    TYPE norm = 0;
    for(unsigned i=0; i<a.size(); i++) {
        norm += (a[i]-b[i])*(a[i]-b[i]);
    }
    return sqrt(norm);
}

vector<TYPE> tridiagonal_algorithm(vector<TYPE> &A, vector<TYPE> &B, vector<TYPE> &C, vector<TYPE> &b) {
    unsigned size = A.size();
    vector<TYPE> c_vector(size), d_vector(size), result(size);

    c_vector[0] = C[0]/B[0];
    for (int i = 1; i < size-1; ++i)
    {
        c_vector[i] = C[i]/(B[i] - A[i]*c_vector[i-1]);
    }

    d_vector[0] = b[0]/B[0];
    for (int i = 1; i < size; ++i)
    {
        d_vector[i] = (b[i]-A[i]*d_vector[i-1])/(B[i] - A[i]*c_vector[i-1]);
    }

    result[size-1] = d_vector[size-1];
    for (int i = size-2; i >= 0; --i)
    {
        result[i] = d_vector[i] - c_vector[i]*result[i+1];
    }

    return result;
}

void run(unsigned N) {
    Matrix Q(N, function_a);
    vector<TYPE> X = create_vector(N, random);
    vector<TYPE> b = Q*X;
    vector<TYPE> A(N);
    vector<TYPE> B(N);
    vector<TYPE> C(N);

    for(int i=0; i<N; i++) {
        B[i] = function_a(i, i);
        if(i-1 >= 0) {
            A[i] = function_a(i, i-1);            
        }
        if(i+1 < N) {
            C[i] = function_a(i, i+1);
        }
    }

    auto t1 = chrono::high_resolution_clock::now();

    vector<TYPE> computed_X = tridiagonal_algorithm(A, B, C, b);

    auto t2 = chrono::high_resolution_clock::now();

    // fractional duration: no duration_cast needed
    chrono::duration<double, milli> fp_ms = t2 - t1;
 
    //cout << N  << ", " << fp_ms.count() << "\n";

    cout << N << "," << calculate_vectors_norm(X, computed_X) << "\n";
}

int main()
{
    srand(time(NULL));

    unsigned N = 100;
    for(int i=2; i<=N; i++) {
        run(i);
    }

    return 0;
}
