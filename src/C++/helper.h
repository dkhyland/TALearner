/*
 * Helper functions for LTA_cpp
*/

#ifndef HELPER_H
#define HELPER_H

#include <vector>

using namespace std;

// Macros
#ifndef VEC
#define	VEC(type)	vector<type>
#endif

#ifndef VEC2D
#define	VEC2D(type)	vector<vector<type>>
#endif

#ifndef VEC3D
#define	VEC3D(type)	vector<vector<vector<type>>>
#endif

#ifndef MAX
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

// Compute the Kronecker Product of two matrices 
static VEC2D(double) KroneckerProduct(VEC2D(double) A, VEC2D(double) B, int r_a, int c_a, int r_b, int c_b){ 
    // printf("Kronecker start\n");
    VEC2D(double) K(r_a * r_b, VEC(double)(c_a * c_b,0.0));
    // i loops till r_a 
    for (int i = 0; i < r_a; i++) {
        // j loops till c_a 
        for (int j = 0; j < c_a; j++) { 
            // k loops till r_b 
            for (int k = 0; k < r_b; k++) { 
                // l loops till c_b 
                for (int l = 0; l < c_b; l++) {
                    // Each element of matrix A is 
                    // multiplied by whole Matrix B 
                    // resp and stored as Matrix K
                    K[(i*(r_b)) + k][(j*(c_b)) + l] = A[i][j] * B[k][l];
                } 
            } 
        } 
    } 
    // printf("Kronecker done\n");
    return K;
} 

// Find largest element in 2D matrix
static int max_val(VEC2D(int) v, int idx){
    int max = v[0][idx];
    int rows = v.size();
    // int cols = v[0].size();
    for(int i = 0; i < rows; i++){
        max = MAX(v[i][idx],max);
        
    }   
    return max;
}

// Return grid coordinate as a single index
static int grid_coord(int s1, int s2, int row_length){
    return s2 + s1 * row_length;
}

// Return n x n identity matrix
static VEC2D(double) Identity(int n){
    VEC2D(double) I(n, VEC(double)(n,0.0));
    for (int i = 0; i < n; i++){
        I[i][i] = 1.0;
    }
    return I;
}


#endif // HELPER_H