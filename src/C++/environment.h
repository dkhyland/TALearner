/*
 * Class for Environment
*/

#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "helper.h"

class Environment{
public:

Environment(){}

VEC2D(double) calcProdMC(VEC2D(double) MC, VEC3D(int) D, VEC(int) L){
    // printf("calcProdMC init\n");
    // Calculate dimensions
    int s_states = MC.size();
    int q_states = D.size();
    int n_prods = s_states * q_states;

    // Initialize product MC
    VEC2D(double) P(n_prods, VEC(double)(n_prods,0.0));
    
    // Calculate product MC
    for(int i=0; i<s_states; i++){
        for(int j=0; j<s_states; j++){
            for(int k=1; k<=q_states; k++){
                for (int l=1; l<=q_states; l++){
                    P[(k-1)*s_states + i][(l-1)*s_states + j] = MC[i][j] * D[k-1][l-1][L[j]-1];
                }
            }
        }
    }
    m_P = P;

    // printf("calcProdMC done\n");
    return P;
}

VEC3D(int) DFA(VEC2D(int) list){
    // Calculate dimensions
    int n_transitions = list.size();
    int n_alph = max_val(list,2); // C++11
    int n_q = MAX(max_val(list,0), max_val(list,1));
    // printf("DFA start, %d, %d\n",n_alph,n_q);
    
    //Initialise D to n_q x n_q x n_alph matrix with zeros
    VEC3D(int) D(n_q, VEC2D(int)(n_q, VEC(int)(n_alph, 0)));
    
    //Calculate D
    for(int i=0; i < n_transitions; i++){
        D[list[i][0]-1][list[i][1]-1][list[i][2]-1] = 1;
    };

    m_D = D;
    // printf("DFA done\n");
    return D;
} 

private:
    // DFA
    VEC3D(int) m_D;
    // Product MC
    VEC2D(double) m_P;
};

#endif // ENVIRONMENT_H