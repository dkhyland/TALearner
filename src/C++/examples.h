/*
 * Class for Examples
*/

#ifndef EXAMPLES_H
#define EXAMPLES_H

#include "environment.h"
#include <random>

class Example{
public:
    Example(int grid_type = 3, int m_DFA_states = 3, int mode = 0):type(grid_type){
        // Set up example depending on which grid size is specified
        DFA_states = m_DFA_states;
        env = new Environment();

        switch(DFA_states){
            case 3:
                init_list = init_list_3;
                break;
            case 4:
                init_list = init_list_4;
                break;
            case 5:
                init_list = init_list_5;
                break;
        }   

        if(mode == 0){
            switch(type){
            case 3:
                Example_3x3();
                break;
            case 4:
                Example_4x4();
                break;
            case 5:
                Example_5x5();
                break;
            case 6:
                Example_Office_World();
                init_list = init_list_4;
                break;
            default:
                return;
            }
            // Initialize environment
            VEC2D(double) P_chain(MC_states, VEC(double)(MC_states,0.0));

            if(type < 6){
                //Make the MC
                for (int i=0; i < height; i++){
                    for (int j=0; j < width; j++){
                        int coord = grid_coord(i,j,width);
                        if(i != height - 1){
                            P_chain[coord][grid_coord(i+1,j,width)] = 0.25;
                        }
                        else{
                            P_chain[coord][coord] += 0.25;
                        }
                        if(j != width - 1){
                            P_chain[coord][grid_coord(i,j+1,width)] = 0.25;
                        }
                        else{
                            P_chain[coord][coord] += 0.25;
                        }
                        if(i != 0){
                            P_chain[coord][grid_coord(i-1,j,width)] = 0.25;
                        }
                        else{
                            P_chain[coord][coord] += 0.25;
                        }
                        if(j != 0){
                            P_chain[coord][grid_coord(i,j-1,width)] = 0.25;
                        }
                        else{
                            P_chain[coord][coord] += 0.25;
                        }
                    }
                }
            }else{
                for(int i=0; i < height; i++){
                    for (int j=0; j < width; j++){
                        int coord = grid_coord(i,j,width);
                        // up
                        // printf("up: %d,%d\n",i,j);
                        switch(i){
                            case 5:
                                if (j != 1 && j != 4 && j != 7 && j != 10){
                                    P_chain[coord][coord] += 0.25;
                                }
                                else{
                                    P_chain[coord][grid_coord(i+1,j,width)] = 0.25;
                                }
                                break;
                            case 2:
                                if (j != 1  && j != 10){
                                    P_chain[coord][coord] += 0.25;
                                }
                                else{
                                    P_chain[coord][grid_coord(i+1,j,width)] = 0.25;
                                }
                                break;
                            case 8:
                                P_chain[coord][coord] += 0.25;
                                break;
                            default:
                                P_chain[coord][grid_coord(i+1,j,width)] = 0.25;
                                break;
                        }
                        // down
                        // printf("down: %d,%d\n",i,j);
                        switch(i){
                            case 0:
                                P_chain[coord][coord] += 0.25;
                                break;
                            case 3:
                                if (j != 1  && j != 10){
                                    P_chain[coord][coord] += 0.25;
                                }
                                else{
                                    P_chain[coord][grid_coord(i-1,j,width)] = 0.25;
                                }
                                break;
                            case 6:
                                if (j != 1 && j != 4 && j != 7 && j != 10){
                                    P_chain[coord][coord] += 0.25;
                                }
                                else{
                                    P_chain[coord][grid_coord(i-1,j,width)] = 0.25;
                                }
                                break;
                            default:
                                P_chain[coord][grid_coord(i-1,j,width)] = 0.25;
                                break;
                        }
                        // left
                        // printf("left: %d,%d\n",i,j);
                        switch(j){
                            case 0:
                                P_chain[coord][coord] += 0.25;
                                break;
                            case 3:
                                if(i != 1 && i != 7){
                                    P_chain[coord][coord] += 0.25;
                                }
                                else{
                                    P_chain[coord][grid_coord(i,j-1,width)] = 0.25;
                                }
                                break;
                            case 6:
                                if(i != 1 && i != 7){
                                    P_chain[coord][coord] += 0.25;
                                }
                                else{
                                    P_chain[coord][grid_coord(i,j-1,width)] = 0.25;
                                }
                                break;
                            case 9:
                                if(i != 1 && i != 7){
                                    P_chain[coord][coord] += 0.25;
                                }
                                else{
                                    P_chain[coord][grid_coord(i,j-1,width)] = 0.25;
                                }
                                break;
                            default:
                                P_chain[coord][grid_coord(i,j-1,width)] = 0.25;
                        }
                        // right
                        // printf("right: %d,%d\n",i,j);
                        switch(j){
                            case 2:
                                if(i != 1 && i != 7){
                                    P_chain[coord][coord] += 0.25;
                                }
                                else{
                                    P_chain[coord][grid_coord(i,j+1,width)] = 0.25;
                                }
                                break;
                            case 5:
                                if(i != 1 && i != 7){
                                    P_chain[coord][coord] += 0.25;
                                }
                                else{
                                    P_chain[coord][grid_coord(i,j+1,width)] = 0.25;
                                }
                                break;
                            case 8:
                                if(i != 1 && i != 7){
                                    P_chain[coord][coord] += 0.25;
                                }
                                else{
                                    P_chain[coord][grid_coord(i,j+1,width)] = 0.25;
                                }
                                break;
                            case 11:
                                P_chain[coord][coord] += 0.25;
                                break;
                            default:
                                P_chain[coord][grid_coord(i,j+1,width)] = 0.25;
                        }
                    }
                }

            }

            VEC(int) L(MC_states,1);

            for(int i=0; i< (int)labels.size(); i++){
                L[grid_coord(labels[i][0],labels[i][1],width)] = labels[i][2];
            }

            printf("Calculating Product MC\n");
            // Calculate initial estimate of product MC based on own policy and guess of number of states in the DFA
            P_init = env->calcProdMC(P_chain, env->DFA(init_list), L);
            // Calculate true product MC
            P = env->calcProdMC(P_chain, env->DFA(list), L);
            n_accept = 1;
            VEC2D(double) I = Identity(MC_states);
            
            VEC2D(double) kron_arg(DFA_states, VEC(double)(2,0.0));
            for(int i=0; i<DFA_states-n_accept; i++){
                kron_arg[i][0] = 1;
            }
            kron_arg[DFA_states-1][1] = 1;

            E_true = KroneckerProduct(kron_arg, I, DFA_states, 2, MC_states, MC_states);
        }
        else{
            printf("Mode 1\n");
            switch(type){
                case 3:
                    Example_3x3_2();
                    break;
                case 4:
                    Example_4x4_2();
                    break;
                case 5:
                    Example_5x5_2();
                    break;
            }
            // Initialize environment
            VEC2D(double) P_chain(MC_states, VEC(double)(MC_states,0.0));
            printf("Creating gridworld %d\n", grid_type);
            //Make the MC
            for (int i=0; i < height; i++){
                for (int j=0; j < width; j++){
                    int coord = grid_coord(i,j,width);
                    if(i != height - 1){
                        P_chain[coord][grid_coord(i+1,j,width)] = 0.25;
                    }
                    else{
                        P_chain[coord][coord] += 0.25;
                    }
                    if(j != width - 1){
                        P_chain[coord][grid_coord(i,j+1,width)] = 0.25;
                    }
                    else{
                        P_chain[coord][coord] += 0.25;
                    }
                    if(i != 0){
                        P_chain[coord][grid_coord(i-1,j,width)] = 0.25;
                    }
                    else{
                        P_chain[coord][coord] += 0.25;
                    }
                    if(j != 0){
                        P_chain[coord][grid_coord(i,j-1,width)] = 0.25;
                    }
                    else{
                        P_chain[coord][coord] += 0.25;
                    }
                }
            }

            VEC(int) L(MC_states,1);

            for(int i=0; i< (int)labels.size(); i++){
                L[grid_coord(labels[i][0],labels[i][1],width)] = labels[i][2];
            }

            P = P_chain;
            n_accept = 1;
            VEC2D(double) I = Identity(MC_states);
            E_true = I;
        }

        DFA_state_bound_guess = DFA_states;
        n_guess_accept = 1;
        // n = MC_states * DFA_state_bound_guess;
        // delete(env);
    }

    /* Conventions for gridworlds
       - (0,0) coordinate is the bottom left corner.
       - Labels are integers starting from 1
    */

    /* 3x3 grid world example
     * 1 = nothing
     * 2 = coffee
     * 3 = tv
     * 4 = couch
     * 5 = lemon
     * 6 = stairs
    */
    void Example_3x3(){
        // Initialize environment
        width = 3;
        height = 3;
        MC_states = width * height;
        N = MC_states * DFA_states;
        M = MC_states * 2;
        num_labels = 6;
        // (start state, end state, label)
        switch(DFA_states){
            case 3:
                list = list_3;
                break;
            case 4:
                list = list_4;
                break;
            case 5:
                list = list_5;
                break;
            default:
                list = list_3;
                break;
        }

        labels = {
            {0,width - 1,2}, // coffee
            {1,0,4}, // couch
            {1,width - 1,3}, // tv
            {2,1,6}, // stairs
            {2,width - 1,5} // lemon
        };
    
        return;
    }

    // 4x4 grid world example
    void Example_4x4(){
        // Initialize environment
        width = 4;
        height = 4;
        MC_states = width * height;
        N = MC_states * DFA_states;
        M = MC_states * 2;
        num_labels = 6;
        // (start state, end state, label)
        switch(DFA_states){
            case 3:
                list = list_3;
                break;
            case 4:
                list = list_4;
                break;
            case 5:
                list = list_5;
                break;
            default:
                list = list_3;
                break;
        }

        labels = {
            {0,width-1,2}, // coffee
            {1,0,4}, // couch
            {1,width-1,3}, // tv
            {height-1,1,6}, // stairs
            {height-1,width-1,5} // lemon
        };

        return;
    }

    // 5x5 grid world example
    void Example_5x5(){
        // Initialize environment
        width = 5;
        height = 5;
        MC_states = width * height;
        N = MC_states * DFA_states;
        M = MC_states * 2;
        num_labels = 6;
        // (start state, end state, label)
        switch(DFA_states){
            case 3:
                list = list_3;
                break;
            case 4:
                list = list_4;
                break;
            case 5:
                list = list_5;
                break;
            default:
                list = list_3;
                break;
        }

        labels = {
            {0,4,2}, // coffee
            {2,0,4}, // couch
            {2,4,3}, // tv
            {3,3,5}, // lemon
            {3,4,5}, // lemon
            {4,1,4}, // couch
            {4,2,6}, // stairs
            {4,3,5}, // lemon
            {4,4,5} // lemon
        };
        return;
    }

    void Example_3x3_2(){
        width = 3;
        height = 3;
        N = width * height;
        M = N;
        MC_states = N;
        num_labels = 6;
        init_list = list_3;
        P_init = Identity(N);
    }

    void Example_4x4_2(){
        width = 4;
        height = 4;
        N = width * height;
        M = N;
        MC_states = N;
        num_labels = 6;
        init_list = list_4;
        P_init = Identity(N);
    }

    void Example_5x5_2(){
        width = 5;
        height = 5;
        N = width * height;
        M = N;
        MC_states = N;
        num_labels = 6;
        init_list = list_5;
        P_init = Identity(N);
    }

    // Office World example
    // 1 = nothing
    // 2 = decoration
    // 3 = coffee
    // 4 = office
    // 5 = mail
    void Example_Office_World(){
        // Initialize environment
        width = 12;
        height = 9;
        MC_states = width * height;
        DFA_states = 4;
        N = MC_states * DFA_states;
        M = MC_states * 2;
        num_labels = 5;

        list = {
            {1,1,1},
            {1,3,2},
            {1,2,3},
            {1,1,4},
            {1,1,5},
            {2,2,1},
            {2,3,2},
            {2,2,3},
            {2,4,4},
            {2,2,5},
            {3,3,1},
            {3,3,2},
            {3,3,3},
            {3,3,4},
            {3,3,5},
            {4,4,1},
            {4,4,2},
            {4,4,3},
            {4,4,4},
            {4,4,5}
        };

        labels = {
            {1,4,2}, // decoration
            {1,7,2},
            {2,8,3}, // coffee
            {4,1,2},
            {4,4,4}, // office
            {4,7,5}, // mail
            {4,10,2},
            {6,3,3}, 
            {7,4,2},
            {7,7,2}
        };
        return;
    }
    
    // Initialise the MC with a random exploration policy. grid_type = 0: plain gridworld, grid_type = 1: office world
    void init_MC(int grid_type = 0){
        printf("Initializing MC\n");
        if(grid_type == 0){     
            //Make the MC
            for (int i=0; i < height; i++){
                for (int j=0; j < width; j++){
                    int coord = grid_coord(i,j,width);
                    if(i != height - 1){
                        P_chain[coord][grid_coord(i+1,j,width)] = 0.25;
                    }
                    else{
                        P_chain[coord][coord] += 0.25;
                    }
                    if(j != width - 1){
                        P_chain[coord][grid_coord(i,j+1,width)] = 0.25;
                    }
                    else{
                        P_chain[coord][coord] += 0.25;
                    }
                    if(i != 0){
                        P_chain[coord][grid_coord(i-1,j,width)] = 0.25;
                    }
                    else{
                        P_chain[coord][coord] += 0.25;
                    }
                    if(j != 0){
                        P_chain[coord][grid_coord(i,j-1,width)] = 0.25;
                    }
                    else{
                        P_chain[coord][coord] += 0.25;
                    }
                }
            }
        }else{
            //Make the MC
            for(int i=0; i < height; i++){
                for (int j=0; j < width; j++){
                    int coord = grid_coord(i,j,width);
                    // up
                    printf("up: %d,%d\n",i,j);
                    switch(i){
                        case 5:
                            if (j != 1 && j != 4 && j != 7 && j != 10){
                                P_chain[coord][coord] += 0.25;
                            }
                            else{
                                P_chain[coord][grid_coord(i+1,j,width)] = 0.25;
                            }
                            break;
                        case 2:
                            if (j != 1  && j != 10){
                                P_chain[coord][coord] += 0.25;
                            }
                            else{
                                P_chain[coord][grid_coord(i+1,j,width)] = 0.25;
                            }
                            break;
                        case 8:
                            P_chain[coord][coord] += 0.25;
                            break;
                        default:
                            P_chain[coord][grid_coord(i+1,j,width)] = 0.25;
                            break;
                    }
                    // down
                    printf("down: %d,%d\n",i,j);
                    switch(i){
                        case 0:
                            P_chain[coord][coord] += 0.25;
                            break;
                        case 3:
                            if (j != 1  && j != 10){
                                P_chain[coord][coord] += 0.25;
                            }
                            else{
                                P_chain[coord][grid_coord(i-1,j,width)] = 0.25;
                            }
                            break;
                        case 6:
                            if (j != 1 && j != 4 && j != 7 && j != 10){
                                P_chain[coord][coord] += 0.25;
                            }
                            else{
                                P_chain[coord][grid_coord(i-1,j,width)] = 0.25;
                            }
                            break;
                        default:
                            P_chain[coord][grid_coord(i-1,j,width)] = 0.25;
                            break;
                    }
                    // left
                    printf("left: %d,%d\n",i,j);
                    switch(j){
                        case 0:
                            P_chain[coord][coord] += 0.25;
                            break;
                        case 3:
                            if(i != 1 && i != 7){
                                P_chain[coord][coord] += 0.25;
                            }
                            else{
                                P_chain[coord][grid_coord(i,j-1,width)] = 0.25;
                            }
                            break;
                        case 6:
                            if(i != 1 && i != 7){
                                P_chain[coord][coord] += 0.25;
                            }
                            else{
                                P_chain[coord][grid_coord(i,j-1,width)] = 0.25;
                            }
                            break;
                        case 9:
                            if(i != 1 && i != 7){
                                P_chain[coord][coord] += 0.25;
                            }
                            else{
                                P_chain[coord][grid_coord(i,j-1,width)] = 0.25;
                            }
                            break;
                        default:
                            P_chain[coord][grid_coord(i,j-1,width)] = 0.25;
                    }
                    // right
                    printf("right: %d,%d\n",i,j);
                    switch(j){
                        case 2:
                            if(i != 1 && i != 7){
                                P_chain[coord][coord] += 0.25;
                            }
                            else{
                                P_chain[coord][grid_coord(i,j+1,width)] = 0.25;
                            }
                            break;
                        case 5:
                            if(i != 1 && i != 7){
                                P_chain[coord][coord] += 0.25;
                            }
                            else{
                                P_chain[coord][grid_coord(i,j+1,width)] = 0.25;
                            }
                            break;
                        case 8:
                            if(i != 1 && i != 7){
                                P_chain[coord][coord] += 0.25;
                            }
                            else{
                                P_chain[coord][grid_coord(i,j+1,width)] = 0.25;
                            }
                            break;
                        case 11:
                            P_chain[coord][coord] += 0.25;
                            break;
                        default:
                            P_chain[coord][grid_coord(i,j+1,width)] = 0.25;
                    }
                }
            }
        }
    }

    // Get coffee, then go up stairs
    VEC2D(int) list_3 = {
                {1,1,1},
                {1,2,2},
                {1,1,3},
                {1,1,4},
                {1,1,5},
                {1,1,6},
                {2,2,1},
                {2,2,2},
                {2,2,3},
                {2,2,4},
                {2,2,5},
                {2,3,6},
                {3,3,1},
                {3,3,2},
                {3,3,3},
                {3,3,4},
                {3,3,5},
                {3,3,6},
                };
    // Get coffee, go to the couch, then go up the stairs
    VEC2D(int) list_4 = {
                {1,1,1},
                {1,2,2},
                {1,1,3},
                {1,1,4},
                {1,1,5},
                {1,1,6},
                {2,2,1},
                {2,2,2},
                {2,2,3},
                {2,3,4},
                {2,2,5},
                {2,2,6},
                {3,3,1},
                {3,3,2},
                {3,3,3},
                {3,3,4},
                {3,3,5},
                {3,4,6},
                {4,4,1},
                {4,4,2},
                {4,4,3},
                {4,4,4},
                {4,4,5},
                {4,4,6},
                };
    // Get coffee, go to the couch, go to the tv, and then go up the stairs.
    VEC2D(int) list_5 = {
                {1,1,1},
                {1,2,2},
                {1,1,3},
                {1,1,4},
                {1,1,5},
                {1,1,6},
                {2,2,1},
                {2,2,2},
                {2,2,3},
                {2,3,4},
                {2,2,5},
                {2,2,6},
                {3,3,1},
                {3,3,2},
                {3,4,3},
                {3,3,4},
                {3,3,5},
                {3,3,6},
                {4,4,1},
                {4,4,2},
                {4,4,3},
                {4,4,4},
                {4,4,5},
                {4,5,6},
                {5,5,1},
                {5,5,2},
                {5,5,3},
                {5,5,4},
                {5,5,5},
                {5,5,6},
                };
    VEC2D(int) init_list_3 =  { 
                {1,1,1},
                {1,1,2},
                {1,1,3},
                {1,1,4},
                {1,1,5},
                {1,1,6},
                {2,2,1},
                {2,2,2},
                {2,2,3},
                {2,2,4},
                {2,2,5},
                {2,2,6},
                {3,3,1},
                {3,3,2},
                {3,3,3},
                {3,3,4},
                {3,3,5},
                {3,3,6}
            };
    VEC2D(int) init_list_4 =  { 
                {1,1,1},
                {1,1,2},
                {1,1,3},
                {1,1,4},
                {1,1,5},
                {1,1,6},
                {2,2,1},
                {2,2,2},
                {2,2,3},
                {2,2,4},
                {2,2,5},
                {2,2,6},
                {3,3,1},
                {3,3,2},
                {3,3,3},
                {3,3,4},
                {3,3,5},
                {3,3,6},
                {4,4,1},
                {4,4,2},
                {4,4,3},
                {4,4,4},
                {4,4,5},
                {4,4,6},
            };
    VEC2D(int) init_list_5 =  { 
                {1,1,1},
                {1,1,2},
                {1,1,3},
                {1,1,4},
                {1,1,5},
                {1,1,6},
                {2,2,1},
                {2,2,2},
                {2,2,3},
                {2,2,4},
                {2,2,5},
                {2,2,6},
                {3,3,1},
                {3,3,2},
                {3,3,3},
                {3,3,4},
                {3,3,5},
                {3,3,6},
                {4,4,1},
                {4,4,2},
                {4,4,3},
                {4,4,4},
                {4,4,5},
                {4,4,6},
                {5,5,1},
                {5,5,2},
                {5,5,3},
                {5,5,4},
                {5,5,5},
                {5,5,6},
            };
    int type;
    int example_no;
    int width;
    int height;
    int N;
    int M;
    int num_labels;
    // int n;
    double tolerance;
	int maxIt;
    int DFA_states;
    int DFA_state_bound_guess;
    int MC_states;
    int n_accept;
    int n_guess_accept;
    VEC2D(int) labels;
    VEC2D(double) P;
    VEC2D(double) P_init;
    VEC2D(double) P_chain;
    VEC2D(double) E_true;
    VEC2D(int) list;
    VEC2D(int) init_list;
    Environment * env;
};

#endif // EXAMPLES_H