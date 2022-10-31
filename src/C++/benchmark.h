/*
 * Benchmark for different sized gridworlds
*/

#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "HMM.h"
#include "environment.h"
#include "examples.h"
#include <random>
#include <iostream>
using namespace std;

class Benchmark{
public:
    Benchmark(){}

    void run(int grid_type=3, int DFA_states=3, int n = 34, int seq = 275, int maxIt = 25000, double tol = 1e-6, double init_P_est = 1.0, string filename = "", int mode = 0){
        printf("Benchmarking gridworld %d\n", grid_type);

        // VEC2D(double) SP_init;
        // VEC2D(double) SP_init(N_init, VEC(double)(N_init,init_P_est/(N_init-4)));

        int N_init = 1;
        
        if(mode == 1){
            grid_MDP = new Example(grid_type, DFA_states, mode);
            N = grid_MDP->N;
            M = grid_MDP->M;
            VEC(int) params_MDP = {N,M,n,seq};

            N_init = N * DFA_states;
            SP_init = (double **) malloc(N_init*sizeof(double *));
            for(int i=0; i < N_init; i++){
                SP_init[i] = (double *) malloc(N_init*sizeof(double));
                for(int j=0; j < N_init; j++){
                    SP_init[i][j] = init_P_est/(N_init-4);
                }
            }

            // Initialise estimates
            SP_est = (double **) malloc(N * sizeof(double *));
            SE_est = (double **) malloc(N * sizeof(double *));

            for(int i=0; i<N; i++){
                SP_est[i] = (double *) malloc(N * sizeof(double));
                SE_est[i] = (double *) malloc(M * sizeof(double));
                for(int j=0; j<N; j++){
                    SP_est[i][j] = 1.0;
                }
                for(int j=0; j<M; j++){
                    SE_est[i][j] = grid_MDP->E_true[i][j];
                }
            }

            hmm_MDP = new HMM(grid_MDP->P,grid_MDP->E_true,tol,maxIt,params_MDP);
            learn_spatial_MDP(grid_type, DFA_states, n, seq, maxIt, tol, init_P_est, hmm_MDP, SP_init);

            // free allocated objects
            for(int i=0; i<N_init; i++){
                free(SP_init[i]);
            }
            free(SP_init);
            // printf("SP_init freed\n");
            for(int i=0; i<N; i++){
                free(SP_est[i]);
                free(SE_est[i]);
            }
            free(SP_est);
            free(SE_est);
            // printf("SP_est freed\n");
            delete(hmm_MDP);
            // printf("hmm_MDP freed\n");
            delete(grid_MDP);
            // printf("grid_MDP freed\n");
        }
        else{

        // Create example
        grid = new Example(grid_type, DFA_states, 0);
        N = grid->N;
        M = grid->M;

        VEC(int) params = {N,M,n,seq};
        
        // printf("Initialising guesses\n");

        printf("DFA states: %d\n", DFA_states);
        printf("init_P_est: %f\n", init_P_est);
        
        // Initialise estimates
        P_est = (double **) malloc(N * sizeof(double *));
        E_est = (double **) malloc(N * sizeof(double *));

        for(int i=0; i<N; i++){
            P_est[i] = (double *) malloc(N * sizeof(double));
            E_est[i] = (double *) malloc(M * sizeof(double));
            double sum = 0.0;
            for(int j=0; j<N; j++){
                if(grid->P_init[i][j] > 0){
                    P_est[i][j] = grid->P_init[i][j];
                }
                else{
                    P_est[i][j] = init_P_est/(N-4);
                }
                sum += P_est[i][j];
            }
            // Normalise probabilities
            for(int j=0; j<N; j++){
                P_est[i][j] /= sum;
            }         
            for(int j=0; j<M; j++){
                E_est[i][j] = grid->E_true[i][j];
            }
        }

        printf("Initializing HMM\n");
        //Run HMM
        printf("tol: %f\n", tol);
        if(mode == 1){
            hmm = new HMM(grid->P,grid->E_true,tol,maxIt,params,"spatial_MDP_est.txt");
        }
        else{
            hmm = new HMM(grid->P,grid->E_true,tol,maxIt,params,filename);
        }
        
        // printf("Running HMM\n");
        int uniform = 0;
        if(init_P_est == 0.0 && mode == 0){
            uniform = 1;
        }
        double l = hmm->run(P_est,E_est,uniform,0);
        printf("Log likelihood: %f\n", l);
        
        delete(hmm);
        for(int i=0; i<N; i++){
            free(P_est[i]);
            free(E_est[i]);
        }
        free(P_est);
        free(E_est);
        delete(grid);

        }

        return;
    }
    
    void learn_spatial_MDP(int grid_type=3, int DFA_states=3, int n = 34, int seq = 275, int maxIt = 25000, double tol = 1e-6, double init_P_est = 1.0, HMM * hmm_MDP = NULL, double ** SP_init = NULL){
        // Create example
        // N = grid_MDP->N;
        // M = grid_MDP->M;

        VEC(int) params = {N,M,n,seq};
        
        printf("Learning spatial MDP\n");

        // Run HMM learning with uniform prior to learn the spatial MDP
        double l_MDP = hmm_MDP->run(SP_est,SE_est,1,1);
        printf("Log likelihood for spatial MDP: %f\n", l_MDP);

        // Copy learnt spatial MDP but add small probabilities to zero entries
        int N_init = N * DFA_states;
        // VEC2D(double) SP_init(N_init, VEC(double)(N_init,init_P_est/(N_init-4)));
        VEC(double) sum(N,0.0);

        // Copy the estimated spatial MDP to an initial estimate of the product MDP
        for(int d=0; d<DFA_states; d++){
            for(int i=0; i<N; i++){
                for(int j=0; j<N; j++){
                    if(SP_est[i][j] != 0.0){
                        SP_init[(d*N)+i][(d*N)+j] = SP_est[i][j];
                    }
                }
            }
        }

        // printf("Estimated spatial MDP copied\n");

        // Normalise entries
        for(int i=0; i<N_init; i++){
            //Sum each row
            for(int j=0; j<N_init; j++){
                sum[i] += SP_init[i][j];
            }
            //Normalise row
            for(int j=0; j<N_init; j++){
                SP_init[i][j] = SP_init[i][j] / sum[i];
            }
        }

        // printf("Entries normalised\n");

        // write P_est to file
        FILE * f = (FILE *) fopen("spatial_MDP_est.txt","w");
        
        // fopen("spatial_MDP_est.txt", "w");
        for(int i=0; i<N_init; i++){
            for(int j=0; j<N_init; j++){
                fprintf(f, "%f ", SP_init[i][j]);
            }
            if(i<N_init-1){
                fprintf(f, "\n");
            }
        }

        // printf("File written\n");

        // //close file
        fclose(f);
        printf("Spatial MDP estimated\n");
    }
private:
    Example * grid;
    int N;
    int M;
    double ** P_est;
    double ** E_est;
    HMM * hmm;
    // Spatial MDP variables
    double ** SP_init;
    double ** SP_est;
    double ** SE_est;
    HMM * hmm_MDP;
    Example * grid_MDP;
};

#endif // BENCHMARK_H