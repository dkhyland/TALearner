/* 
 * Class for HMM
 *
 * Based on code by Aurelien Garivier, CNRS & Telecom Paristech
 * (small bug corrected by Julius Su, Caltech)
 *
 * Baum-Welch algorithm for discrete Hidden Markov models
 * see http://www.telecom-paristech.fr/~garivier/code/index.html
 * 
*/

#ifndef HMM_H
#define HMM_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <memory.h>
#include <string>
#include <sstream>
#include <random>
#include <fstream>
#include "environment.h"

using namespace std;

class HMM{
public:
	/*
	Constructor for HMM class
	in:   P_est = pointer to beginning of N x N array of initial transition matrix estimate
		  E_est = pointer to beginning of N x M array of initial emission matrix estimate
		  P = pointer to beginning of N x N array of true transition matrix
		  E = pointer to beginning of N x M array of true emission matrix
		  tol = tolerance for the stopping criterion
		  maxIt = maximal number of iterations
		  params = pointer to a vector of length 4 containing (N,M,n,seq)
	*/
	HMM(VEC2D(double) m_P, VEC2D(double) m_E, double m_tol, int m_maxIt, VEC(int) m_params, string m_filename = ""):\
	P(m_P),E(m_E),tol(m_tol),maxIt(m_maxIt),filename(m_filename){
		// P_est = P_e;
		// E_est = E_e;
		ll = 0;
		N = m_params[0];
		M = m_params[1];
		n = m_params[2];
		seq = m_params[3];
		printf("N: %d\n", m_params[0]);
        printf("M: %d\n", m_params[1]);
        printf("n: %d\n", m_params[2]);
		printf("sequences: %d\n", m_params[3]);
		printf("maxIt: %d\n", maxIt);
		printf("tol: %f\n", tol);
	}

	/*
	sample from discrete distribution
	in:   p = vector of probabilities,assumed to sum to 1
	*/
	int randm(VEC(double) p){
		int res=0;
		double q=p[0];
		double u = (rand()+0.0)/RAND_MAX;;
		while(u>q){
			res++;
			q+=p[res];
		}
		return(res);
	}

	/*
	Compute maximum absolute row sum of a matrix
	in:   P1 = transition matrix of size N x N
		  P2 = transition matrix of size N x N
	out:  max = maximum absolute row sum of P1 - P2
	*/
	double norm_P(double ** P1, double ** P2){
		double max = 0;
		double sum = 0;
		int i,j;
		for(i=0;i<N;i++){
			sum=0;
			for(j=0;j<N;j++){
				sum += fabs(P1[i][j] - P2[i][j]);
			}
			max = MAX(sum,max);
		}
		return max;
	}

	double diff_P(double ** P1, VEC2D(double) P2){
		double sum = 0;
		int i,j;
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				sum += fabs(P1[i][j] - P2[i][j]);
			}
		}
		return sum;
	}

	/*
	Compute maximum absolute row sum of a matrix
	in:   E1 = emission matrix of size N x M
		  E2 = emission matrix of size M x N (transposed)
	out:  max = maximum absolute row sum of E1 - E2
	*/
	double norm_E(double ** E1, double ** E2){
		double max = 0;
		double sum = 0;
		int i, j;
		for(i=0;i<N;i++){
			sum=0;
			for(j=0;j<M;j++){
				sum += fabs(E1[i][j] - E2[i][j]);
			}
			max = MAX(sum,max);
		}
		return max;
	}

	/*
	print out transition and emission matrices
	in: P = transition matrix of size N x N
		E = emission matrix of size N x M
	*/
	void print_matrices(VEC2D(double) P, VEC2D(double) E){
		int i,j;
		// print out each element of P and E
		printf("P = \n");
		for (i = 0; i < N; i++){
			for (j = 0; j < N; j++){
				printf("%f ", P[i][j]);
			}
			printf("\n");
		}

		printf("E = \n");
		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
				printf("%f ", E[i][j]);
			}
			printf("\n");
		}
	}

	void print_matrices_ptr(double ** P, double ** E){
		int i,j;
		// print out each element of P and E
		printf("P = \n");
		for (i = 0; i < N; i++){
			for (j = 0; j < N; j++){
				printf("%f ", P[i][j]);
			}
			printf("\n");
		}

		printf("E = \n");
		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
				printf("%f ", E[i][j]);
			}
			printf("\n");
		}
	}

	/*
	sample a trajectory from a hidden markov chain
	in:   nu = initial distribution as vector of size N
		  P = transition matrix of size N x N
		  E = emission matrix of size N x M
	out:  (x,y) = sample trajectory of size n of a HMM defined by (nu, P, E):
		  x = sample trajectory of size n of a Markov Chain with initial distribution nu and transition matrix P
		  y = observations such that the conditional distribution of y[k]
		  given x[k] is E(x[k], :)
	*/
	void HMMsample(VEC(double) nu, VEC2D(double) P, VEC2D(double) E, int * x, int * y){
		int k;
		// Sample first transition
		// (*x)[0] = randm(nu);
		// (*y)[0] = randm(E[(*x)[0]]);
		// (*x)[0] = randm(P[0]);
		// (*y)[0] = randm(E[(*x)[0]]);
		x[0] = randm(P[0]);
		y[0] = randm(E[x[0]]);

		// Sample rest of transitions
		for(k=1; k<n; k++){
			x[k] = randm(P[x[k-1]]);
			y[k] = randm(E[x[k]]);
		}
	} 

	/*
	HMM filtering of an observation sequence, given hmm parameters
	in:   y = vector of observations, assumed to be in range(g.shape[1])
		  nu = initial distribution as vector of size N
		  P_e = transition matrix of size N x N
		  E_e = emission matrix with size N x M
	out:  fs = filter: P(x[t]=x | y[0:t]=y[0:t]) for 0<=x<k and 0<=t<n (size n x N)
		  c(t) = conditional likelihood: P(Y[t] = y[t]| Y[0:t-1]=y[0:t-1]) (size n)
		  ll = log sum of logs of c(t)
	*/
	double HMMfilter(int * y, VEC(double) nu, double ** P_e, double ** E_e, double ** fs,  double * c){
		int i,j,t;
		double ll = 0;
		VEC(double) alpha(N,0.0);
		
		// set fs to 0
		for(i=0;i<n+1;i++){
			for(j=0;j<N;j++){
				fs[i][j] = 0.0;
			}
		}
		
		// Assume we are starting in state 0
		c[0]= 1.0;
		fs[0][0] = 1.0;

		for(t=1; t<n+1; t++){
			c[t]=0;
			for(j=0; j<N; j++){
				alpha[j]=0;
				for(i=0; i<N; i++) {
					alpha[j] += fs[t-1][i] * P_e[i][j];
				}
				fs[t][j] = alpha[j] * E_e[j][y[t]];
				c[t] += fs[t][j];
			}
			 
			for(j=0; j<N; j++){
				fs[t][j] = fs[t][j]/c[t];
			}
			ll += log(c[t]);
		}
		return ll;
	}


	double HMMfilter_opt(int * y, VEC(double) nu, double ** P_e, double ** E_e, double ** fs,  double * c, int * nzc, int ** nzc_row, int opt){
		int i,j,k,t;
		double ll = 0;
		VEC(double) alpha(N,0.0);
		// set fs to 0
		for(t=0;t<n+1;t++){
			for(j=0;j<N;j++){
				fs[t][j] = 0.0;
			}
		}
		
		// Assume we are starting in state 0
		c[0]= 1.0;
		fs[0][0] = 1.0;

		// Run optimised filtering procedure if opt flag is set to 1 or 2 (default is 0)
		if(opt > 0){
			for(t=1; t<n+1; t++){
				c[t]=0;
				for(j=0; j<N; j++){
					alpha[j]=0;
					for(i=0; i<nzc[j]; i++) {
						k = nzc_row[j][i];
						alpha[j] += fs[t-1][k] * P_e[k][j];
					}
					fs[t][j] = alpha[j] * E_e[j][y[t]];
					c[t] += fs[t][j];
				}
				
				for(j=0; j<N; j++){
					fs[t][j] = fs[t][j]/c[t];
				}
				ll += log(c[t]);
			}
		}
		else{
			for(t=1; t<n+1; t++){
				c[t]=0;
				for(j=0; j<N; j++){
					alpha[j]=0;
					for(i=0; i<N; i++) {
						alpha[j] += fs[t-1][i] * P_e[i][j];
					}
					fs[t][j] = alpha[j] * E_e[j][y[t]];
					c[t] += fs[t][j];
				}
				
				for(j=0; j<N; j++){
					fs[t][j] = fs[t][j]/c[t];
				}
				ll += log(c[t]);
			}
		}
		return ll;
	}

	/*
	HMM smoothing of an observation sequence, given hmm parameters
	in:   y = vector of observations, assumed to be in range(P.shape[0])
		  P_e = estimate of transition matrix of size N x N
		  E_e = estimate emission matrix of size N x M
		  c = conditional likelihoods, computed by HMMfilter
	out:  bs = smoothing factors: P(y[t+1:n]=y[t+1:n] | X[t]=x) / P(Y[t+1:n]=y[t+1:n] | Y[0:t]=y[1:t]) for 0<=x<N and 1<=t<n
		  permits to compute the posterior distribution of the hidden states 
		  P(X[t]=x | Y[0:n]=y[0:n])  as post = fs .* bs
	*/
	void HMMsmoother(int * y, double ** P_e, double ** E_e, double * c, double ** bs){
		int i,j,t;
		double z;
		for(t=0;t<n+1;t++){
			for(j=0;j<N;j++){
				bs[t][j] = 1;
			}
		}

		for(t=n-1; t>=0; t--){
			for(i=0; i<N; i++){
				z = 0;
				for(j=0; j<N; j++){
					z += P_e[i][j] * E_e[j][y[t+1]] * bs[t+1][j];
				}
				bs[t][i] = z / c[t+1];
			}
		}
	}

	void HMMsmoother_opt(int * y, double ** P_e, double ** E_e, double * c, double ** bs, int * nzr, int ** nzr_col, int opt){
		int i,j,k,t;
		double z;
		for(t=0;t<n+1;t++){
			for(j=0;j<N;j++){
				bs[t][j] = 1;
			}
		}

		if(opt > 0){
			for(t=n-1; t>=0; t--){
				for(i=0; i<N; i++){
					z = 0;
					for(j=0; j<nzr[i]; j++){
						k = nzr_col[i][j];
						z += P_e[i][k] * E_e[k][y[t+1]] * bs[t+1][k];
					}
					bs[t][i] = z / c[t+1];
				}
			}
		}
		else{
			for(t=n-1; t>=0; t--){
				for(i=0; i<N; i++){
					z = 0;
					for(j=0; j<N; j++){
						z += P_e[i][j] * E_e[j][y[t+1]] * bs[t+1][j];
					}
					bs[t][i] = z / c[t+1];
				}
			}
		}
		
	}

	/*
	utility functions: sample random transition and emission kernels
	*/
	void randomTransitionKernel(double ** K){
		int i,j;
		double s;
		for(i=0; i<N; i++){
			s=0;
			for(j=0; j<N; j++)
				s+=(K[i][j] = (rand()+0.0)/RAND_MAX);
			for(j=0; j<N; j++)
				K[i][j] /= s;
		}
	}

	void randomEmissionKernel(double ** K){
		int i,j;
		double s;
		for(i=0; i<N; i++){
			s=0;
			for(j=0; j<M; j++)
				s+=(K[i][j] = (rand()+0.0)/RAND_MAX);
			for(j=0; j<M; j++)
				K[i][j] /= s;
		}
	}

	void uniformTransitionKernel(double ** K){
		int i,j;
		for(i=0; i<N; i++){
			for(j=0; j<N; j++){
				K[i][j] = 1.0/N;
			}
		}
	}

	void copyEmissionKernel(double ** K, double ** K2){
		int i,j;
		for(i=0; i<N; i++){
			for(j=0; j<M; j++){
				K[i][j] = K2[i][j];
			}
		}
	}

	/*
	compute maximum likehood estimate using Expectation-Maximization
	iterations
	in:   y = vector of observations 
		  nu = initial distribution of the hidden chain
		  tol = tolerance for the stopping criterion
		  maxIt = maximal number of iterations
	out:  P = estimate of the transition matrix of the hidden markov process
		  E = estimated probabilities of transition: E(x,y) = estimate of P(Y=y | X=x) for 0<=x<k
		  l = log-likelihood of y for parameters P and E
	*/
	double HMMtrain(int ** y, VEC(double) nu, double ** P_est, double ** E_est, double tol, int maxIt){
		int i, j, r, it, t;
		// int k;
		r = 0;
		double z = 0;
		double loglik = 0;
		double ll = 0;
		double ll_old = 1;
		double change = 1;
		
		// Estimates for gamma and xi in each episode
		// VEC2D(double) P_tmp(N, VEC(double)(N,0.0)); 
		// VEC2D(double) E_tmp(N, VEC(double)(M,0.0));
		

		double ** fs = (double **) malloc((n+1) * sizeof(double *));
		double ** bs = (double **) malloc((n+1) * sizeof(double *));
		double * c = (double *) malloc((n+1) * sizeof(double));

		// Initialise fs, bs, c
		for(i=0; i<n+1; i++){
			fs[i] = (double *) malloc(N * sizeof(double));
			bs[i] = (double *) malloc(N * sizeof(double));
			c[i] = 0;
		}

		// gamma
		VEC2D(double) A(N, VEC(double)(M,0.0));
		VEC(double) s(N,0.0);
		VEC(double) sum_E(M,0.0);
		// xi
		VEC2D(double) B(N, VEC(double)(N,0.0));
		VEC(double) s2(N,0.0);
		VEC(double) sum_P(N,0.0);

		printf("Initialising estimates\n");
		// Initialise P_tmp, E_tmp, sum variables, and print out initial estimates
		for (i = 0; i < N; i++){
			for (j = 0; j < N; j++){
				// printf("%f ", P[i][j]);
				P_tmp[i][j] = P_est[i][j];
			}
			// printf("\n");
		}

		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
				// printf("%f ", E[i][j]);
				E_tmp[i][j] = E_est[i][j];
			}
			// printf("\n");
		}    

		printf("Training commenced\n");
		// Outer loop: keep iterating until convergence or maxIt reached
		for(it=0; (change > tol) && (it<maxIt); it++){
			change = 0;
			ll_old = loglik;
			loglik = 0;
			// printf("Calculating updates\n");
			// For each episode r compute gamma_{ir}(t) and xi_{ijr}(t)
			for(r = 0; r < seq; r++){
				// Compute P(x(k) | y(1:t)) for all k < t
				ll = HMMfilter(y[r], nu, P_tmp, E_tmp, fs, c);
				HMMsmoother(y[r], P_tmp, E_tmp, c, bs);
				loglik += ll;

				// Initialise values
				fill(s.begin(), s.end(), 0.0);
				fill(s2.begin(), s2.end(), 0.0);

				for(i=0; i<N; i++){
					for(j=0; j<M; j++){
						for(t=0; t<n+1; t++){
							// Find positions in the sequence where the emission == j
							if(y[r][t] == j){
								z = log(fs[t][i]) + log(bs[t][i]);
								A[i][j] += exp(z);
							}
						}
						s[i] += A[i][j];
					}
				}

				// Compute xi_{ij}(t) = P( X_t = i, X_{t+1} = j | Y, theta)
				for(i=0; i<N; i++){
					for(j=0; j<N; j++){
						for(t=0; t<n; t++){			
							B[i][j] += 1/(c[t+1]) * exp(log(fs[t][i]) + log(P_tmp[i][j]) + log(E_tmp[j][y[r][t+1]]) + log(bs[t+1][j]));
						}
						s2[i] += B[i][j];
					}
				}
			}

			// Compute b_i^*(v_k)
			for(i=0; i<N; i++){
				for(j=0; j<M; j++){
					E_tmp[i][j] = A[i][j] / s[i];
				}
			}

			// Compute a_{ij}^*
			for(i=0; i<N; i++){
				for(j=0; j<N; j++){
					P_tmp[i][j] = B[i][j] / s2[i];
				}
			}

			// clean up matrices to deal with zeros and nans
			for(i=0; i<N; i++){	
				// if a row has zero transitions, assume that there are no transitions out of the state
				sum_P[i] = 0;
				for(j=0; j<N; j++){
					sum_P[i] += P_tmp[i][j];
					// Set nan values to zero
					if(isnan(P_tmp[i][j])){
						P_tmp[i][j] = 0;
					}
				}
				if(sum_P[i] == 0){
					for(j=0; j<N; j++){
						P_tmp[i][j] = 0;
					}
					P_tmp[i][i] = 1;
				}

				for(j=0; j < M; j++){
					if(isnan(E_tmp[i][j])){
						E_tmp[i][j] = 0;
					}
				}
			}

			// printf("Updates Computed\n");

			// Compute changes in estimates
			// change = MAX(change,norm_E(E_tmp, E_est)/M);
			change = MAX(change,norm_P(P_tmp, P_est)/N);

			//Compute final estimates by combining estimates from all episodes
			for(i=0; i<N; i++){
				for(j=0; j<N; j++){
					P_est[i][j] = P_tmp[i][j];
				}
				for(j=0; j<M; j++){
					E_est[i][j] = E_tmp[i][j];
				}
			}

			// printf("Updates Applied\n");

			change = MAX(change,fabs(loglik-ll_old)/(1+fabs(ll_old)));
			
			printf("Iteration %d: loglik = %.9f, change = %.9f \n", it, fabs(loglik-ll_old)/(1+fabs(ll_old)), change);
			// printf("Change updated\n");
			// printf("Change: %f\n",change);

			for(i=0; i<N; i++){
				fill(A[i].begin(), A[i].end(), 0);
				fill(B[i].begin(), B[i].end(), 0);
			}
		}
		printf("After\n");
		print_matrices_ptr(P_est, E_est);

		// Free memory

		// printf("done\n");
		
		free(c);
		for(i=0; i<n+1; i++){
			free(fs[i]);
			free(bs[i]);
		}
		free(fs);
		free(bs);

		return loglik;
	}

	
	/*
	compute maximum likehood estimate using Expectation-Maximization - optimised version
	iterations
	in:   y = vector of observations 
		  nu = initial distribution of the hidden chain
		  tol = tolerance for the stopping criterion
		  maxIt = maximal number of iterations
	out:  P = estimate of the transition matrix of the hidden markov process
		  E = estimated probabilities of transition: E(x,y) = estimate of P(Y=y | X=x) for 0<=x<k
		  l = log-likelihood of y for parameters P and E
	*/
	double HMMtrain_opt(int ** y, VEC(double) nu, double ** P_est, double ** E_est, double tol, int maxIt){
		int i, j, k, r, it, t, opt;
		// int k;
		r = 0;
		opt = 0;
		double z = 0;
		double loglik = 0;
		double ll = 0;
		double ll_old = 1;
		double change = 1;
		
		
		// Estimates for gamma and xi in each episode
		// VEC2D(double) P_tmp(N, VEC(double)(N,0.0)); 
		// VEC2D(double) E_tmp(N, VEC(double)(M,0.0));
		double ** fs = (double **) malloc((n+1) * sizeof(double *));
		double ** bs = (double **) malloc((n+1) * sizeof(double *));
		double * c = (double *) malloc((n+1) * sizeof(double));

		// Initialise fs, bs, c
		for(i=0; i<n+1; i++){
			fs[i] = (double *) malloc(N * sizeof(double));
			bs[i] = (double *) malloc(N * sizeof(double));
			c[i] = 0;
		}

		// gamma
		VEC2D(double) A(N, VEC(double)(M,0.0));
		VEC(double) s(N,0.0);
		VEC(double) sum_E(M,0.0);
		// xi
		VEC2D(double) B(N, VEC(double)(N,0.0));
		VEC(double) s2(N,0.0);
		VEC(double) sum_P(N,0.0);

		/*
		* nzr : Count number of non-zero entries in each row of P_est
		* nzc : Count number of non-zero entries in each column of P_est
		* nzr_col: Store the column indices of the non-zero entries in each row of P_est
		* nzc_row: Store the row indices of the non-zero entries in each column of P_est
		*/
		int * nzr = (int *) malloc(N * sizeof(int));
		int * nzc = (int *) malloc(N * sizeof(int));
		int ** nzr_col = (int **) malloc(N * sizeof(int *));
		int ** nzc_row = (int **) malloc(N * sizeof(int *));

		printf("Initialising estimates\n");
		// Initialise P_tmp, E_tmp, sum variables, and print out initial estimates
		for (i = 0; i < N; i++){
			nzr[i] = 0;
			nzc[i] = 0;
			nzr_col[i] = (int *) malloc(N * sizeof(int));
			nzc_row[i] = (int *) malloc(N * sizeof(int));
			for (j = 0; j < N; j++){
				// printf("%f ", P[i][j]);
				P_tmp[i][j] = P_est[i][j];
				nzr_col[i][j] = 0;
				nzc_row[i][j] = 0;
			}
			// printf("\n");
			for (j = 0; j < M; j++){
				// printf("%f ", E[i][j]);
				E_tmp[i][j] = E_est[i][j];
			}
		}

		
		
		printf("Training commenced\n");
		// Outer loop: keep iterating until convergence or maxIt reached
		for(it=0; (change > tol) && (it<maxIt); it++){
			ll_old = loglik;
			loglik = 0;
			// When change is sufficiently small, only compute updates for non-zero entries of transition matrix to save time
			if(change < 1e-3 || opt > 0){
				change = 0;
				// If this is the first time change is < 1e-3, set opt flag for quicker computations
				if(opt == 0){
					opt = 1;
					for(i=0; i<N; i++){
						for(j=0; j<N; j++){
							if(P_est[i][j] > 0){
								nzr_col[i][nzr[i]] = j;
								// printf("%d ", nzr_col[i][nzr[i]]);
								nzr[i]++;
							}
						}
						// printf("\n");
					}
					for(j=0; j<N; j++){
						for(i=0; i<N; i++){
							if(P_est[i][j] > 0){
								nzc_row[j][nzc[j]] = i;
								// printf("%d ", nzc_row[j][nzc[j]]);
								nzc[j]++;
							}
						}
						// printf("\n");
					}
					// printf("\n");
					// // print contents of nzr vector
					// for(i=0; i<N; i++){
					// 	printf("%d ", nzc[i]);
					// }
					// printf("\n");
				}

				for(r = 0; r < seq; r++){
					ll = HMMfilter_opt(y[r], nu, P_tmp, E_tmp, fs, c, nzc, nzc_row, opt);
					HMMsmoother_opt(y[r], P_tmp, E_tmp, c, bs, nzr, nzr_col, opt);
					loglik += ll;

					// Initialise values
					fill(s2.begin(), s2.end(), 0.0);

					// Only update P_tmp for non-zero entries
					for(i=0; i<N; i++){
						for(j=0; j<nzr[i]; j++){
							k = nzr_col[i][j];
							for(t=0; t<n; t++){
								B[i][k] += 1/(c[t+1]) * exp(log(fs[t][i]) + log(P_tmp[i][k]) + log(E_tmp[k][y[r][t+1]]) + log(bs[t+1][k]));
							}
							s2[i] += B[i][k];				
						}
					}
				}
				
				if(opt == 1) opt = 2;

				for(i=0; i<N; i++){
					for(j=0; j<nzr[i]; j++){
						k = nzr_col[i][j];
						P_tmp[i][k] = B[i][k] / s2[i];
					}
				}
				

				// clean up matrices to deal with zeros and nans
				for(i=0; i<N; i++){	
					// if a row has zero transitions, assume that there are no transitions out of the state
					sum_P[i] = 0;
					for(j=0; j<N; j++){
						sum_P[i] += P_tmp[i][j];
						// Set nan values to zero
						if(isnan(P_tmp[i][j])){
							P_tmp[i][j] = 0;
						}
					}
					if(sum_P[i] == 0){
						for(j=0; j<N; j++){
							P_tmp[i][j] = 0;
						}
						P_tmp[i][i] = 1;
					}
				}
				// change = MAX(change,norm_P(P_tmp, P_est)/N);
				// P_est = P_tmp;
			}
			else{
				change = 0;
				// printf("Calculating updates\n");
				// For each episode r compute gamma_{ir}(t) and xi_{ijr}(t)
				for(r = 0; r < seq; r++){
					// Compute P(x(k) | y(1:t)) for all k < t
					ll = HMMfilter(y[r], nu, P_tmp, E_tmp, fs, c);
					HMMsmoother(y[r], P_tmp, E_tmp, c, bs);
					loglik += ll;

					// Initialise values
					fill(s.begin(), s.end(), 0.0);
					fill(s2.begin(), s2.end(), 0.0);

					// MATLAB
					for(i=0; i<N; i++){
						// int same = 0;
						for(j=0; j<M; j++){
							for(t=0; t<n+1; t++){
								// Find positions in the sequence where the emission == j
								if(y[r][t] == j){
									z = log(fs[t][i]) + log(bs[t][i]);
									A[i][j] += exp(z);
								}
							}
							s[i] += A[i][j];
						}
					}

					// Compute xi_{ij}(t) = P( X_t = i, X_{t+1} = j | Y, theta)
					for(i=0; i<N; i++){
						for(j=0; j<N; j++){
							for(t=0; t<n; t++){			
								B[i][j] += 1/(c[t+1]) * exp(log(fs[t][i]) + log(P_tmp[i][j]) + log(E_tmp[j][y[r][t+1]]) + log(bs[t+1][j]));
							}
							s2[i] += B[i][j];
						}
					}
				}
				// printf("Updates Computed\n");

				// Compute b_i^*(v_k)
				for(i=0; i<N; i++){
					for(j=0; j<M; j++){
						E_tmp[i][j] = A[i][j] / s[i];
					}
				}

				// Compute a_{ij}^*
				for(i=0; i<N; i++){
					for(j=0; j<N; j++){
						P_tmp[i][j] = B[i][j] / s2[i];
					}
				}
				
				// clean up matrices to deal with zeros and nans
				for(i=0; i<N; i++){	
					// if a row has zero transitions, assume that there are no transitions out of the state
					sum_P[i] = 0;
					for(j=0; j<N; j++){
						sum_P[i] += P_tmp[i][j];
						// Set nan values to zero
						if(isnan(P_tmp[i][j])){
							P_tmp[i][j] = 0;
						}
					}
					if(sum_P[i] == 0){
						for(j=0; j<N; j++){
							P_tmp[i][j] = 0;
						}
						P_tmp[i][i] = 1;
					}

					for(j=0; j < M; j++){
						if(isnan(E_tmp[i][j])){
							E_tmp[i][j] = 0;
						}
					}
				}
				// printf("Updates Applied\n");	
			}

			// Compute changes in estimates
			change = MAX(change,norm_E(E_tmp, E_est)/M);
			change = MAX(change,norm_P(P_tmp, P_est)/N);

			//Compute final estimates by combining estimates from all episodes
			// E_est = E_tmp;
			// P_est = P_tmp;
			for(i=0; i<N; i++){
				for(j=0; j<N; j++){
					P_est[i][j] = P_tmp[i][j];
				}
				for(j=0; j<M; j++){
					E_est[i][j] = E_tmp[i][j];
				}
			}

			change = MAX(change,fabs(loglik-ll_old)/(1+fabs(ll_old)));
			
			printf("Iteration %d: loglik = %.9f, change = %.9f\n", it, fabs(loglik-ll_old)/(1+fabs(ll_old)), change);
			// printf("Change updated\n");
			// printf("Change: %f\n",change);

			// reset A and B matrices
			for(i=0; i<N; i++){
				fill(A[i].begin(), A[i].end(), 0);
				fill(B[i].begin(), B[i].end(), 0);
			}
		}
		// printf("After\n");
		// print_matrices_ptr(P_est, E_est);

		// Free memory

		// printf("done\n");
		
		free(c);
		for(i=0; i<n+1; i++){
			free(fs[i]);
			free(bs[i]);
		}
		for(i=0; i<N; i++){
			free(nzr_col[i]);
		}
		free(nzr_col);
		free(nzr);
		free(fs);
		free(bs);
		return loglik;
	}


	/*
	Method that executes the Baum-Welch algorithm
	out: ll = double containing the final log-likelihood value
	*/
	double run(double ** P_e, double ** E_e, int uniform_guess = 0, int mode = 0){
		int i,j;
		
		P_est = P_e;
		E_est = E_e;

		// If filename is not empty, read entries of P_est from the file
		if(strcmp(filename.c_str(), "") != 0){
			printf("Reading from file: %s\n", filename.c_str());
			string line;
			ifstream read_file(filename);
			i = 0;
			while(getline(read_file, line)){
				if(i < N){
					double sum = 0.0;
					// Entries separated by space
					stringstream ss(line);
					string s;
					for(j = 0; j < N; j++){
						getline(ss, s, ' ');
						P_est[i][j] = (double) stod(s.c_str());
						sum += P_est[i][j];
					}
					// Normalise P_Est
					for(j = 0; j < N; j++){
						P_est[i][j] /= sum;
					}
				}
				i++;
			}
		}
		else if(uniform_guess == 1){
			// randomTransitionKernel(&P_est);
			uniformTransitionKernel(P_est);
		}
		
		if(N < 400){
			print_matrices_ptr(P_est, E_est);
		}
		
		printf("Running Baum-Welch algorithm\n");

		// Multiple trials
		// VEC2D(int) x_m(seq, VEC(int)(n, 0));
		x_m = (int **) malloc(seq * sizeof(int *));
		for(i=0; i<seq; i++){
			x_m[i] = (int *) malloc(n * sizeof(int));
			for(j=0; j<n; j++){
				x_m[i][j] = 0;
			}
		}

		// VEC2D(int) y_m(seq, VEC(int)(n, 0));
		y_m = (int **) malloc(seq * sizeof(int *));
		for(i=0; i<seq; i++){
			y_m[i] = (int *) malloc(n * sizeof(int));
			for(j=0; j<n; j++){
				y_m[i][j] = 0;
			}
		}
		// VEC2D(int) y_m_ext(seq, VEC(int)(n+1, 0));
		y_m_ext = (int **) malloc(seq * sizeof(int *));
		for(i=0; i<seq; i++){
			y_m_ext[i] = (int *) malloc((n+1) * sizeof(int));
			for(j=0; j<n+1; j++){
				y_m_ext[i][j] = 0;
			}
		}
		

		clock_t begin, end;
		double tps;
		VEC(double) nu(N, 0.0); 
		nu[0] = 1.0;

		P_tmp = (double **) malloc(N * sizeof(double *));
		E_tmp = (double **) malloc(N * sizeof(double *));
		for(i=0; i<N; i++){
			P_tmp[i] = (double *) malloc(N * sizeof(double));
			E_tmp[i] = (double *) malloc(M * sizeof(double));
			for(j=0; j<N; j++){
				P_tmp[i][j] = 0;
			}
			for(j=0; j<M; j++){
				E_tmp[i][j] = 0;
			}
		}

		// Initialise guesses
		// printf("Init\n");
		
		// Sample observations
		// printf("HMM sampling\n");
		int round = 0;
		for(round = 0; round < seq; round++){
			HMMsample(nu, P, E, x_m[round], y_m[round]);
		}
		
		int wins = 0;
		int won = 0;

		// Add an aditional element to the front of each observation set
		for(int i=0; i<seq; i++){
			y_m_ext[i][0] = 0;
			won = 0;
			for(int j=1; j<n+1; j++){
				y_m_ext[i][j] = y_m[i][j-1];
				if(won == 0 && y_m[i][j-1] > M/2 ){
					wins++;
					won = 1;
				}
			}
		}
		
		printf("Win rate: %f\n", (double) wins/seq);

		// //print contents of y_m
		printf("y_m:\n");
		for(int i=0; i<seq; i++){
			for(int j=0; j<n+1; j++){
				printf("%d, ", y_m_ext[i][j]);
			}
			printf("\n");
		}

		printf("HMM training commence\n");
		//Run Baum-Welch algorithm
		begin = clock();
		ll = HMMtrain_opt(y_m_ext, nu, P_est, E_est, tol, maxIt);
		end = clock();

		tps = (end-begin)/(CLOCKS_PER_SEC);
		
		printf("Final estimates\n");
		print_matrices_ptr(P_est, E_est);

		// printf("True Matrices\n");
		// print_matrices(P,E);

		printf("Difference from true P: %.9f\n",diff_P(P_est,P));

		// Print time
		printf("Time: %f s\n", tps);
		
		// free allocated memory
		for(i=0;i<N;i++){
			free(P_tmp[i]);
			P_tmp[i] = NULL;
			free(E_tmp[i]);
			E_tmp[i] = NULL;
		}
		free(P_tmp);
		P_tmp = NULL;
		free(E_tmp);
		E_tmp = NULL;

		// free x_m, y_m, and y_m_ext
		for(i=0;i<seq;i++){
			free(x_m[i]);
			free(y_m[i]);
			free(y_m_ext[i]);
			x_m[i] = NULL;
			y_m[i] = NULL;
			y_m_ext[i] = NULL;
		}
		free(x_m);
		free(y_m);
		free(y_m_ext);
		x_m = NULL;
		y_m = NULL;
		y_m_ext = NULL;
		printf("n = %d\n",n);
		printf("HMM done\n");
		return ll;
	}

protected:
	// VEC2D(double) P_est;
	// VEC2D(double) E_est;
	double ** P_est;
	double ** E_est;
	double ** P_tmp;
	double ** E_tmp;
	int ** x_m;
	int ** y_m;
	int ** y_m_ext;
	VEC2D(double) P;
	VEC2D(double) E;
	double ll;
	double tol;
	int maxIt;
	int N;
	int M;
	int n;
	int seq;
	string filename;
};

#endif // HMM_H