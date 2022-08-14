/* 
 * Main file for running benchmarks
*/

#include "benchmark.h"

using namespace std;

/*
Optional arguments to main:
	-h, --help: Prints help message
	-s, --size: Specifies grid size
	-n, --num_trials: Specifies number of trials
*/

int main(int argc, char ** argv){
	Benchmark * bench = new Benchmark();
	// Default arguments
	int grid_size = 3;
	int DFA_states = 3;
	int n = 34;
	int seq = 275;
	int maxIt = 25000;
	double tol = 1e-6;
	double init_P_est = 1.0;
	int mode = 0;
	string filename = "";

	if(argc < 2){
		// Run with default arguments
		bench->run();
	}
	else{
		for(int i=1; i < argc; i += 2){
			if(strcmp(argv[i], "-size") == 0 && i+1 < argc){
				grid_size = atoi(argv[i+1]);
				// Adjust episode length and number of episodes for different grid sizes
				if(grid_size == 4){
					n = 80;
					seq = 500;
				}
				else if(grid_size == 5){
					n = 85;
					seq = 2000;
				}
				else if(grid_size == 6){
					n = 750;
					seq = 4000;
				}
			}
			else if(strcmp(argv[i], "-d") == 0 && i+1 < argc){
				DFA_states = atoi(argv[i+1]);
				if(DFA_states != 3 && DFA_states != 4 && DFA_states != 5){
					printf("Invalid number of DFA states: %d\n",DFA_states);
					return 0;
				}
			}
			else if(strcmp(argv[i], "-n") == 0 && i+1 < argc){
				n = atoi(argv[i+1]);
			}
			else if(strcmp(argv[i], "-seq") == 0 && i+1 < argc){
				seq = atoi(argv[i+1]);
			}
			else if(strcmp(argv[i], "-maxIt") == 0 && i+1 < argc){
				maxIt = atoi(argv[i+1]);
			}
			else if(strcmp(argv[i], "-tol") == 0 && i+1 < argc){
				tol = (double) stod(argv[i+1]);
			}
			else if(strcmp(argv[i], "-init") == 0 && i+1 < argc){
				init_P_est = (double) stod(argv[i+1]);
			}
			else if(strcmp(argv[i], "-file") == 0 && i+1 < argc){
				filename = argv[i+1];
			}
			else if(strcmp(argv[i], "-mode") == 0 && i+1 < argc){
				mode = atoi(argv[i+1]);
				if(mode > 1){
					printf("Invalid mode: %d\n",mode);
					return 0;
				}
			}
			else{
				printf("Invalid argument: %s\n", argv[i]);
				return 1;
			}
		}
		if(mode == 0){
			bench->run(grid_size,DFA_states,n,seq,maxIt,tol,init_P_est,filename.c_str(), mode);
		}
		else{
			if(init_P_est == 0){
				init_P_est = 1.0;
			}
			bench->run(grid_size,DFA_states,n,seq,maxIt,tol,init_P_est,filename.c_str(), mode);
		}
		
	}
	printf("Done\n");
	return 0;
}