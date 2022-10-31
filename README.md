# TALearner
Codebase for learning task automata

## Requirements

1. g++ compiler

## Instructions for Running Benchmarks

1.	Clone the repository and change directory into it.
2.	Grant permission to your system to run the script using ‘chmod u+x ./run.sh’ 
3.	To run a test on a 3x3 grid world with 3 TA states, simply run ‘./run.sh’. There are a number of additional parameters you can pass in (running ‘./run.sh -h’ will also give you the list of options):\
    -s <size> (default = 3): size of the gridworld (3 = 3x3, 4 =4x4, 5 = 5x5)\
    -d <DFA_states> (default = 3): number of DFA states to use (options are 3,4,5). The 4 state and 5 state ones are essentially just adding additional tasks like ‘pick up coffee, then go to couch, then go up stairs’.\
    -n <number of steps per episode> (default = 34 for 3x3, 80 for 4x4, 85 for 5x5)\
    -q <number of episodes to train on> (default = 275 for 3x3, 500 for 4x4, 2000 for 5x5)\
    -i <maxIterations> (default = 25000)\
    -t <convergence tolerance> (default = 1e-6)\
    -p <parameter that changes the small probability added to each non-zero entry of the guess> (default = 1.0)\
    -f <file containing transition matrix to initialise to> (default = “”): specify a file containing a transition matrix with spaces between entries and newline for each row.\
    -m <mode> (default = 0): Flag for choosing to directly learn the full product MDP (-m 0) or to use the two-stage pipeline (-m 1).\
    -o <output format> (default = 0): 0: write to log only, 1: print to console + write to log, 2: print to console only\

## Example

To run the 5x5 with 4 DFA states, 125 steps per episode, and 3000 sequences, enter the command ‘./run.sh -s 5 -d 4 -n 125 -q 3000'. By default, the output will be written to a file in the logs folder, named by the date and time when the command was executed.

