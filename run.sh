#!/bin/bash
cd src/C++
g++ -Werror -Wall -O3 main.cpp -o main
if [[ $? == 0 ]]; then
    if [ ! -d ../../logs/$(date "+%Y-%m-%d") ]; then
        mkdir ../../logs/$(date "+%Y-%m-%d")
    fi

    num_logs=$(ls ../../logs/$(date "+%Y-%m-%d") | wc -l)

    if [ "$#" -eq 0 ]; then
        ./main > ../../logs/"$(date "+%Y-%m-%d")"/"out_${num_logs}_$(date "+%Y-%m-%d_%H:%M:%S")".log
    else
        out=0
        mode=0
        args=""
        while getopts s:d:n:q:i:t:p:o:i:f:m:hl flag
        do
            case "${flag}" in
                s) args="${args} -size ${OPTARG}";;
                d) args="${args} -d ${OPTARG}";;
                n) args="${args} -n ${OPTARG}";;
                q) args="${args} -seq ${OPTARG}";;
                i) args="${args} -maxIt ${OPTARG}";;
                t) args="${args} -tol ${OPTARG}";;
                p) args="${args} -init ${OPTARG}";;
                f) args="${args} -file ${OPTARG}";;
                m) args="${args} -mode ${OPTARG}";;
                o) out=${OPTARG};;
                h) echo "Usage: run.sh [-s <size>] [-d <DFA_states>] [-n <n>] [-q <seq>] [-i <maxIt>] [-t <tol>] [-p <init_P_est>][-f <filename>] [-m <mode>] [-o <output format> (0: write to log, 1: print + write to log, 2: print only)]"; exit;;
            esac
        done
        if [ ${out} -eq 0 ]; then
            ./main ${args} > ../../logs/$(date "+%Y-%m-%d")/"out_${num_logs}_$(date "+%Y-%m-%d_%H:%M:%S")".log
        elif [ ${out} -eq 1 ]; then
            ./main ${args} | tee ../../logs/$(date "+%Y-%m-%d")/"out_${num_logs}_$(date "+%Y-%m-%d_%H:%M:%S")".log
        elif [ ${out} -eq 2 ]; then
            ./main ${args}
        fi
    fi
fi
echo "Done"
cd ../..