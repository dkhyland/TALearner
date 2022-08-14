#!/bin/bash
filename="$1"
data=$(grep 'change' $filename)

# Get the right folder and unique filename
if [ ! -d data/$(date "+%Y-%m-%d") ]; then
    mkdir data/$(date "+%Y-%m-%d")
fi
num_files=$(ls data/$(date "+%Y-%m-%d") | wc -l)

write_data=""
# Extract change in transition matrix
while IFS= read -r line; do
    # Transition matrix updates
    text=$(echo $line |  sed -n -e 's/^.*change = //p')
    write_data="${write_data}${text}
"
done <<< "$data"

echo "$write_data" >>  data/"$(date "+%Y-%m-%d")"/"P_${num_files}_$(date "+%Y-%m-%d_%H:%M:%S")".log