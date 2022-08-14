import matplotlib.pyplot as plt
import argparse
import math

# Read data from input file
def read_data(filename):
    data = []
    with open(filename, 'r') as f_in:
        raw_data = f_in.readlines()
    for i in range(1,len(raw_data)):
        try:
            data.append(math.log(float(raw_data[i].strip("\n"))))
        except:
            continue
    return data

# Create plot of data and save as png file
def plot_data(data, filename):
    x = []
    y = []
    for i in range(len(data)):
        x.append(i)
        y.append(data[i])
    plt.plot(x,y)
    plt.title("Change in estimate of transition matrix over BW iterations")
    plt.xlabel("Iteration")
    plt.ylabel(r'log($\Delta P_{est}$)')
    plt.savefig(filename)
    plt.close()
    return

# Main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot data from input file")
    parser.add_argument("-i", "--input", help="Input file consisting of updates to the transition matrix per iteration - one on each line. First entry will be skipped.", required=True)
    args = parser.parse_args()
    # Read data from input file
    data = read_data(args.input)
    # print(data)
    # Create plot of data and save as png file
    plot_data(data, args.input + ".png")