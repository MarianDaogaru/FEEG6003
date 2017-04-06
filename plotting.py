import numpy
import matplotlib.pyplot as plt


d1 = {}

d2 = {}

l1 = [1, 2, 4, 8, 16, 32, 64]

data_graph_points = [["b", "bo"],
                     ["r", "ro"],
                     ["g", "go"],
                     ["k", "ko"],
                        ["b", "bx"],
                        ["r", "rx"]]

names = ["init", "static", "auto", "static_n", "dynamic_n", "guided_n"]
def plot_draph(d, t):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    keys = d.keys
    for i in range(len(keys)):
        ax.plot(l1, d[keys[i]], data_graph_points[i][0])
        ax.plot(l1, d[keys[i]], data_graph_points[i][1], label=names[i])

    ax.set_xlabel("chunksize")
    ax.set_ylabel("time (s)")
    ax.set_title("time vs chunksize for {}".format(t))
    ax.grid(which="both", axis="both")
    plt.legend()
    plt.savefig("{}.jpeg".format(t))
    plt.close()

