import numpy
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['figure.figsize'] = (15, 7.5)
rcParams['font.size'] = 15

l1 = [[0.200321,0.200321,0.200321,0.200321,0.200321,0.200321,0.200321],
      [0.015195,0.015195,0.014767,0.015195,0.014770,0.014765,0.014773],
       [0.015064,0.015064,0.014804,0.015064,0.014805,0.014798,0.014800],
       [0.033237,0.033231,0.033176,0.034635,0.035951,0.039639,0.047092],
       [0.033477,0.033449,0.032876,0.033504,0.032999,0.033024,0.036094],
        [0.032929,0.032868,0.032347, 0.032921,0.032563,0.033550,0.032518]
    ]

l2 = [[0.617023,0.617023,0.617023,0.617023,0.617023,0.617023,0.617023],
      [0.405484,0.405484,0.405484,0.405484,0.405484,0.405484,0.405484],
        [0.405642,0.405642,0.405642,0.405642,0.405642,0.405642,0.405642],
        [0.266813,0.193197,0.126849,0.138048,0.194990,0.330853,0.389248],
        [0.110584,0.110625,0.108558,0.110520,0.155836,0.291869,0.369729],
        [0.339318,0.339361,0.333666,0.3392613,0.333779,0.333601,0.340930]]



lx = [1, 2, 4, 8, 16, 32, 64]

data_graph_points = [["b", "bo"],
                     ["r", "ro"],
                     ["g", "go"],
                     ["k", "ko"],
                        ["b", "bx"],
                        ["r", "rx"]]

names = ["init", "static", "auto", "static_n", "dynamic_n", "guided_n"]


def plot_draph(l, t):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    n = 3
    for i in range(len(l)-n):
        print(lx, l[i+n])
        ax.plot(lx, l[i+n], data_graph_points[i][0])
        ax.plot(lx, l[i+n], data_graph_points[i][1], label=names[i+n])

    ax.set_xlabel("chunksize")
    ax.set_ylabel("time (s)")
    ax.set_title("time vs chunksize for {}".format(t))
    ax.grid(which="both", axis="both")
    plt.legend()
    plt.savefig("{}.jpeg".format(t))
    plt.close()

