import numpy
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['figure.figsize'] = (15, 10)
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



l3 = numpy.array([0.045269, 0.034604, 0.026104, 0.014761, 0.007732, 0.004024])
l4 = numpy.array([0.615633, 0.301664, 0.202422, 0.108431, 0.054432, 0.038987])

lx_2 = [1, 2, 3, 6, 12, 24]

l5 = numpy.array([0.191396,0.109869,0.082297,0.057153,0.042846,0.066724])
l6 = numpy.array([0.518554,0.383360,0.329376,0.190079,0.062579,0.032091])

l8 = [19.195229, 6.705055]
l9 = [52.099636, 3.262320]

data_graph_points = [["b", "bo"],
                     ["r", "ro"],
                     ["g", "go"],
                     ["k", "ko"],
                        ["b", "bx"],
                        ["r", "rx"]]

names = ["init", "static", "auto", "static_n", "dynamic_n", "guided_n"]

def plot_graph():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(lx_2, l3[0]/l3, "b")
    ax.plot(lx_2, l3[0]/l3, "bo", label="loop1 static")
    ax.plot(lx_2, l4[0]/l4, "g")
    ax.plot(lx_2, l4[0]/l4, "go", label="loop2 dynamic(4)")

    ax.set_xlabel("number of threads")
    ax.set_ylabel("speed up T1/T")
    ax.set_title("Speed up for different number of threads")
    ax.grid(which="both", axis="both")
    plt.legend(loc=4)
    plt.savefig("speed_up.jpeg")
    plt.close()

def plot_graph2():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(lx_2, l5[0]/l5, "b")
    ax.plot(lx_2, l5[0]/l5, "bo", label="loop1")
    ax.plot(lx_2, l6[0]/l6, "g")
    ax.plot(lx_2, l6[0]/l6, "go", label="loop2")

    ax.set_xlabel("number of threads")
    ax.set_ylabel("speed up T1/T")
    ax.set_title("Speed up for different number of threads using affinity")
    ax.grid(which="both", axis="both")
    plt.legend(loc=4)
    plt.savefig("speed_up_af.jpeg")
    plt.close()

def plot_graph3():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(lx_2, l3[0]/l3, "b")
    ax.plot(lx_2, l3[0]/l3, "bo", label="loop1")
    ax.plot(lx_2, l5[0]/l5, "g")
    ax.plot(lx_2, l5[0]/l5, "go", label="loop1 with affinity")

    ax.set_xlabel("number of threads")
    ax.set_ylabel("speed up T1/T")
    ax.set_title("Speed up for different number of threads for loop1")
    ax.grid(which="both", axis="both")
    plt.legend(loc=4)
    plt.savefig("speed_up_l1.jpeg")
    plt.close()

def plot_graph4():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(lx_2, l4, "b")
    ax.plot(lx_2, l4, "bo", label="loop2")
    ax.plot(lx_2, l6, "g")
    ax.plot(lx_2, l6, "go", label="loop2 with affinity")

    ax.set_xlabel("number of threads")
    ax.set_ylabel("speed up T1/T")
    ax.set_title("Execution times for different number of threads for loop2")
    ax.grid(which="both", axis="both")
    plt.legend(loc=3)
    plt.savefig("speed_up_l22.jpeg")
    plt.close()


def plot_draph(l, t, loc):
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
    plt.legend(loc=loc)
    plt.savefig("{}.jpeg".format(t))
    plt.close()

