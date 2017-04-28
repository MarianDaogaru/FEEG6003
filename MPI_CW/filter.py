import os
import numpy
import scipy
from matplotlib import pyplot as plt

#-----VARS
val_space="="


def get_pbs_files():
    paths = numpy.array(os.listdir(os.getcwd()))
    return paths[[i for i, item in enumerate(paths) if "pbs" in item]]


def rename_files():
    paths = get_pbs_files()
    main_path = os.getcwd()+"/"

    for path in paths:
        data = open_pbs_file(path)
        data = [line for line in data if line[0] == "init"][0]
        name = "{}x{}_{}x{}_{}_{}_{}_{}.pbs".format(data[4].split(val_space)[1],
                                        data[5].split(val_space)[1],
                                        data[2].split(val_space)[1],
                                        data[3].split(val_space)[1],
                                        data[1].split(val_space)[1],
                                        data[6].split(val_space)[1],
                                        data[7].split(val_space)[1],
                                        data[8].split(val_space)[1])
        os.rename(main_path+path, main_path+name)
    return None


def open_pbs_file(path, how="r"):
    with open(path, how) as data:
        dt = data.read()
    del(data)
    return numpy.array([numpy.array([j for j in i.split(" ")]) for i in dt.split("\n")[:-1]])


def plot_delta(path):
    data = open_pbs_file(path)
    data = numpy.array([a for a in data if a[0].split("=")[0]=="max_delta"])

    max_delta = numpy.array([eval(val.split(val_space)[1]) for val in data[:,0]])
    iterations = numpy.array([eval(val.split(val_space)[1]) for val in data[:,1]])

    plt.plot(iterations, max_delta, "b")
    plt.plot(iterations, max_delta, "bo", label=r"max $\Delta$")
    plt.xlabel("iterations")
    plt.ylabel("Maximum $\Delta$")
    plt.title(r"Maximum $\Delta$ using {} processes for $\Delta$ limit {}.".format(data[0,3].split(val_space)[1], data[0,4].split(val_space)[1]))
    plt.grid(which="both", axis="both")
    plt.ylim(ymax=2)
    plt.axhline(y=eval(data[0,4].split(val_space)[1]), xmin=0, xmax=1000, color="k")
    plt.legend()
    plt.savefig("max_delta_{}_{}.jpeg".format(data[0,3].split(val_space)[1], data[0,4].split(val_space)[1]))
    return max_delta


def plot_avg(path):
    data = open_pbs_file(path)
    data = numpy.array([a for a in data if a[0].split(val_space)[0]=="local_avg"])

    local_avg = numpy.array([eval(val.split(val_space)[1]) for val in data[:,0]])
    iterations = numpy.array([eval(val.split(val_space)[1]) for val in data[:,1]])

    plt.plot(iterations, local_avg, "b")
    plt.plot(iterations, local_avg, "bo", label="local average")
    plt.xlabel("iterations")
    plt.ylabel("Local Average")
    plt.title(r"Average average at using {} processes for $\Delta$ limit {}.".format(data[0,2].split(val_space)[1], data[0,3].split(val_space)[1]))
    plt.grid(which="both", axis="both")
    plt.legend(loc=2)
    plt.savefig("local_avg_{}_{}.jpeg".format(data[0,2].split(val_space)[1], data[0,3].split(val_space)[1]))
    return local_avg


def get_times(path):
    data = open_pbs_file(path)
    data = numpy.array([line for line in data if line[0].split(val_space)[0]=="avg_time"])

    times = numpy.zeros((data.shape[0], data.shape[1]-1))
    for line in data:
        pos = eval(line[1].split(val_space)[1])
        times[pos, 0] = eval(line[0].split(val_space)[1])
        for i in range(2, line.shape[0]):
            times[pos, i-1] = eval(line[i].split(val_space)[1])
    return times


def plot_times_one_file(path):
    times = get_times(path)
    overall = times[:, 1]
    loops = times[:, 2]
    times = scipy.delete(times, 1, 1)
    times = scipy.delete(times, 1, 1)
    lgth = times.shape[1]
    terr = times[times.argmax(axis=0), numpy.arange(0, lgth)] - \
           times[times.argmin(axis=0), numpy.arange(0, lgth)]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(numpy.arange(0, lgth), times.mean(axis=0), fmt="x", yerr=terr, label="timing", ms=10, mew=2)
    ax.set_xlim(xmin=-1, xmax=lgth)
    ax.set_xlabel("Different timings: avg loop, make MP, neighboors, buf, reconstruct, barrier")
    ax.set_title("Max operating time {:1.6f}, average operating time {:1.6f}".format(max(overall), overall.mean()))
    ax.set_ylabel("Times(s)")
    ax.grid(which="both", axis="both")
    plt.legend()
#    plt.grid(which="both", axis="both")
    plt.savefig("{}.jpeg".format(path[:-4]), format="jpeg")
    plt.close()


def get_final_res(path):
    data = open_pbs_file(path)
    data = numpy.array([a for a in data if a[0].split("=")[0]=="global_avg"])

    print(r"Program finished with a global average of {}, at maximum $\Delta$={}, after {} iterations".format(data[0,0].split(val_space)[1], data[0,1].split(val_space)[1], data[0,2].split(val_space)[1]))

if __name__ == "__main__":
    rename_files()
    path = get_pbs_files()[0]
#    dt = open_pbs_file(path)
#    avg = plot_avg(path)
#    delta = plot_delta(path)
    plot_times_one_file(path)
    get_final_res(path)
