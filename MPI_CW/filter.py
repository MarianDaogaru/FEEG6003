import os
import numpy
from matplotlib import pyplot as plt

#-----VARS
val_space="="


def get_pbs_files():
    paths = numpy.array(os.listdir(os.getcwd()))
    return paths[[i for i, item in enumerate(paths) if "pbs" in item]]


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
    data = numpy.array([a for a in data if a[0].split("=")[0]=="local_avg"])

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
    pass


def get_final_res(path):
    data = open_pbs_file(path)
    data = numpy.array([a for a in data if a[0].split("=")[0]=="global_avg"])

    print(r"Program finished with a global average of {}, at maximum $\Delta$={}, after {} iterations".format(data[0,0].split(val_space)[1], data[0,1].split(val_space)[1], data[0,2].split(val_space)[1]))

if __name__ == "__main__":
    path = get_pbs_files()[0]
#    dt = open_pbs_file(path)
#    avg = plot_avg(path)
#    delta = plot_delta(path)

    get_final_res(path)
