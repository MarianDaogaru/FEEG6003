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
    plt.plot(iterations, max_delta, "bo", label="max delta")
    plt.xlabel("iterations")
    plt.ylabel("Maximum Delta")
    plt.title("Maximum delta using {} processes for limit .".format(data[0,3].split(val_space)[1], data[0,4].split(val_space)[1]))
    plt.grid(which="both", axis="both")
    plt.ylim(ymax=2)
    plt.axhline(y=eval(data[0,4].split(val_space)[1]), xmin=0, xmax=1000, color="k")
    plt.legend()
    plt.savefig("max_delta_{}_{}.jpeg".format(data[0,3].split(val_space)[1], data[0,4].split(val_space)[1]))
    return max_delta



if __name__ == "__main__":
    dt = open_pbs_file(get_pbs_files()[0])
    data = plot_delta(get_pbs_files()[0])

