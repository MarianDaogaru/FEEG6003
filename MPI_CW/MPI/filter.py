import os
import numpy
import scipy
import time
from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['figure.figsize'] = (12, 12)
rcParams['font.size'] = 15
#-----VARS
val_space="="


def get_pbs_files():
    paths = numpy.array(os.listdir(os.getcwd()))
    return paths[[i for i, item in enumerate(paths) if "pbs" in item]]


def file_splitter(path):
    with open(path, 'r') as dt:
        data = dt.read()
    data = data.split("Application")


    print(len(data))
    with open(path, 'w') as dt:
        dt.write(data[0])

    if len(data) > 2: #1 will be the initial, 2 will be the residual
        for i in range(1, len(data)-1):
            with open(str(time.time())+".pbs", 'w') as dt:
                dt.write(data[i])


def rename_files():
    paths = get_pbs_files()
    main_path = os.getcwd()+"/"


    for path in paths:
        data = open_pbs_file(path)
        data = [line for line in data if line[0] == "init"][0]
        print(path, data)
        name = "{}x{}_{}x{}_{}_{}_{}_{}.pbs".format(data[4].split(val_space)[1],
                                        data[5].split(val_space)[1],
                                        data[2].split(val_space)[1],
                                        data[3].split(val_space)[1],
                                        data[1].split(val_space)[1],
                                        str(eval(data[6].split(val_space)[1])),
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

    print(data.shape, path)
    times = numpy.zeros((data.shape[0], data.shape[1]-1))
    for line in data:
        pos = eval(line[1].split(val_space)[1])
        times[pos, 0] = eval(line[0].split(val_space)[1])
        for i in range(2, line.shape[0]):
            times[pos, i-1] = eval(line[i].split(val_space)[1])
    return times

def get_times_serial(path):
    data = open_pbs_file(path)
    data = numpy.array([line for line in data if line[0].split(val_space)[0]=="iter"])

    print(data.shape, path, data)
    return eval(data[0,1].split(val_space)[1])

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
    plt.savefig("{}.jpeg".format(path[:-4]), format="jpeg")
    plt.close()


def plot_execution_times(var, which_val=0):
    which_a = ["768x768_256x384_6_", "768x768_-1x-1_-1_"]
    which = which_a[which_val]
    titles = ["number of processes", "maximum $\Delta$","$\Delta$ frequency", "average calculation frequency"]
    endings = ['', '_100_200.pbs', '_', '.pbs']
    def_val = ['','{}_100_200.pbs','0.05_{}_200.pbs','0.05_100_{}.pbs']
    if var == 2:
        pass
    elif (var == 4 or var == 5):
        ranking = [1, 2, 5, 10, 25, 50, 75, 100, 150, 200, -1]
    else:
        #ranking = [1.0, 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.025, 0.01]
        ranking = [0.1, 0.075, 0.05, 0.025, 0.01]

    paths = []
    for r in ranking:
        paths.append(which+def_val[var-2].format(r))
    paths = numpy.array(paths)
#    print(paths, type(paths))
    t_d = get_times(paths[0])

    times = numpy.zeros((paths.shape[0], t_d.shape[0], t_d.shape[1]))
    for i in range(paths.shape[0]):
        times[i] = get_times(paths[[pa for pa, item in enumerate(paths) if "_"+str(ranking[i])+endings[var-2] in item]][0])
#    print(times)
#
    lgth = paths.shape[0]
    terr = numpy.zeros(lgth)
    avg = numpy.zeros(lgth)
    for i in range((lgth)):
        terr[i] = max(times[i, :, 1]) - min(times[i, :, 1])
        avg[i] = times[i, :, 1].mean()


    terr =  terr[-1] / terr[:-1]
    avg_0 = avg[-1]
    avg = avg[-1] / avg[:-1]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(ranking[:-1], avg, fmt="x", yerr=terr, label="timing", ms=10, mew=2)
    ax.set_xlim(xmax=max(ranking)*1.01)
    ax.set_xlabel(r"$\Delta$ values")
    ax.set_ylabel("Time factor T0/T, T0={:.4f}s".format(avg_0))
    ax.grid(which="both", axis="both")
    plt.legend()
    plt.savefig("iter vs {} 2.jpeg".format(var), format="jpeg")
    plt.close()



def plot_serial(var):
    which = "768x768_-1x-1_-1_"
    titles = ["number of processes", "maximum $\Delta$","$\Delta$ frequency", "average calculation frequency"]
    endings = ['', '_100_200.pbs', '_', '.pbs']
    def_val = ['','{}_100_200.pbs','0.05_{}_200.pbs','0.05_100_{}.pbs']

    if var == 2:
        pass
    elif (var == 4 or var == 5):
        ranking = [1, 2, 5, 10, 25, 50, 75, 100, 150, 200, -1]
    else:
        ranking = [1.0, 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.025, 0.01]

    paths = []
    for r in ranking:
        paths.append(which+def_val[var-2].format(r))
    paths = numpy.array(paths)
    print(paths)

    times = numpy.zeros(paths.shape[0])
    for i in range(paths.shape[0]):
        times[i] = get_times_serial(paths[i])
    print(times)

    times_0 = times[-1]
    times = times[:-1] / times[-1]
    plt.plot(ranking[:-1], times, "bx", label="timing", ms=10, mew=2)
    plt.xlabel(r"iterations")
    plt.ylabel("Time speed-up T0/T, T0={:.4f}s".format(times_0))
    plt.grid(which="both", axis="both")
    plt.xlim(xmax=max(ranking)*1.01)
    plt.legend()
    plt.savefig("serial vs {}.jpeg".format(var), format="jpeg")
    plt.close()



def get_final_res(path):
    data = open_pbs_file(path)
    data = numpy.array([a for a in data if a[0].split("=")[0]=="global_avg"])

    print(r"Program finished with a global average of {}, at maximum $\Delta$={}, after {} iterations".format(data[0,0].split(val_space)[1], data[0,1].split(val_space)[1], data[0,2].split(val_space)[1]))


def plot_multiple(var, N):
    all_paths = get_pbs_files()
    paths = []
    for path in all_paths:
        with open(path) as dt:
            data = dt.read()
        if "M={}".format(var) in data:
            paths.append(path)

    times = numpy.zeros(30)
    processes = numpy.zeros(30)

    for path in paths:
        if "serial" in path:
            with open(path) as dt:
                data = dt.readlines()
            del(dt)
            for line in data:
                if line.split(" ")[0].split("=")[0] == "iter":
                    time_single = eval(line.split(" ")[1].split("=")[1])
        else:
            with open(path) as dt:
                data = dt.read()
            del(dt)
            data = data.split("Application")

            for i in range(30):
                dt = data[i].split("\n")
                max_overall_time = 0
                for line in dt:
                    if line.split(" ")[0].split("=")[0] == "avg_time":
                        max_local = eval(line.split(" ")[2].split("=")[1])
                        if max_local > max_overall_time:
                            max_overall_time = max_local
                    elif line.split(" ")[0].split("=")[0] == "init":
                        processes[i] = eval(line.split(" ")[1].split("=")[1])
                times[i] = max_overall_time

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
    print(time_single/times)
    ax0.plot(processes, times, "bx", label="timing", ms=10, mew=2)
    ax0.set_xlabel("processes")
    ax0.set_ylabel("Time (s)")
    ax0.set_title("Overall execution time for the picture M={}, N={}".format(var, N))
    ax0.grid(which="both", axis="both")

    ax1.plot(processes, time_single/times, "bx", label="timing", ms=10, mew=2)
    ax1.set_xlabel("processes")
    ax1.set_ylabel("Time speed-up T0/T, T0={:.4f}s".format(time_single))
    ax1.set_title("Relative execution time for the picture M={}, N={}".format(var, N))
    ax1.grid(which="both", axis="both")
    plt.legend()
    plt.savefig("exec_{}x{}.jpeg".format(var, N), format="jpeg")
    plt.close()

if __name__ == "__main__":
    pics = [[192, 128], [256, 192], [512, 384], [768, 768]]
    for pic in pics:
        plot_multiple(pic[0], pic[1])
#    rename_files()
#    plot_execution_times(5)
#    paths = get_pbs_files()
##    dt = open_pbs_file(path)
##    avg = plot_avg(path)
##    delta = plot_delta(path)
##    plot_times_one_file(path)
#    for path in paths:
#        print(path)
#        get_final_res(path)
