import os
import numpy
import scipy
import time


def get_pgm_files():
    paths = numpy.array(os.listdir(os.getcwd()))
    return paths[[i for i, item in enumerate(paths) if "pgm" in item]]


def open_pgm(path):
    print(path)
    with open(path) as dt:
        data = dt.readlines()
    del(dt)
    data = data[4:]
    #data_2 = numpy.zeros((len(data), len(data[0].split("\n")[0].split(" ")[:-1])))
    data_2 = []
    cont = len(data[0].split("\n")[0].split(" ")[:-1])
    for i in range(len(data)):
        for j in range(len(data[i].split("\n")[0].split(" ")[:-1])):
            if data[i].split("\n")[0].split(" ")[j] != "":
                data_2.append(float(data[i].split("\n")[0].split(" ")[j]))
    return data_2


if __name__=="__main__":
    p = get_pgm_files()
    d = open_pgm(p[25])