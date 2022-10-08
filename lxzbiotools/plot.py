#!/usr/bin/env python

from re import X
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
#import numpy as np

def stat_distribution(data, window, x_max, x_min=0, mode="freq"):
    num = int((x_max-x_min)/window)+1
 #   print(x_max, x_min,num)
    x = [x_min+i*window for i in range(num)]
    y = []
    for d in data:
        dis_dict = {i: 0 for i in x}
        value = []
        for n in d:
            # smaller < n <= larger
            if n in dis_dict:  
                dis_dict[n] += 1
                continue 
            pos = (int(n/window)+1)*window
            if pos in dis_dict:
                dis_dict[pos] += 1
        y_num = len(d)*1.0
        for i in x:
            if mode == "freq":            
                value.append(100*dis_dict[i]/y_num)
            elif mode == "num":
                value.append(dis_dict[i])
        y.append(value)

    return x, y


def plotDistribution(data, window, label=[], prefix="out", x_max=None, x_min=0, title=None, x_label=None, y_label=None):
    assert len(data) == len(label)

    color = ['r-','b-','g-','k-','m-','c-','y-','r--','b--','g--','k--','m--','c--','y--','r.','b.','g.','k.','m.','c.','y.']
    out = open("%s.xls" % prefix, "w")
    out.write("x\t%s\n" % "\t".join(map(str, label)))
    x, ys = stat_distribution(data, window, x_max, x_min)
    for i in range(len(ys)):
        plt.plot(x, ys[i],color[i], label=label[i])
        
    for i in range(len(x)):
        tmpx = [x[  -++-i]]
        for j in range(len(label)):
            tmpx.append(ys[j][i])
        out.write("\t".join(map(str, tmpx))+"\n")

    out.close()

    plt.xlim(x_min, x_max)
    plt.legend(loc='upper right', shadow=False)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig("%s.pdf" % prefix, format='pdf')
    #plt.savefig("%s.svg" % prefix, format='svg')
    plt.savefig("%s.png" % prefix, format='png', dpi=600 )
    plt.close()


def read_file(fn):
    r = []
    with open(fn) as fh:
        for line in fh.readlines():
            line = line.strip()

            if line.startswith("#"):
                continue
            if line:
                r.append(float(line))

    return r
            

def main():
    import sys
    print("Reading data from file")
    data = read_file(sys.argv[1])
    print("plot image")
    plotDistribution([data], 1, label=[""], prefix="out",title="",x_label="GC content (%)",y_label="Percent (%)",x_min=0, x_max=80)


if __name__ == "__main__":
    main()