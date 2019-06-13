#!/usr/bin/env python

import csv
import numpy as np
NAME= 'kurde.csv'

handle = open(NAME, 'r')
def parse():
    rows = csv.reader(handle, delimiter=',')
    domains = next(rows)[1:]
    rows = list(rows)
    proteins = [el[0] for el in rows]
    rows = np.array([el[1:] for el in rows])

if __name__ == '__main__':
    parse()