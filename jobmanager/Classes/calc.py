import pickle
import os
import numpy as np
import json
from abc import ABC, abstractmethod


def strip_new_line(string):
    """Tries to strip string of new line.

        Parameters
        ----------
            string : str
                Input string.

        Returns
        -------
            output : str
                Output string with newline characters removed..

    """
    if string[-1] == '\n':
        return string[:-1]
    else:
        return string


class Calc(ABC):
    """
    Interface with software/write inputs
    Contains information about a single calculation
        - type (single point, minimization, TS search)
    ###for now just copy-pasted from Readme
    """
    def __init__(self):

        self.run = None
        self.coordinates = None
        self.spinmult = 1
        self.charge = 0
        self.method = None
        self.basis = None
        self.opt_dict = {}
        # ### props ####
        # self.ml_prop = 0
        # self.poptype = 0
        # self.bond_order_list = 0

    @abstractmethod
    def read_input(self,input_path):
        ...


class TerachemCalc(Calc):

    def read_input(self,input_path):
        '''
        reads terachem input file

            Parameters
            ----------
                input_path : str
                    absolute path to the input file
        '''
        with open(input_path, 'r', errors='ignore') as fin:
            input_textfile = fin.readlines()

        input_textfile = [strip_new_line(i) for i in input_textfile]

        track_list = []
        ref = ['run', 'coordinates','spinmult','charge','method','basis']
        req_dict = {} #dictionary of required properties
        lines = input_textfile
        for line in lines:
            # remove whitespace
            line = line.strip()
            if not line: continue  # skip blank lines
            if line.startswith('#'): continue  # skip comments
            if line.startswith('end'): break  # end
            # skip block values
            if line.startswith('$'):
                while not line.startswith('$'):
                    line = next(lines)
            else:
                key, val = line.split(None, 1)
                if key in track_list:
                    print(f'{key:s} specified multiple times. Ignoring value {val:s}')
                else:
                    if key in ref:
                        req_dict[key] = val
                    else:
                        self.opt_dict[key] = val
                    track_list.append(key)

        # second read through for block options
        # lines = iter(input_textfile.lines)
        for line in lines:
            # remove whitespace
            line = line.strip()
            if not line.startswith('$'):
                continue
            else:
                key = line
                val = []
                line = next(lines)
                while not line.startswith('$'):
                    val.append(str(line))
                    line = next(lines)
                if key in ref:
                    req_dict[key] = val
                else:
                    self.opt_dict[key] = val
                track_list.append(key)

        self.run = req_dict['run']
        self.coordinates = req_dict['coordinates']
        self.spinmult = req_dict['spinmult']
        self.charge = req_dict['charge']
        self.method = req_dict['method']
        self.basis = req_dict['basis']