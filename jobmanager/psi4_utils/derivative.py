#Moved to derivative_utils.py

import os
import operator
from typing import List


def get_subfolders(path: str = "./") -> list:
    files = os.listdir(path)
    folders = []
    for file in files:
        if os.path.isdir(path + "/" + file):
            folders.append(file)
    return folders


def sanity_check(folders: list, trigger: str = "_derivNo_") -> dict:
    trigger_appear = False
    basename: List[str] = list()
    folders_ignored = list()
    jobs = dict()
    for folder in folders:
        if trigger in folder:
            trigger_appear = True
            base = folder.split(trigger)[0]
            # print(base, folder)
            if not len(basename):
                first_base = base
                basename.append(base)
            if base not in basename:
                basename.append(base)
            if base == first_base:
                jobs[folder] = int(folder.split(trigger)[-1])
        else:
            folders_ignored.append(folder)
    if not trigger_appear:
        raise KeyError("<trigger>: %s not appeared in any folders." % trigger)
    if len(basename) > 1:
        raise ValueError("multiple base jobs occur: %s" % basename)
    if len(folders_ignored):
        print("Warning: The following folders are ignored: %s" % str(folders_ignored))
    return jobs


def derivative_tree(path: str = "./", trigger: str = "_derivNo_") -> list:
    '''
    Make a sequence structure for derivative jobs in a path.
    '''
    folders = get_subfolders(path)
    jobs = sanity_check(folders, trigger=trigger)
    jobs = dict(sorted(jobs.items(), key=operator.itemgetter(0)))
    return list(jobs.keys())

#Altered in derivative_utils so that derivative calculations are executed from the parent directory
def get_wfn_path(jobs, ii):
    assert ii > 0
    return "../../" + jobs[ii - 1] + "/b3lyp/wfn.180.npy"
