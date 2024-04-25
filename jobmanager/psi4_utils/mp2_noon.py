from jobmanager.psi4_utils.run import lacvps, get_molecule
from numpy.linalg import eigh
import psi4
import shutil
import os
import json

def run_mp2_noon(psi4_config):
    basedir =   os.getcwd()
    rundir  =   "mp2_noon"
    if not os.path.isdir(rundir):
        os.makedirs(rundir)
    os.chdir(rundir)

    ## Set psi4 options
    with open(psi4_config["charge-spin-info"]) as fil:
        psi4_config.update(json.load(fil))
    psi4.core.be_quiet()
    psi4.qcdb.libmintsbasisset.basishorde['LACVPS'] = lacvps
    psi4.set_options({
        'puream': False,
        'basis': psi4_config["basis"],
        'reference':'uhf',
        'maxiter':500
    })
    psi4.set_memory(psi4_config["memory"])
    psi4.set_num_threads(psi4_config["num_threads"])
    psi4.core.set_output_file(psi4_config["output"], False)
    molecule = get_molecule(
        psi4_config["xyzfile"],
        psi4_config["charge"],
        psi4_config["spin"]
    )
    kwargs = {'return_wfn':True,'molecule':molecule}
    if os.path.exists(psi4_config["wfnfile"]):
        # kwargs.update({'restart_file':psi4_config["wfnfile"],'guess':'read'})  # only works with psi4 v1.5+
        filename = psi4_config["output"].split('.')[0]
        pid = str(os.getpid())
        targetfile = filename + '.default.' + pid + '.180.npy'
        shutil.copyfile(psi4_config["wfnfile"],targetfile)
        psi4.set_options({'guess':'read'})

    ## Perform computation
    G, wfn = psi4.gradient('mp2',**kwargs)
    wfn.to_file(psi4_config["wfnfile"])
    rdm = wfn.Da_subset('MO').to_array()
    NOONs = sorted(eigh(rdm)[0],reverse=True)

    num_alpha = wfn.nalpha()
    n_hono = NOONs[num_alpha-1]
    n_luno = NOONs[num_alpha]

    labels = [f'n_{i}' for i in range(len(NOONs))]
    str_NOONs = [str(NOON) for NOON in NOONs]

    print('n_hono,n_luno,' + ','.join(labels))
    print(f'{n_hono},{n_luno},' + ','.join(str_NOONs))

    os.chdir(basedir)
    return n_hono, n_luno

