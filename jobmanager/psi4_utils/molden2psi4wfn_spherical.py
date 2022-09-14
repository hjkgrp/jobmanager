import numpy as np


def internal_ao_mapping(shell_type):
    if abs(shell_type) == 0:
        return {0: 0}
    elif abs(shell_type) == 1:
        return {0: 1, 1: 1, 2: -2}
    elif abs(shell_type) == 2:
        return {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    else:
        raise KeyError("only s, p, d shell is allowed.")


def shell_sequence_mapping(atom_shell_types):
    '''
    only work for def2-sv(p) now
    '''
    d = {}
    for ii, _ in enumerate(atom_shell_types):
        d.update({ii: 0})
    return d


def mocoeff_c2s(mocoeffs_c, shell_types):
    '''
    Covert Cartesian MO coefficients in TC to spherical MO coefficients.
    Work for the case where the internal calculation is done in spherical basis but the molden file is written in Cartesian basis.
    mocoeffs_c: MO coefficients at Cartesian basis (but with "sphericalbasis yes")
    shell_types: list of shell types ranked in the same order of MO coefficients (0 for s, 1 for p, and 2 for d)
    c2s: A matrix for conversion (for d shell)
    '''
    c2s = np.array([
        [-0.5*2/3, -0.5*2/3, 1*2/3, 0, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1],
        [np.sqrt(3)/2 * 2/3, -np.sqrt(3)/2*2/3, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0]
    ])
    mocoeffs_s = []
    start = 0
    end = 0
    for t in shell_types:
        if t == 0:
            end += 1
            mocoeffs_s += mocoeffs_c[start:end].tolist()
        elif t == 1:
            end += 3
            mocoeffs_s += mocoeffs_c[start:end].tolist()
        elif t == 2:
            end += 6
            tmp_c = mocoeffs_c[start:end]
            tmp_s = np.matmul(c2s, tmp_c)
            mocoeffs_s += tmp_s.tolist()
        else:
            raise KeyError("only support for s, p, and d for now.")
        start = end
    mocoeffs_s = np.array(mocoeffs_s)
    return mocoeffs_s


def mocoeff_mapping(d_molden, mapping, key='alpha'):
    coeff = np.copy(d_molden['orb_%s_coeffs' % key])
    coeff = mocoeff_c2s(coeff, d_molden['obasis']['shell_types'])
    mooeff_psi4 = np.zeros(shape=coeff.shape)
    for ii in range(len(mapping.keys())):
        #         print(ii, mapping[ii])
        mooeff_psi4[mapping[ii], :] = coeff[ii, ]
    return mooeff_psi4


def tcmolden2psi4wfn_ao_mapping_spherical(d_molden, restricted=False):
    mapping = {}
    start = 0
    ao_start = 0
    for atom_i in np.unique(d_molden['obasis']['shell_map']):
        num_shell = np.count_nonzero(d_molden['obasis']['shell_map'] == atom_i)
        atom_shell_types = d_molden['obasis']['shell_types'][start: start+num_shell]
        start += num_shell
        shell_mapping = shell_sequence_mapping(atom_shell_types)
#         print(shell_mapping)
        for ii in range(len(shell_mapping.keys())):
            shell_type = atom_shell_types[ii]
            ao_mapping_local = internal_ao_mapping(shell_type)
#             print(ii, shell_type, ao_mapping_local)
            for jj in range(len(ao_mapping_local.keys())):
                mapping.update(
                    {ao_start: ao_start+ao_mapping_local[jj]+shell_mapping[ii]})
                ao_start += 1
    Ca = mocoeff_mapping(d_molden, mapping, key='alpha')
    if not restricted:
        Cb = mocoeff_mapping(d_molden, mapping, key='beta')
    else:
        Cb = False
    return Ca, Cb, mapping
