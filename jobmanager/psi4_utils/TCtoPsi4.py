import numpy as np

class TCtoPsi4:
    """
    Contains functions useful for converting the Cartesian molden files generated by TeraChem
    to Psi4 .wfn files.

    Note that these conversions (currently) only work for 6-31g* or LACVP* for Cartesian, def2 for Spherical.
    """

    def __init__(self, basis, coordinates):
        """
        Parameters:
            basis: str
                Basis set used for the calculations.
                Must be the same between TeraChem and Psi4, currently only support
                LACVP*, 6-31g(*) and def2 basis sets (lacvps, 6-31g, 6-31g*, and def2-sv(p) are supported inputs)
            coordinates: str
                Coordinates that will be used for the Psi4 calculation.
                Currently only support Cartesian for LACVP* and Spherical for def2 basis sets.
                (use Cartesian for puream=False, Spherical for puream=True)
        """
        assert coordinates in ['Cartesian', 'Spherical']
        supported_basis = ['lacvps', 'def2', '6-31g']
        assert any([supported in basis for supported in supported_basis])
        if 'def2' in basis:
            self.def2 = True
            self.lacvps = False
        elif 'lacvps' in basis or '6-31g' in basis:
            #lacvps and 6-31g basis sets can be handled analogously
            self.def2 = False
            self.lacvps = True
        if coordinates == 'Cartesian':
            self.cartesian = True
            self.spherical = False
        elif coordinates == 'Spherical':
            self.cartesian = False
            self.spherical = True
        assert (self.cartesian and self.lacvps) or (self.spherical and self.def2)

    def internal_ao_mapping(self, shell_type):
        """
        Informs how to shift the individual AOs within a subshell (same n, l)
        Read as i:j where i is the position of the original (TC) AO,
        and j is the position of the AO for Psi4.
        e.g., for Cartesian p (shell_type==1), the px, py, and pz orderings are the same,
        so a zero shift for 3 of the AOs.
        Note: each subshell is not a full subshell, but rather a basis function for that subshell.
        So, there could be multiple calls of this function for one subshell for multiple-zeta.

        Parameters:
            shell_type: int
                0 is s, 1 is p, 2 is d.

        Outputs:
            dict
                maps TC AO indexes to Psi4 AO indexes within a subshell.
        """
        if self.cartesian:
            if abs(shell_type) == 0:
                return {0: 0}
            elif abs(shell_type) == 1:
                return {0: 0, 1: 0, 2: 0}
            elif abs(shell_type) == 2:
                return {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
            else:
                raise KeyError("only s, p, d shell is allowed.")
        if self.spherical:
            if abs(shell_type) == 0:
                return {0: 0}
            elif abs(shell_type) == 1:
                return {0: 1, 1: 1, 2: -2}
            elif abs(shell_type) == 2:
                return {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
            else:
                raise KeyError("only s, p, d shell is allowed.")


    def shell_sequence_mapping(self, atom_shell_types):
        """
        Informs how to shift the subshells (groups of AOs with same n, l)
        Read as i:j where i is the position of the original (TC) set of AOs
        and j is the position of the set of AOs in Psi4
        Note: tells how to shift the AOs in the subshell, not how the subshells related.
        The number of elements in the dictionary is not the number of subshells, but rather the
        number of basis functions used to represent the atom.
        See below for examples. The notation 1s(6) means the basis function for the 1s subshell
        composed of 6 primitives. So, for valence subshells, have 2 basis functions per subshell,
        one with 3 primitives and the other with 1 (in 6-31g*).

        For Cartesian, this function has only been designed for LACVP* (6-31g* for H - Ar, LANL2DZ for heavier atoms)
        For Spherical, this function has only been designed for def2-SV(P), where the orderings are the same

        Parameters:
            atom_shell_types: list
                List of the shell types associated with a particular atom.
        Outputs:
            dict
                maps subshell indices in TC to those in Psi4.
        """
        if self.cartesian and self.lacvps:
            # 6-31g* on C, O, N...
            #Molden is 1s(6), 2s(3), 2s(1), 2p(3), 2p(1), d
            #Generates 1,     1,     1,     3,     3,     6 orbitals
            #Psi4 is   1s(6), 2s(3), 2p(3), 2s(1), 2p(1), d
            #Generates 1,     1,     3,     1,     3,     6 orbitals
            #So, need to shift 2s(1) orbitals right 3, 2p(3) orbitals left one.
            if "-".join([str(x) for x in atom_shell_types]) == "-".join([str(x) for x in [0, 0, 0, 1, 1, 2]]):
                return {0: 0, 1: 0, 2: 3, 3: -1, 4: 0, 5: 0}  # for TC

            # 6-31g on C, O, N...
            #Same as above but with no d function.
            elif "-".join([str(x) for x in atom_shell_types]) == "-".join([str(x) for x in [0, 0, 0, 1, 1]]):
                return {0: 0, 1: 0, 2: 3, 3: -1, 4: 0}

            # 6-31g* on P, S, Cl...
            #Molden is 1s(6), 2s(6), 3s(3), 3s(1), 2p(6), 3p(3), 3p(1), d
            #Generates 1,     1,     1,     1,     3,     3,     3,     6 orbitals
            #Psi4   is 1s(6), 2s(6), 2p(6), 3s(3), 3p(3), 3s(1), 3p(1), d
            #Generates 1,     1,     3,     1,     3,     1,     3,     6 orbitals
            #So, need to shift 3s(3) orbitals right 3, 3s(1) right 6, 2p(6) left 2, 3p(3) left 1.
            elif "-".join([str(x) for x in atom_shell_types]) == "-".join([str(x) for x in [0, 0, 0, 0, 1, 1, 1, 2]]):
                return {0: 0, 1: 0, 2: 3, 3: 6, 4: -2, 5: -1, 6: 0, 7: 0}

            # 6-31g* on metals
            elif "-".join([str(x) for x in atom_shell_types]) == "-".join([str(x) for x in [0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2]]):
                return {0: 0, 1: 0, 2: 3, 3: 6, 4: 9, 5: -3, 6: -2, 7: -1, 8: 0, 9: 0, 10: 0}

            else:
                #In LACVP*, the LANL2DZ ordering is the same in TeraChem and Psi4
                #so just return 0 shift for all of the non-(6-31g) atoms
                d = {}
                for ii, _ in enumerate(atom_shell_types):
                    d.update({ii: 0})
                return d
        elif self.spherical and self.def2:
            #In def2-SV(P), the ordering is the same between TeraChem and Psi4
            d = {}
            for ii, _ in enumerate(atom_shell_types):
                d.update({ii: 0})
            return d


    def mocoeff_c2s(self, mocoeffs_c, shell_types):
        '''
        Psi4 normalizes basis functions according to (for d-orbitals):
        For Cartesian functions we need to add a basis function normalization constant of
        //      _______________________________
        //     / (2lx-1)!! (2ly-1)!! (2lz-1)!!
        //    /  -----------------------------
        //  \/             (2l-1)!!
        //

        d_conv, d_conv_inv: conversion metrix.
        source: https://github.com/psi4/psi4/blob/master/psi4/src/psi4/libmints/writer.cc, FCHKWriter

        For Cartesian: converts the molecular orbital coefficients
        (as represented by number of basis functions, where each s is counted once, p thrice, d six times)
        from TeraChem normalization to Psi4 normalization, by the above scale factor

        For Spherical: converts the Cartesian TeraChem MO coefficients into a spherical basis

        Parameters:
            moeffs_c: array
                Gives the cartesian MO coefficients from the TC molden.
            shell_types: array
                Gives the shell types associated with a particular atom.

        Outputs:
            mocoeffs_s: array
                MO coefficients to use in Psi4 (converted from their TC values); not reordered
        '''  # noqa W605
        if self.cartesian:
            #to normalize as above for the d-orbitals
            pf1 = 1.0
            pf2 = np.sqrt(1.0 / 3.0)
            d_conv = np.array([
                [pf1, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, pf1, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, pf1],
                [0.0, pf2, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, pf2, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, pf2, 0.0],
            ])
            tc2psi4 = np.linalg.inv(d_conv)
        if self.spherical:
            #to convert the cartesian basis to a spherical one
            tc2psi4 = np.array([
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
                tmp_s = np.matmul(tc2psi4, tmp_c)
                mocoeffs_s += tmp_s.tolist()
            else:
                raise KeyError("only support for s, p, and d for now.")
            start = end
        mocoeffs_s = np.array(mocoeffs_s)
        return mocoeffs_s


    def mocoeff_mapping(self, d_molden, mapping, key='alpha'):
        """
        Maps the molecular orbital coefficients (again represented with number of basis functions
        with s counted once, p thrice, and d six times) in TeraChem to that in Psi4, using the
        mapping developed in mapping.

        Parameters:
            d_molden: dict
                Dictionary containing information from the molden file, loaded from io.load_molden()
            mapping: dict
                Maps the MO coefficient indexes in TC format (keys) to those in Psi4 (values)
            key: str
                Use alpha/beta, denotes which coefficients should be copied from the molden.

        Outputs:
            moeff_psi4: array
                MO coefficients to use in Psi4 (after reordering)
        """
        coeff = np.copy(d_molden['orb_%s_coeffs' % key])
        coeff = self.mocoeff_c2s(coeff, d_molden['obasis']['shell_types'])
        mooeff_psi4 = np.zeros(shape=coeff.shape)
        for ii in range(len(mapping.keys())):
            #         print(ii, mapping[ii])
            mooeff_psi4[mapping[ii], :] = coeff[ii, :]
        return mooeff_psi4


    def tcmolden2psi4wfn_ao_mapping(self, d_molden, restricted=False):
        """
        mapping returns a dictionary mapping the MO coefficients in the TC molden to those in a Psi4 wfn.
        Ca and Cb are the coefficients for the alpha and beta orbitals, again read from the molden.

        Combines the above functions to get the alpha and beta coefficient matrix in Psi4 from the TC molden.

        Parameters:
            d_molden: dict
                Dictionary containing information from the molden file, loaded from io.load_molden()
            restricted: bool
                Whether or not the calculation was / will be run restricted or not. If True,
                the beta coefficients are identical to the alpha coefficients.

        Outputs:
            Ca, Cb: array
                MO coefficients for the alpha and beta electrons, respectively.
            mapping: dict
                Maps the MO coefficient indices from those in TC to those in Psi4.
        """
        mapping = {}
        start = 0 #iterate over the number of basis functions
        ao_start = 0 #iterate over the atomic orbitals (1 for s basis, 3 for p, 6 for d)
        for atom_i in np.unique(d_molden['obasis']['shell_map']): #for each atom in the molecule
            num_shell = np.count_nonzero(d_molden['obasis']['shell_map'] == atom_i) #get the number of basis functions associated with the atom
            atom_shell_types = d_molden['obasis']['shell_types'][start: start+num_shell] #get the identity of associated shells
            start += num_shell
            shell_mapping = self.shell_sequence_mapping(atom_shell_types) #get a dictionary mapping how shells are ordered for the atom
            #print(shell_mapping)
            for ii in range(len(shell_mapping.keys())): #for each basis function
                shell_type = atom_shell_types[ii] #get the shell type
                ao_mapping_local = self.internal_ao_mapping(shell_type) #reorder AOs in the shell as needed
                #print(ii, shell_type, ao_mapping_local)
                for jj in range(len(ao_mapping_local.keys())): #for each AO
                    #map the MO position to what it would be in Psi4
                    #adjust the MO by the rearrangement in the basis function, and by rearrangements in basis functions for each atom
                    mapping.update(
                        {ao_start: ao_start+ao_mapping_local[jj]+shell_mapping[ii]})
                    ao_start += 1
        #assign the coefficients
        Ca = self.mocoeff_mapping(d_molden, mapping, key='alpha')
        if not restricted:
            Cb = self.mocoeff_mapping(d_molden, mapping, key='beta')
        else:
            Cb = False
        return Ca, Cb, mapping
