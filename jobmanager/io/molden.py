import numpy as np
from typing import Dict, Tuple, List, Any
from jobmanager.constants import angstrom2bohr
from molSimplify.Classes.globalvars import amassdict
from warnings import warn


_shell_labels = ["s", "p", "d", "f", "sp", "g"]


def load_molden(filename: str) -> Dict[str, Any]:
    """Follows the specification defined in:
    https://www.theochem.ru.nl/molden/molden_format.html
    """
    with open(filename, "r") as fin:
        lines = fin.readlines()

    # Split the lines into "sections"
    sections: List[List[str]] = []  # List of all sections
    section: List[str] = []  # Current section
    for line in lines:
        if line.startswith("["):
            # Add current section to the list
            if section:  # Check if the section is empty
                sections.append(section)
            # Start new section
            section = []
        section.append(line.rstrip("\n"))
    # Append last section
    sections.append(section)
    # Parse the individual sections into a dictionary
    results: Dict[str, Any] = {}

    for section in sections:
        name = section[0].lower()
        if name.startswith("[title]"):
            results["title"] = "\n".join(section[1:])
        elif name.startswith("[atoms]"):
            length_unit = 1.0
            if "angs" in name:
                length_unit = angstrom2bohr
            (coordinates, numbers,
             pseudo_numbers) = read_geometry_section(section[1:])
            results["coordinates"] = coordinates * length_unit
            results["numbers"] = numbers
            results["pseudo_numbers"] = pseudo_numbers
        elif name.startswith("[gto]"):
            # Not added to the results dict right here because it also needs
            # the "centers" field that will come from the [atoms] section
            obasis = read_gto_section(section[1:])
        elif name.startswith("[mo]"):
            results.update(read_mo_section(section[1:]))
    # Add centers to obasis
    obasis["centers"] = results["coordinates"]
    results["obasis"] = obasis
    # I think the permutation field is only for internal use in iodata and
    # set to None here
    results["permutation"] = None

    if results["title"] != "Written by TeraChem":
        warn("This function has only been tested for TeraChem molden files."
             "Please double-check that your file was parsed correctly.")
    return results


def read_geometry_section(
        lines: List[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    coordinates = []
    numbers = []
    pseudo_numbers = []

    for line in lines:
        sp = line.split()
        if len(sp) != 6:
            break
        numbers.append(int(amassdict[sp[0]][1]))
        pseudo_numbers.append(int(sp[2]))
        coordinates.append([float(x) for x in sp[3:6]])
    # Convert to numpy arrays
    return np.array(coordinates), np.array(numbers), np.array(pseudo_numbers)


def read_gto_section(lines: List[str]) -> Dict[str, Any]:
    obasis: Dict[str, Any] = {}

    shell_map = []
    nprims = []
    shell_types = []
    alphas = []
    con_coeffs: List[float] = []  # Contraction coefficients

    line_iter = iter(lines)
    for line in line_iter:
        sp = line.split()
        if len(sp) == 0:  # Empty line
            continue
        elif len(sp) == 2 and sp[1] == "0":
            # Start of new list of shells, save current center
            center = int(sp[0]) - 1
        elif sp[0] in _shell_labels:
            # Start of shell, parse shell type and number of primitives
            shell_types.append(_shell_labels.index(sp[0]))
            nprims.append(int(sp[1]))
            shell_map.append(center)
            # Loop over the primitives and save alphas and con_coeffs
            for _ in range(nprims[-1]):
                line = next(line_iter)
                alpha, coeff = line.split()
                alphas.append(float(alpha))
                con_coeffs.append(float(coeff))

    obasis["shell_map"] = shell_map
    obasis["nprims"] = nprims
    obasis["shell_types"] = shell_types
    obasis["alphas"] = alphas
    obasis["con_coeffs"] = con_coeffs
    return obasis


def read_mo_section(lines: List[str]) -> Dict[str, Any]:
    orb_alpha_energies: List[float] = []
    orb_beta_energies: List[float] = []
    orb_alpha_occs: List[float] = []
    orb_beta_occs: List[float] = []
    orb_alpha_coeffs = []
    orb_beta_coeffs = []

    current_spin = "None"
    current_energy = 0.0
    current_occ = 0.0
    current_coeffs: List[float] = []
    for line in lines:
        sp = line.split()
        if len(sp) == 0:
            continue
        elif "=" in line:
            if current_coeffs:  # Not empty
                if current_spin.lower() == "alpha":
                    orb_alpha_energies.append(current_energy)
                    orb_alpha_occs.append(current_occ)
                    orb_alpha_coeffs.append(current_coeffs)
                elif current_spin.lower() == "beta":
                    orb_beta_energies.append(current_energy)
                    orb_beta_occs.append(current_occ)
                    orb_beta_coeffs.append(current_coeffs)
            if sp[0].lower() == "ene=":
                current_energy = float(sp[1])
            elif sp[0].lower() == "spin=":
                current_spin = sp[1].lower()
            elif sp[0].lower() == "occup=":
                current_occ = float(sp[1])
            # Ensure that the coeffs are reset
            current_coeffs = []
        else:
            current_coeffs.append(float(sp[1]))
    # Append last MO
    if current_spin.lower() == "alpha":
        orb_alpha_energies.append(current_energy)
        orb_alpha_occs.append(current_occ)
        orb_alpha_coeffs.append(current_coeffs)
    else:
        orb_beta_energies.append(current_energy)
        orb_beta_occs.append(current_occ)
        orb_beta_coeffs.append(current_coeffs)

    results = {
        "orb_alpha": np.shape(orb_alpha_coeffs),
        "orb_alpha_energies": np.array(orb_alpha_energies),
        "orb_alpha_occs": np.array(orb_alpha_occs),
        "orb_alpha_coeffs": np.array(orb_alpha_coeffs).T,
    }
    if orb_beta_energies:  # Not empty
        results["orb_beta"] = np.shape(orb_beta_coeffs)
        results["orb_beta_energies"] = np.array(orb_beta_energies)
        results["orb_beta_occs"] = np.array(orb_beta_occs)
        results["orb_beta_coeffs"] = np.array(orb_beta_coeffs).T
    else:
        # If the calculation is restricted ensure that the
        # occupation is only 1.0 instead of 2.0 for
        # compatibility with iodata
        results["orb_alpha_occs"] /= 2.0  # type: ignore
    return results
