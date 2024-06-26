import os

class RunUtils():
    """
    Contains functions useful for running Psi4 calculations.
    """

    def __init__(self):
        pass

    def ensure_dir(self, dirpath):
        #Checks if a directory exists. If not, makes that directory.
        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)

    def check_sucess(self, path='./', output_name='output.dat'):
        """
        Checks if a Psi4 calculation is completed.
        Assumes it is being called in the same directory as the Psi4 output file,
        and that the output file is named output.dat.
        """
        success = False
        with open(path + output_name, "r") as fo:
            txt = "".join(fo.readlines())
        #Psi4 will write Computation Completed once the calculation finishes.
        if 'Computation Completed' in txt:
            success = True
        return success