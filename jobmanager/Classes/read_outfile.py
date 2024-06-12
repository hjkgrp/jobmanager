from itertools import islice
import json

class read_outfile():
    def __init__(self,fname,jsonfile = None):
        self.fname = fname
        if jsonfile is None:
            self.current_scf = 0
            self.linenum = 0
            self.currently_running = False
            self.energies = []
        else:
            with open(jsonfile, 'r') as f:
                temp = json.load(f)
            self.current_scf = temp['current_scf']
            self.linenum = temp['linenum']
            self.currently_running = temp['currently_running']
            self.energies = temp['energies']

    def read_from(self, end_line=None):
        if self.current_scf == 0:
            start = 0
            starting_point = False
        else:
            start = self.linenum+1
            starting_point = self.currently_running
            if starting_point:
                energy_this_scf = self.energies[0]
        reached_final_energy = False
        with open(self.fname) as lines:
            iter = islice(enumerate(lines),start,end_line,1)
            for linenum, line in iter:
                line_list_no_space = line.split()
                if "Start SCF Iterations" in line:
                    self.current_scf+=1
                    starting_point = True
                    energy_this_scf = []
                elif starting_point and len(line_list_no_space) == 11 and line_list_no_space[0].isdigit():
                    energy_this_scf.append(float(line_list_no_space[-2]))
                elif 'FINAL ENERGY:' in line and starting_point:
                    reached_final_energy = True
                    starting_point = False
                    final_energy = float(line_list_no_space[2])
                    # self.energies.append(final_energy)
                self.linenum = linenum
            if self.current_scf != 0 and starting_point:
                self.energies = [energy_this_scf]
                self.currently_running = True
            elif self.current_scf != 0 and reached_final_energy:
                self.energies = final_energy
                self.currently_running = False
            else:
                pass

        dict = {'linenum' : self.linenum, 'current_scf':self.current_scf,
                'currently_running':self.currently_running,'energies':self.energies}
        with open(self.fname.rsplit('.',1)[0]+'.json', 'w') as json_file:
            json.dump(dict, json_file)
        return dict
