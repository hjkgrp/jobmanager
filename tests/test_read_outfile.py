from jobmanager.read_outfile import read_outfile

def test_read_outfile(resource_path_root):
    file_1 = str(resource_path_root/'inputs'/'read_outfile'/'first.out') # full 3 scf cycles
    file_2 = str(resource_path_root/'inputs'/'read_outfile'/'second.out') # starts 4th scf but has only 'Start SCF Iterations'
    file_3 = str(resource_path_root/'inputs'/'read_outfile'/'third.out') # 0-9 energies in 4th scf (not full)
    file_4 = str(resource_path_root/'inputs'/'read_outfile'/'fourth.out') # 0-15 energies in 4th scf (not full)
    file_5 = str(resource_path_root/'inputs'/'read_outfile'/'fifth.out') # complete 4th scf cycle with 'FINAL ENERGY'

    sample_1 = read_outfile(file_1)
    sample_1.read_from()
    bool_1 = (sample_1.energies == -602.4607426882)

    sample_2 = read_outfile(file_2, jsonfile=file_1.rsplit('.',1)[0]+'.json')
    sample_2.read_from()
    bool_2 = (sample_2.energies == [[]])

    sample_3 = read_outfile(file_3, jsonfile=file_2.rsplit('.',1)[0]+'.json')
    sample_3.read_from()
    bool_3 = (sample_3.energies == [[-602.456973586, -602.2724649804, -602.4452898263,
                                    -602.4615364022, -602.4617684258, -602.4618023219,
                                    -602.461810432, -602.4618403048, -602.4618436049]])

    sample_4 = read_outfile(file_4, jsonfile=file_3.rsplit('.',1)[0]+'.json')
    sample_4.read_from()
    bool_4 = (sample_4.energies == [[-602.456973586, -602.2724649804, -602.4452898263,
                                     -602.4615364022, -602.4617684258, -602.4618023219,
                                     -602.461810432, -602.4618403048, -602.4618436049,
                                     -602.4618444624, -602.4618456718, -602.4618466261,
                                     -602.4618471806, -602.4618473329, -602.4618475092]])

    sample_5 = read_outfile(file_5, jsonfile=file_4.rsplit('.',1)[0]+'.json')
    sample_5.read_from()
    bool_5 = (sample_5.energies == -602.4618478648)

    assert (bool_1 and bool_2 and bool_3 and bool_4 and bool_5)
