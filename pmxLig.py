#!/usr/bin/env python3

import os
import sys
import zipfile
import time
import argparse
import subprocess

# gropy library to manipulate gro files
file_path = '/orozco/homes/adam/biobbdev/gropy/gropy-master/'
sys.path.append(os.path.dirname(file_path))
from gropy.Gro import Gro

# biobb common modules
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

# pycompss: biobb pmx modules
from biobb_adapters.pycompss.biobb_pmx.pmx.mutate_pc import mutate_pc
from biobb_adapters.pycompss.biobb_pmx.pmx.gentop_pc import gentop_pc
from biobb_adapters.pycompss.biobb_pmx.pmx.analyse_pc import analyse_pc

# pycompss: biobb md modules
from biobb_adapters.pycompss.biobb_md.gromacs.pdb2gmx_pc import pdb2gmx_pc
from biobb_adapters.pycompss.biobb_md.gromacs.make_ndx_pc import make_ndx_pc
from biobb_adapters.pycompss.biobb_md.gromacs.grompp_pc import grompp_pc
from biobb_adapters.pycompss.biobb_md.gromacs.grompp_cpt_pc import grompp_cpt_pc
from biobb_adapters.pycompss.biobb_md.gromacs.mdrun_pc import mdrun_pc
from biobb_adapters.pycompss.biobb_md.gromacs.mdrun_cpt_pc import mdrun_cpt_pc
from biobb_adapters.pycompss.biobb_md.gromacs_extra.append_ligand_pc import append_ligand_pc

# pycompss: biobb analysis modules
from biobb_adapters.pycompss.biobb_analysis.gromacs.gmx_image_pc import gmx_image_pc
from biobb_adapters.pycompss.biobb_analysis.gromacs.gmx_trjconv_str_ens_pc import gmx_trjconv_str_ens_pc


def main(config, system=None):
    start_time = time.time()

    
    conf = settings.ConfReader(config, system)
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    dhdl_paths_listA = []
    dhdl_paths_listB = []
    for ensemble, mutation in conf.properties['mutations'].items():

        ensemble_prop = conf.get_prop_dic(prefix=ensemble, global_log=global_log)
        ensemble_paths = conf.get_paths_dic(prefix=ensemble)

        # step0_image
        global_log.info(ensemble+" Step 0: gmx image: Imaging trajectories to remove PBC issues")
        ensemble_paths['step0_image']['input_traj_path'] = conf.properties['input_trajs'][ensemble]['input_traj_path']
        ensemble_paths['step0_image']['input_top_path'] = conf.properties['input_trajs'][ensemble]['input_tpr_path']
        gmx_image_pc(**ensemble_paths["step0_image"], properties=ensemble_prop["step0_image"])

        # step0_trjconv
        global_log.info(ensemble+" Step 0: gmx trjconv: Extract snapshots from equilibrium trajectories")
        ensemble_paths['step0_trjconv']['input_top_path'] = conf.properties['input_trajs'][ensemble]['input_tpr_path']
        gmx_trjconv_str_ens_pc(**ensemble_paths["step0_trjconv"], properties=ensemble_prop["step0_trjconv"])

        with zipfile.ZipFile(ensemble_paths["step0_trjconv"]["output_str_ens_path"]) as zip_f:
            zip_f.extractall()
            state_pdb_list = zip_f.namelist()

        for pdb_path in state_pdb_list:
            pdb_name = os.path.splitext(pdb_path)[0]
            prop = conf.get_prop_dic(prefix=os.path.join(ensemble, pdb_name), global_log=global_log)
            paths = conf.get_paths_dic(prefix=os.path.join(ensemble, pdb_name))

            # TODO: Check if this should be removed, replaced biobb_utils
            if(pdb_name == "frame98" or pdb_name == "frame12"):
                continue

            # step1_pmx_mutate
            global_log.info(ensemble+" Step 1: pmx mutate: Generate Hybrid Structure")
            paths['step1_pmx_mutate']['input_structure_path'] = pdb_path
            prop['step1_pmx_mutate']['mutation_list'] = mutation
            mutate_pc(**paths["step1_pmx_mutate"], properties=prop["step1_pmx_mutate"])

            # TODO: Check if this should be removed replaced biobb_utils
            if ensemble == 'stateA':
                mut = "L2R"
            elif ensemble == 'stateB':
                mut = "R2L"

            # TODO: Check if this should be removed,  replaced biobb_utils
            # Grep for dummy atoms
            cmd = "grep "+ mut +" wf_pmx/" + ensemble +"/"+pdb_name+"/step1_pmx_mutate/mut.gro  | cut -c12-15 |  sed 's/ //g' | sed -n '/^D/p'"
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            out, err = process.communicate()
            process.wait()
            dummy = False
            if len(out.decode("utf-8")) >= 2 :
                dummy = True
            print("CMD:" + cmd)
            print("DUMMY:\n" + str(out.decode("utf-8")))

            # TODO: Replace by biobb_utils
            # Remove ligand
            nolig_gro = mut +"_" + ensemble +"_"+pdb_name+".nolig.gro"
            system_gro = Gro()
            system_gro.read_gro_file(paths['step1_pmx_mutate']['output_structure_path'])
            system_gro.remove_residue_entry(292, 'IRE')
            system_gro.write_gro_file(nolig_gro)

            # step2_gmx_pdb2gmx
            paths['step2_gmx_pdb2gmx']['input_pdb_path'] = nolig_gro
            global_log.info(ensemble + " Step 2: gmx pdb2gmx: Generate Topology")
            pdb2gmx_pc(**paths["step2_gmx_pdb2gmx"], properties=prop["step2_gmx_pdb2gmx"])

            # TODO: Replace by biobb_utils
            # WARNING: It's using step1 input not step2
            # Re-ordering gro file, first ions, then water molecules
            ordered_gro = mut +"_" + ensemble +"_"+pdb_name+".sorted.gro"
            system_gro = Gro()
            system_gro.read_gro_file(paths['step1_pmx_mutate']['output_structure_path'])
            system_gro.sort_residues2(['NA', 'CL', 'SOL'])
            system_gro.renumber_atoms()
            system_gro.write_gro_file(ordered_gro)

            # step2_lig_gmx_appendLigand
            global_log.info(ensemble+" Step 2_lig: gmx appendLigand: Append a ligand to a GROMACS topology")
            append_ligand_pc(**paths["step2_lig_gmx_appendLigand"], properties=prop["step2_lig_gmx_appendLigand"])

            # step3_pmx_gentop
            global_log.info(ensemble+" Step 3: pmx gentop: Generate Hybrid Topology")
            gentop_pc(**paths["step3_pmx_gentop"], properties=prop["step3_pmx_gentop"])

            if not dummy:
                paths['step7_gmx_grompp']['input_gro_path'] = ordered_gro
            else:
                # step4_gmx_makendx
                global_log.info(ensemble+" Step 4 (Dummies): gmx make_ndx: Generate Gromacs Index file to select atoms to freeze")
                paths['step4_gmx_makendx']['input_structure_path'] = ordered_gro
                make_ndx_pc(**paths["step4_gmx_makendx"], properties=prop["step4_gmx_makendx"])

                # step5_gmx_grompp
                global_log.info(ensemble+" Step 5 (Dummies): gmx grompp: Creating portable binary run file for energy minimization")
                paths['step5_gmx_grompp']['input_gro_path'] = ordered_gro
                grompp_pc(**paths["step5_gmx_grompp"], properties=prop["step5_gmx_grompp"])

                # step6_gmx_mdrun
                global_log.info(ensemble+" Step 6 (Dummies): gmx mdrun: Running energy minimization")
                mdrun_pc(**paths["step6_gmx_mdrun"], properties=prop["step6_gmx_mdrun"])

            # step7_gmx_grompp
            global_log.info(ensemble+" Step 7: gmx grompp: Creating portable binary run file for system equilibration")
            grompp_pc(**paths["step7_gmx_grompp"], properties=prop["step7_gmx_grompp"])

            # step8_gmx_mdrun
            global_log.info(ensemble+" Step 8: gmx mdrun: Running system equilibration")
            mdrun_pc(**paths["step8_gmx_mdrun"], properties=prop["step8_gmx_mdrun"])

            # step9_gmx_grompp
            global_log.info(ensemble+" Step 9: Creating portable binary run file for thermodynamic integration (ti)")
            grompp_pc(**paths["step9_gmx_grompp"], properties=prop["step9_gmx_grompp"])

            # step10_gmx_mdrun
            global_log.info(ensemble+" Step 10: gmx mdrun: Running thermodynamic integration")
            mdrun_pc(**paths["step10_gmx_mdrun"], properties=prop["step10_gmx_mdrun"])

            if ensemble == "stateA":
                dhdl_paths_listA.append(paths["step10_gmx_mdrun"]["output_dhdl_path"])
            elif ensemble == "stateB":
                dhdl_paths_listB.append(paths["step10_gmx_mdrun"]["output_dhdl_path"])

    # Creating zip file containing all the dhdl files
    dhdlA_path = 'dhdlA.zip'
    dhdlB_path = 'dhdlB.zip'
    fu.zip_list(dhdlA_path, dhdl_paths_listA)
    fu.zip_list(dhdlB_path, dhdl_paths_listB)

    # step11_pmx_analyse
    global_log.info(ensemble+" Step 11: pmx analyse: Calculate free energies from fast growth thermodynamic integration simulations")
    global_paths["step11_pmx_analyse"]["input_A_xvg_zip_path"]=dhdlA_path
    global_paths["step11_pmx_analyse"]["input_B_xvg_zip_path"]=dhdlB_path
    analyse_pc(**global_paths["step11_pmx_analyse"], properties=global_prop["step11_pmx_analyse"])

    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % config)
    if system:
        global_log.info('  System: %s' % system)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Based on the official PMX tutorial")
    parser.add_argument('--config', required=True)
    parser.add_argument('--system', required=False)
    args = parser.parse_args()
    main(args.config, args.system)
