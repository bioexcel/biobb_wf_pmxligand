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

# biobb pmx modules
from biobb_pmx.pmx.mutate import Mutate
from biobb_pmx.pmx.gentop import Gentop
from biobb_pmx.pmx.analyse import Analyse

# biobb md modules
from biobb_md.gromacs.pdb2gmx import Pdb2gmx
from biobb_md.gromacs.make_ndx import MakeNdx
from biobb_md.gromacs.grompp import Grompp
from biobb_md.gromacs.mdrun import Mdrun
from biobb_md.gromacs_extra.append_ligand import AppendLigand

# biobb analysis module
from biobb_analysis.gromacs.gmx_image import GMXImage
from biobb_analysis.gromacs.gmx_trjconv_str_ens import GMXTrjConvStrEns

def main(config, system=None):
    start_time = time.time()
    conf = settings.ConfReader(config, system)
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Added to avoid computing again what is already computed...
    # Get Current Working Directory
    #cwd = os.getcwd()

    dhdl_paths_listA = []
    dhdl_paths_listB = []
    for ensemble, mutation in conf.properties['mutations'].items():

        ensemble_prop = conf.get_prop_dic(prefix=ensemble, global_log=global_log)
        ensemble_paths = conf.get_paths_dic(prefix=ensemble)

        global_log.info(ensemble+" Step 0: gmx image: Imaging trajectories to remove PBC issues")
        ensemble_paths['step0_image']['input_traj_path'] = conf.properties['input_trajs'][ensemble]['input_traj_path']
        ensemble_paths['step0_image']['input_top_path'] = conf.properties['input_trajs'][ensemble]['input_tpr_path']
        GMXImage(**ensemble_paths["step0_image"], properties=ensemble_prop["step0_image"]).launch()

        global_log.info(ensemble+" Step 0: gmx trjconv: Extract snapshots from equilibrium trajectories")
        ensemble_paths['step0_trjconv']['input_top_path'] = conf.properties['input_trajs'][ensemble]['input_tpr_path']
        GMXTrjConvStrEns(**ensemble_paths["step0_trjconv"], properties=ensemble_prop["step0_trjconv"]).launch()

        with zipfile.ZipFile(ensemble_paths["step0_trjconv"]["output_str_ens_path"], 'r') as zip_f:
            zip_f.extractall()
            state_pdb_list = zip_f.namelist()


        for pdb_path in state_pdb_list:
            pdb_name = os.path.splitext(pdb_path)[0]
            prop = conf.get_prop_dic(prefix=os.path.join(ensemble, pdb_name), global_log=global_log)
            paths = conf.get_paths_dic(prefix=os.path.join(ensemble, pdb_name))

            # Added to avoid computing again what is already computed...
            #global_log.info("CHECKING: " + os.path.join(cwd,ensemble,pdb_name,'ti.dhdl.xvg'))
            #print (os.path.join(cwd,ensemble,pdb_name,'ti.dhdl.xvg'))
            #if os.path.isfile(os.path.join(cwd,ensemble,pdb_name,'ti.dhdl.xvg')):
            #    print("Frame " + pdb + " already computed. Jumping to next frame.")
            #else:

            #if(pdb_name == "frame19" or pdb_name == "frame39" or pdb_name == "frame73"):
            if(pdb_name == "frame98" or pdb_name == "frame12"):
                continue

            #Create and launch bb
            global_log.info(ensemble+" Step 1: pmx mutate: Generate Hybrid Structure")
            paths['step1_pmx_mutate']['input_structure_path'] = pdb_path
            prop['step1_pmx_mutate']['mutation_list'] = mutation
            Mutate(**paths["step1_pmx_mutate"], properties=prop["step1_pmx_mutate"]).launch()

            # Find out if generated gro file has any dummy atoms!! If so, run minimization (steps 5 and 6). If not, skip it!
            # PATH: wf_pmx/stateA/frame0/step1_pmx_mutate/mut.gro

            if ensemble == 'stateA':
                mut = "L2R"
            elif ensemble == 'stateB':
                mut = "R2L"

            cmd = "grep "+ mut +" wf_pmx/" + ensemble +"/"+pdb_name+"/step1_pmx_mutate/mut.gro  | cut -c12-15 |  sed 's/ //g' | sed -n '/^D/p'"
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            out, err = process.communicate()
            process.wait()
            dummy = False
            if len(out.decode("utf-8")) >= 2 :
                dummy = True
            #dummy = out.decode("utf-8")

            print("CMD:" + cmd)
            print("DUMMY:\n" + str(out.decode("utf-8")))

            # Step 2: gmx pdb2gmx: Generate Topology
            # From pmx tutorial:
            # gmx pdb2gmx -f mut.pdb -ff amber99sb-star-ildn-mut -water tip3p -o pdb2gmx.pdb
            global_log.info(ensemble+" Step 2: gmx pdb2gmx: Generate Topology")

            # First of all, remove ligand from the GRO structure
            #cmd = "grep -v AQ4 wf_pmx/" + ensemble +"/"+pdb_name+"/step1_pmx_mutate/mut.gro > structure.nolig.gro"
            #process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            #out, err = process.communicate()
            #process.wait()
            #print (out.decode("utf-8"))

            # First of all, remove ligand from the GRO structure
            # Using library for gro files manipulation (gropy)
            # https://github.com/caizkun/gropy
            nolig_gro = mut +"_" + ensemble +"_"+pdb_name+".nolig.gro"
            system_gro = Gro()
            system_gro.read_gro_file(paths['step1_pmx_mutate']['output_structure_path'])
            system_gro.remove_residue_entry(292,'IRE')
            #system_gro.renumber_atoms()
            system_gro.write_gro_file(nolig_gro)

            # Now use this apo structure to generate the GROMACS topology
            paths['step2_gmx_pdb2gmx']['input_pdb_path'] = nolig_gro
            Pdb2gmx(**paths["step2_gmx_pdb2gmx"], properties=prop["step2_gmx_pdb2gmx"]).launch()

            # Re-ordering gro file, first ions, then water molecules
            # Using library for gro files manipulation (gropy)
            # https://github.com/caizkun/gropy
            ordered_gro = mut +"_" + ensemble +"_"+pdb_name+".sorted.gro"
            system_gro = Gro()
            #system_gro.read_gro_file(paths['step2_gmx_pdb2gmx']['output_gro_path'])
            system_gro.read_gro_file(paths['step1_pmx_mutate']['output_structure_path'])
            system_gro.sort_residues2(['NA', 'CL', 'SOL'])
            system_gro.renumber_atoms()
            system_gro.write_gro_file(ordered_gro)

            # And now add the ligand to the previously generated topology
            # Step 2_lig: gmx appendLigand: Append a ligand to a GROMACS topology
            global_log.info(ensemble+" Step 2_lig: gmx appendLigand: Append a ligand to a GROMACS topology")
            AppendLigand(**paths["step2_lig_gmx_appendLigand"], properties=prop["step2_lig_gmx_appendLigand"]).launch()

            # Step 3: pmx gentop: Generate Hybrid Topology
            # From pmx tutorial:
            # python generate_hybrid_topology.py -itp topol_Protein.itp -o topol_Protein.itp -ff amber99sb-star-ildn-mut
            global_log.info(ensemble+" Step 3: pmx gentop: Generate Hybrid Topology")
            Gentop(**paths["step3_pmx_gentop"], properties=prop["step3_pmx_gentop"]).launch()

            #if ensemble == 'stateA':
            if not dummy:

                # In stateA (lamdda=0), if the residue to be mutated is smaller than the new one, no dummy atom is generated in the hybrid topology.
                # This means we don't need the energy minimization step, so simply get the output
                # from the step2 (pdb2gmx) as output from the step6 (energy minimization)
                # From pmx tutorial:
                # There are no dummies in this state at lambda=0, therefore simply convert mut.pdb to emout.gro
                #paths['step7_gmx_grompp']['input_gro_path'] = paths['step2_gmx_pdb2gmx']['output_gro_path']
                paths['step7_gmx_grompp']['input_gro_path'] = ordered_gro

            #elif ensemble == 'stateB':
            else:

                # Step 4 (Dummies): gmx make_ndx: Generate Gromacs Index File to select atoms to freeze
                # From pmx tutorial:
                # echo -e "a D*\n0 & ! 19\nname 20 FREEZE\nq\n" | gmx make_ndx -f frame0/pdb2gmx.pdb -o index.ndx
                global_log.info(ensemble+" Step 4 (Dummies): gmx make_ndx: Generate Gromacs Index file to select atoms to freeze")
                paths['step4_gmx_makendx']['input_structure_path'] = ordered_gro
                MakeNdx(**paths["step4_gmx_makendx"], properties=prop["step4_gmx_makendx"]).launch()

                # Step 5 (Dummies): gmx grompp: Creating portable binary run file for energy minimization
                # From pmx tutorial:
                # gmx grompp -c pdb2gmx.pdb -p topol.top -f ../../mdp/em_FREEZE.mdp -o em.tpr -n ../index.ndx
                global_log.info(ensemble+" Step 5 (Dummies): gmx grompp: Creating portable binary run file for energy minimization")
                paths['step5_gmx_grompp']['input_gro_path'] = ordered_gro
                Grompp(**paths["step5_gmx_grompp"], properties=prop["step5_gmx_grompp"]).launch()

                # Step 6: gmx mdrun: Running energy minimization
                # From pmx tutorial:
                # gmx mdrun -s em.tpr -c emout.gro -v
                global_log.info(ensemble+" Step 6 (Dummies): gmx mdrun: Running energy minimization")
                Mdrun(**paths["step6_gmx_mdrun"], properties=prop["step6_gmx_mdrun"]).launch()

                #paths['step7_gmx_grompp']['input_ndx_path'] = paths['step4_gmx_makendx']['output_ndx_path']

            # Step 7: gmx grompp: Creating portable binary run file for system equilibration
            # From pmx tutorial:
            # gmx grompp -c emout.gro -p topol.top -f ../../mdp/eq_20ps.mdp -o eq_20ps.tpr -maxwarn 1
            global_log.info(ensemble+" Step 7: gmx grompp: Creating portable binary run file for system equilibration")
            Grompp(**paths["step7_gmx_grompp"], properties=prop["step7_gmx_grompp"]).launch()

            # Step 8: gmx mdrun: Running system equilibration
            # From pmx tutorial:
            # gmx mdrun -s eq_20ps.tpr -c eqout.gro -v
            global_log.info(ensemble+" Step 8: gmx mdrun: Running system equilibration")
            Mdrun(**paths["step8_gmx_mdrun"], properties=prop["step8_gmx_mdrun"]).launch()

            # Step 9: gmx grompp: Creating portable binary run file for thermodynamic integration (ti)
            # From pmx tutorial:
            # gmx grompp -c eqout.gro -p topol.top -f ../../mdp/ti.mdp -o ti.tpr -maxwarn 1
            global_log.info(ensemble+" Step 9: Creating portable binary run file for thermodynamic integration (ti)")
            Grompp(**paths["step9_gmx_grompp"], properties=prop["step9_gmx_grompp"]).launch()

            # Step 10: gmx mdrun: Running thermodynamic integration
            # From pmx tutorial:
            # gmx mdrun -s ti.tpr -c eqout.gro -v
            global_log.info(ensemble+" Step 10: gmx mdrun: Running thermodynamic integration")
            Mdrun(**paths["step10_gmx_mdrun"], properties=prop["step10_gmx_mdrun"]).launch()
            if ensemble == "stateA":
                dhdl_paths_listA.append(paths["step10_gmx_mdrun"]["output_dhdl_path"])
            elif ensemble == "stateB":
                dhdl_paths_listB.append(paths["step10_gmx_mdrun"]["output_dhdl_path"])

    #Creating zip file containing all the dhdl files
    dhdlA_path = 'dhdlA.zip'
    dhdlB_path = 'dhdlB.zip'
    fu.zip_list(dhdlA_path, dhdl_paths_listA)
    fu.zip_list(dhdlB_path, dhdl_paths_listB)

    # Step 11: pmx analyse: Calculate free energies from fast growth thermodynamic integration simulations
    # From pmx tutorial:
    # python analyze_dhdl.py -fA ../stateA/frame*/dhdl*.xvg -fB ../stateB/frame*/dhdl*.xvg --nbins 25 -t 293 --reverseB
    global_log.info(ensemble+" Step 11: pmx analyse: Calculate free energies from fast growth thermodynamic integration simulations")
    global_paths["step11_pmx_analyse"]["input_A_xvg_zip_path"]=dhdlA_path
    global_paths["step11_pmx_analyse"]["input_B_xvg_zip_path"]=dhdlB_path
    Analyse(**global_paths["step11_pmx_analyse"], properties=global_prop["step11_pmx_analyse"]).launch()

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
