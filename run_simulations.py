import os
import subprocess
import multiprocessing
import shutil
import MDAnalysis as md
import h5py
import numpy as np
import math

# A function to log everything we do
def log_message(message, log_file="simulation_log.txt"):
    with open(log_file, "a") as f:
        f.write(message + "\n")
    print(message)

# This function makes sure that the exp.mdp file has all the right parameters we want
def change_mdp(file_name, energy, time, peak, width, gen_log):
    with open(file_name, 'r') as file:
        file_data = file.read()
    file_data = file_data.replace('XXX', str(time))
    file_data = file_data.replace('600;', f"{str(energy)};")
    file_data = file_data.replace('0.005;', f"{str(peak)};")
    file_data = file_data.replace('0.00125;', f"{str(width)};")
    if gen_log == 'yes':
        file_data = file_data.replace('0; Enable logging of electronic dynamics.', '1; Enable logging of electronic dynamics.')
    with open(file_name, 'w') as file:
        file.write(file_data)

# This giant function is what actually initiates the simulation
def sim_run(parameters):

    # First, we load all of the input parameters here:
    num_cores, gro, cwd_project, cwd_structures, conformation_name, time_steps, photon_energy, pulse_peak, pulse_width, simulation_number, gen_log = parameters

    #We make a folder for this simulation
    cwd_simulation = f"{cwd_project}/sim{simulation_number}"
    os.makedirs(cwd_simulation) 
    os.chdir(cwd_simulation)
    os.makedirs(f"simulation_output")       # This is an import file that is needed for the simulation to work

    # Here we get an exp.mdp file copied to our simulation folder
    shutil.copy('/home/simon/mds_params/exp.mdp', '.')

    # Here we define some paths to important files we will need
    top = f"{cwd_structures}/tops/{conformation_name}.top"
    pdb = f"{cwd_structures}/pdbs/{conformation_name}.pdb"
    tpr = f"{conformation_name}.tpr"
    trr = f"{conformation_name}.trr"
    exp_pdb = f"exp.pdb"
    edr = f"ener.edr"
    xtc = f"{conformation_name}.xtc"
    h5file = f"data.h5"  # This is what the h5-file will be named
    exp_mdp = f"exp.mdp"

    # Here are some paths to important functions
    moldstruct_path = "/home/spidocstester/MolDStruct/bin"  # Path to MolDStruct executables
    grompp = f"{moldstruct_path}/grompp"  # Path to grompp executable
    mdrun = f"{moldstruct_path}/mdrun"  # Path to mdrun executable
    g_energy = f"{moldstruct_path}/g_energy" 

    # Now we change the mdp-file so that it includes the right parameters
    change_mdp(exp_mdp, photon_energy, time_steps, pulse_peak, pulse_width, gen_log)

    # Now we copy the important files we need to intialize the atomic parameters ...
    shutil.copytree('/home/simon/mds_params/Atomic_model', f"Atomic_model/")
    shutil.copy('/home/simon/mds_params/generate_atomic_parameters.py', '.')
    # And then we actually calculate them
    os.system(f"python3 generate_atomic_parameters.py {photon_energy} {pdb} . .")

    # Run grompp command to prepare the simulation input
    grompp_cmd = [grompp, '-f', exp_mdp, '-c', gro, '-p', top, '-o', tpr, '-maxwarn', '5']
    log_message(f"Running grompp command: {' '.join(grompp_cmd)}")
    try:
        result = subprocess.run(grompp_cmd, check=True, capture_output=True, text=True)
        log_message(f"grompp command output: {result.stdout}")
        log_message("grompp command completed successfully")
    except subprocess.CalledProcessError as e:
        log_message(f"Error in grompp command: {e}")
        log_message(f"Error output: {e.stderr}")

    # Run mdrun command to execute the simulation
    mdrun_cmd = [mdrun, '-s', tpr, '-o', trr, '-x', xtc, '-c', exp_pdb, '-v', '-nt', '1', '-ionize']
    log_message(f"Running mdrun command: {' '.join(mdrun_cmd)}")
    try:
        result = subprocess.run(mdrun_cmd, check=True, capture_output=True, text=True)
        log_message(f"mdrun command output: {result.stdout}")
        log_message("mdrun command completed successfully")
    except subprocess.CalledProcessError as e:
        log_message(f"Error in mdrun command: {e}")
        log_message(f"Error output: {e.stderr}")
    
    # Run energy command to get the energy data
    for i in [['pot.xvg', '9'], ['kin.xvg', '10'], ['tot.xvg', '11']]:
        energy_cmd = [g_energy, '-f', edr, '-o', f"{i[0]}"]
        log_message(f"Running energy command: {' '.join(mdrun_cmd)}")
        try:
            result = subprocess.run(energy_cmd, input = f"{i[1]}\n", check=True, capture_output=True, text=True)
            os.system('10')
            os.system(' ')
            log_message(f"mdrun command output: {result.stdout}")
            log_message("mdrun command completed successfully")
        except subprocess.CalledProcessError as e:
            log_message(f"Error in energy command: {e}")
            log_message(f"Error output: {e.stderr}")

    # Detector implementation
    try:
        universe = md.Universe(gro, trr)  # Load the trajectory data using MDAnalysis
        ag = universe.atoms.select_atoms("all")  # Select all atoms in the trajectory
    except Exception as e:
        log_message(f"Error loading trajectory: {e}")
        return
    
    # Collect mass data for the atoms (only once per simulation set)
    if simulation_number == 1:
        mass_data = ag.masses.tolist()  # Get the mass data of all atoms
        mass_data = np.array(mass_data)
        log_message(f"Collected mass data for atoms: {mass_data}")
    
    universe.trajectory[0]  # Move to the first frame of the trajectory
    vel_i = ag.velocities.copy()  # Copy initial velocities
    pos_i = ag.positions.copy()  # Copy initial positions
    universe.trajectory[-1]  # Move to the last frame of the trajectory
    vel_f = ag.velocities.copy()  # Copy final velocities
    pos_f = ag.positions.copy()  # Copy final positions
    vel_data = [(x / np.linalg.norm(x)) if np.linalg.norm(x) != 0 else x for x in vel_f]  # Normalize velocities
    pos_data = [(x / np.linalg.norm(x)) if np.linalg.norm(x) != 0 else x for x in (pos_f - pos_i)]  # Normalize displacements

  # Save data to the HDF5 file
    os.chdir('..')
    with h5py.File(h5file, 'a') as file:
        group_path = f"sim{simulation_number}"
        group = file.require_group(group_path)  # Create or get the group for this simulation

        # Save mass data (only once per set of simulations)
        if simulation_number == 1:
            group.create_dataset("mass", data=mass_data)
            log_message(f"Mass data saved to {group_path}")
        # Save other simulation data
        group.create_dataset("unit_velocity", data=vel_data)
        group.create_dataset("unit_displacement", data=pos_data)
        group.create_dataset("initial_position", data=pos_i)
        group.create_dataset("final_position", data=pos_f)
        group.create_dataset("initial_velocity", data=vel_i)
        group.create_dataset("final_velocity", data=vel_f)
        log_message(f"Simulation data saved to {group_path}")




if __name__ == "__main__":
    try:   
        # What should the project name be?
        project_name = 'ubi48_quick_pulse_unaligned_2000ev'

        # Type here what the name of the conformation is (i.e. ubiWT, ubi6, ...)
        conformation_name = 'ubi48'

        # Fill in these simulation parameters accordingly
        number_of_simulation = 100
        time_steps = 200000
        photon_energy = 2000
        pulse_peak = 0.005
        pulse_width = 5 / (2 * math.sqrt(2 * math.log(2))) * 0.001

        # How many cores should we use?
        num_cores = 2

        # Should the results be loged? ('yes' or 'no')
        gen_log = 'yes'

        # Are we going to be using the gros at different equilibrium time steps? ('yes' or 'no')
        use_time_steps = 'yes'



        ########################################################################
        #### The rest of the things here are not things you should input :) ####
        ########################################################################

        # Create and switch to the project folder
        cwd_project = f"/home/simon/results/{project_name}"
        os.makedirs(cwd_project)
        os.chdir(cwd_project)

        # This is where all the necessary structure files are located
        cwd_structures = '/home/simon/structure_files'



        # If we want to use the equilibrated gros, then we append all of the gros to a list and input that as the gro parameter
        if use_time_steps == 'yes':
            t_end = 250
            t_start = 50
            gros = []
            for t in np.linspace(t_start, t_end, number_of_simulation):
                gro = f"{cwd_structures}/time_step_gros/{conformation_name}/{conformation_name}_{int(t)}.gro"
                gros.append(gro)
            parameters = [[num_cores, gros[simulation_number], cwd_project, cwd_structures, conformation_name, time_steps, photon_energy, pulse_peak, pulse_width, simulation_number + 1, gen_log] for simulation_number in range(len(gros))]

        # Otherwise we just use the normal gro-file for the conformation
        if use_time_steps == 'no':
            gro = f"{cwd_structures}/gros/{conformation_name}.gro"
            parameters = [[num_cores, gro, cwd_project, cwd_structures, conformation_name, time_steps, photon_energy, pulse_peak, pulse_width, simulation_number + 1, gen_log] for simulation_number in range(number_of_simulation)]


        # Use multiprocessing to run simulations with limited cores
        log_message(f"Starting multiprocessing pool with {num_cores} cores")
        with multiprocessing.Pool(num_cores) as pool:
            pool.map(sim_run, parameters)
    
    except Exception as e:
        log_message(f"Error in main: {e}")