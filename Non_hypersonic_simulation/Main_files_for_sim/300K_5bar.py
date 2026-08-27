import numpy as np
import subprocess

#INPUTS----------------------------------------------------------------------------------------------------------------
# scaling factors in mean free paths for the simulation box
scale_x = 1
scale_y = 1
scale_z = 2

#this structure below can change, it was an idea of trying to control pressure in the box by varying the z direction independently
# reference conditions for fixed mean free path, sets the box in x and y direction
T_ref = 300 #K
P_ref_bar = 5 #bar

#sets the box in z direction, for consistency  T_ref = temps and P_ref_bar = pressures_bar
temps = [300] #K
pressures_bar = [5] #bar

#FOR SAICHEN==========================================================================================================
lammps_exe = "/home/lkm41/lammps-30Mar2026/build/lmp" # change your lammps file location
#=====================================================================================================================

#only inputs needed are above this line, everything below are constansts and outputs!

#--------------------------------------------------------------------------------------------------------------------

# constants
kB = 1.380649e-23
d = 3.46e-10

P_ref = P_ref_bar * 1e5
# mean free path at reference condition
lambda_mfp_ref = (kB * T_ref) / (np.sqrt(2) * np.pi * d**2 * P_ref)

# fixed box dimensions
Lx = scale_x * lambda_mfp_ref
Ly = scale_y * lambda_mfp_ref

# convert to Å
Lx_A = Lx * 1e10
Ly_A = Ly * 1e10

# snap to Ti lattice spacing
Lx_A = Lx_A - Lx_A % 2.95
Ly_A = Ly_A - Ly_A % 2.95

for P_bar in pressures_bar:

    P = P_bar * 1e5

    for T in temps:

        # variable mean free path
        lambda_mfp = (kB * T) / (np.sqrt(2) * np.pi * d**2 * P)

        # only Lz changes
        Lz = scale_z * lambda_mfp
        Lz_A = Lz * 1e10

        # fixed x/y dimensions in meters
        Lx_m = Lx_A * 1e-10
        Ly_m = Ly_A * 1e-10

        # volume
        V = Lx_m * Ly_m * Lz

        # molecule count
        N = int((P * V) / (kB * T))
        N = max(N, 1)

        print(f"T={T}K | P={P_bar}bar | "
              f"Lx={Lx_A:.1f}Å Ly={Ly_A:.1f}Å "
              f"Lz={Lz_A:.1f}Å | N={N}")

        command = [
            "srun", #FOR HCP
            lammps_exe,
            "-var", "T", str(T),
            "-var", "P", str(P_bar),
            "-var", "Lx", str(Lx_A),
            "-var", "Ly", str(Ly_A),
            "-var", "Lz", str(Lz_A),
            "-var", "Lz_minus_one", str(Lz_A - 1),
             "-var", "N", str(N),
            "-in", "300K_5bar.in"
        ]

        subprocess.run(command, check=True)

print("\nAll simulations finished.")
