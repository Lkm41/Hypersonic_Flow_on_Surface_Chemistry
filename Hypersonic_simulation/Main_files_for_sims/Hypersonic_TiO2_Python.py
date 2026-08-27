import numpy as np
import subprocess

#inputs-----------------------------------------------------------------------------------
#Chosen Mach Number
MachN = 5

#Chosen box scaling factors to mean free path
scale_x = 16
scale_y = 7
scale_z = 7

#z and y direction Ti scaling factor to box mean free path box scale
Ti_scale = 0.5714285714285714 #value between 0 and 1 ## 1 means ti covers 100% of box dimension. This value gives a slab 4 Å in z and y
#scale of slab in y and z /scale of box in y and z = 4/7 = 0.5714285714285714

#This could be changed, its carry over from old code!
# T_ref and P_ref_bar are your inputs to calculate the mean free path in x and y direction
#can change or just ensure T_ref = temps and P_ref_bar = pressures_bar to create a symmetrical box
T_ref = 300 #K
P_ref_bar = 5 #bar
#this is for the z direction
temps = [300] #K
pressures_bar = [5] #bar

#this controls the thickness of the slab in x direction by controlling each of the three layers
dx_upper  = 8   # Å
dx_middle = 4   # Å
dx_bottom = 4   # Å

#FOR SAICHEN==========================================================================================================
lammps_exe = "/home/lkm41/lammps-30Mar2026/build/lmp" # saichen and hpc
#======================================================================================================================

#everything below here are constants and outputs, the only inputs to run the simulation are above!
#---------------------------------------------------------------------------------------------------------------------
# constants
kB = 1.380649e-23
d = 3.46e-10

#Mean free path calulcation
P_ref = P_ref_bar * 1e5

# mean free path at reference condition
lambda_mfp_ref = (kB * T_ref) / (np.sqrt(2) * np.pi * d**2 * P_ref)

# fixed box dimensions
Lx = scale_x * lambda_mfp_ref
Ly = scale_y * lambda_mfp_ref

#-----------------------------------------------------------------------------------------------------------------

#Ti box dimesnions

# Box centre coordinates
x_centre = (scale_x * lambda_mfp_ref * 1e10) / 2
y_centre = (scale_y * lambda_mfp_ref * 1e10) / 2
z_centre = (scale_z * lambda_mfp_ref * 1e10) / 2

# Titanium block dimensions
Lx_box = scale_x * lambda_mfp_ref * 1e10
Ly_box = scale_y * lambda_mfp_ref * 1e10
Lz_box = scale_z * lambda_mfp_ref * 1e10

# ============================================================
# Titanium slab dimensions
# ============================================================

dx = dx_upper + dx_middle + dx_bottom   # Total = 16 Å

dy = Ly_box * Ti_scale
dz = Lz_box * Ti_scale

# ============================================================
# Titanium block (entire slab)
# ============================================================

Lx_Ti_zone1 = x_centre - dx/2
Ly_Ti_zone1 = y_centre - dy/2
Lz_Ti_zone1 = z_centre - dz/2

Lx_Ti_zone2 = x_centre + dx/2
Ly_Ti_zone2 = y_centre + dy/2
Lz_Ti_zone2 = z_centre + dz/2

# ============================================================
# Upper layer (8 Å)
# ============================================================

Lx_Ti_upper1 = Lx_Ti_zone1
Lx_Ti_upper2 = Lx_Ti_upper1 + dx_upper

Ly_Ti_upper1 = Ly_Ti_zone1
Ly_Ti_upper2 = Ly_Ti_zone2

Lz_Ti_upper1 = Lz_Ti_zone1
Lz_Ti_upper2 = Lz_Ti_zone2

# ============================================================
# Middle layer (4 Å)
# ============================================================

Lx_Ti_middle1 = Lx_Ti_upper2
Lx_Ti_middle2 = Lx_Ti_middle1 + dx_middle

Ly_Ti_middle1 = Ly_Ti_zone1
Ly_Ti_middle2 = Ly_Ti_zone2

Lz_Ti_middle1 = Lz_Ti_zone1
Lz_Ti_middle2 = Lz_Ti_zone2

# ============================================================
# Bottom layer (4 Å)
# ============================================================

Lx_Ti_bottom1 = Lx_Ti_middle2
Lx_Ti_bottom2 = Lx_Ti_zone2

Ly_Ti_bottom1 = Ly_Ti_zone1
Ly_Ti_bottom2 = Ly_Ti_zone2

Lz_Ti_bottom1 = Lz_Ti_zone1
Lz_Ti_bottom2 = Lz_Ti_zone2

# ============================================================
# Observation region (3 Å in front of slab)
# ============================================================

Lx_Ti_obs1 = Lx_Ti_zone1 - 3
Lx_Ti_obs2 = Lx_Ti_zone1

Ly_Ti_obs1 = Ly_Ti_zone1
Ly_Ti_obs2 = Ly_Ti_zone2

Lz_Ti_obs1 = Lz_Ti_zone1
Lz_Ti_obs2 = Lz_Ti_zone2
#---------------------------------------------------------------------------------------------

# convert to Å
Lx_A = Lx * 1e10
Ly_A = Ly * 1e10

# snap to Ti lattice spacing
Lx_A = Lx_A - Lx_A % 2.95
Ly_A = Ly_A - Ly_A % 2.95

#velocities
Mach_speed = np.sqrt(T_ref*259.8*1.286)*(MachN) # m/s
Vxmin = Mach_speed/100 # Å/ps - needed to convert m/s to Å/ps for lammps units metal #1 Å/ps = 100 m/s
Vxmax = (Mach_speed + 0.1)/100 # Å/ps


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

        # volume - current changes with temperature and pressure
        V = Lx_m * Ly_m * Lz

        # molecule count
        N = int((P * V) / (kB * T))
        N = max(N, 1)

        print(f"T={T}K | P={P_bar}bar | "
              f"Lx={Lx_A:.1f}Å Ly={Ly_A:.1f}Å "
              f"Lz={Lz_A:.1f}Å | N={N}")

        command = [
            "srun",
            lammps_exe,

            "-var", "T", str(T),
            "-var", "P", str(P_bar),

            "-var", "Lx", str(Lx_A),
            "-var", "Ly", str(Ly_A),
            "-var", "Lz", str(Lz_A),

            "-var", "Lx_outlet",
            str((Lx_A / 2) - ((dx / 2) + 100)),

            "-var", "Lz_minus_one", str(Lz_A - 1),
            "-var", "Lx_ten", "10",

            "-var", "N", str(N),
            "-var", "N_more", str(N * 10),

            "-var", "Vxmin", str(Vxmin),
            "-var", "Vxmax", str(Vxmax),

            # ========================================================
            # Entire titanium slab
            # ========================================================

            "-var", "Lx_Ti_zone1", str(Lx_Ti_zone1),
            "-var", "Ly_Ti_zone1", str(Ly_Ti_zone1),
            "-var", "Lz_Ti_zone1", str(Lz_Ti_zone1),

            "-var", "Lx_Ti_zone2", str(Lx_Ti_zone2),
            "-var", "Ly_Ti_zone2", str(Ly_Ti_zone2),
            "-var", "Lz_Ti_zone2", str(Lz_Ti_zone2),

            # ========================================================
            # Upper titanium layer: 8 Å
            # ========================================================

            "-var", "Lx_Ti_upper1", str(Lx_Ti_upper1),
            "-var", "Ly_Ti_upper1", str(Ly_Ti_upper1),
            "-var", "Lz_Ti_upper1", str(Lz_Ti_upper1),

            "-var", "Lx_Ti_upper2", str(Lx_Ti_upper2),
            "-var", "Ly_Ti_upper2", str(Ly_Ti_upper2),
            "-var", "Lz_Ti_upper2", str(Lz_Ti_upper2),

            # ========================================================
            # Middle titanium layer: 4 Å
            # ========================================================

            "-var", "Lx_Ti_middle1", str(Lx_Ti_middle1),
            "-var", "Ly_Ti_middle1", str(Ly_Ti_middle1),
            "-var", "Lz_Ti_middle1", str(Lz_Ti_middle1),

            "-var", "Lx_Ti_middle2", str(Lx_Ti_middle2),
            "-var", "Ly_Ti_middle2", str(Ly_Ti_middle2),
            "-var", "Lz_Ti_middle2", str(Lz_Ti_middle2),

            # ========================================================
            # Bottom titanium layer: 4 Å
            # ========================================================

            "-var", "Lx_Ti_bottom1", str(Lx_Ti_bottom1),
            "-var", "Ly_Ti_bottom1", str(Ly_Ti_bottom1),
            "-var", "Lz_Ti_bottom1", str(Lz_Ti_bottom1),

            "-var", "Lx_Ti_bottom2", str(Lx_Ti_bottom2),
            "-var", "Ly_Ti_bottom2", str(Ly_Ti_bottom2),
            "-var", "Lz_Ti_bottom2", str(Lz_Ti_bottom2),

            # ========================================================
            # Physisorption observation region
            # ========================================================

            "-var", "Lx_Ti_obs1", str(Lx_Ti_obs1),
            "-var", "Ly_Ti_obs1", str(Ly_Ti_obs1),
            "-var", "Lz_Ti_obs1", str(Lz_Ti_obs1),

            "-var", "Lx_Ti_obs2", str(Lx_Ti_obs2),
            "-var", "Ly_Ti_obs2", str(Ly_Ti_obs2),
            "-var", "Lz_Ti_obs2", str(Lz_Ti_obs2),

            "-in", "Hypersonic_TiO2.in",
        ]


        subprocess.run(command, check=True)

print("\nAll simulations finished.")