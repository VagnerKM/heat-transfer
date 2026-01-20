"""
SINTERING SIMULATION - 2D FINITE DIFFERENCE METHOD
--------------------------------------------------
This script models the heat transfer in a mold containing Nickel and Copper.
It produces:
1. A static contour map of the final state.
2. A temperature history plot for specific points.
3. A GIF animation of the heat propagation.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib.animation as animation
import os
import time


# Start execution timer
start_time = time.time()

# ==========================================
# 1. GEOMETRIC PARAMETERS
# ==========================================
# Dimensions converted to meters [m]
total_height = (15 + 30) / 1000    # [0.045 m]
total_width = ((2 * 15) + (2 * 30)) / 1000 # [0.090 m]

# ==========================================
# 2. INITIAL & BOUNDARY CONDITIONS
# ==========================================
T_initial = 20.0       # Initial temperature of all materials [°C]
T_bottom_surf = 1000.0 # Prescribed temperature at the bottom surface [°C]
T_ambient = 15.0       # Ambient temperature for convection at top [°C]
h = 125.0              # Convection heat transfer coefficient [W/(m²·K)]

# ==========================================
# 3. MATERIAL PROPERTIES
# ==========================================

# Stainless Steel (Mold)
rho_ss = 8000.0    # Density [kg/m³]
cp_ss = 500.0      # Specific Heat [J/(kg·K)]
k_ss = 16.0        # Thermal Conductivity [W/(m·K)]

# Nickel (Sample 1)
rho_ni = 8900.0    # Density [kg/m³]
cp_ni = 440.0      # Specific Heat [J/(kg·K)]
k_ni = 91.0        # Thermal Conductivity [W/(m·K)]

# Copper (Sample 2)
rho_cu = 8960.0    # Density [kg/m³]
cp_cu = 385.0      # Specific Heat [J/(kg·K)]
k_cu = 401.0       # Thermal Conductivity [W/(m·K)]

# ==========================================
# 4. DISCRETIZATION (Spatial & Temporal)
# ==========================================

# Spatial Discretization (dx = dy)
S = 2          # Mesh subdivision factor
M = 6 * S      # Number of divisions in X (Must be multiple of 6)
N = 3 * S      # Number of divisions in Y (Must be multiple of 3)
dx = total_width / M
dy = total_height / N
l = dx         # Characteristic length for equations

# Validate Square Mesh
if not np.isclose(dx, dy):
    raise ValueError("The formulation requires dx = dy. Please adjust M and N.")

# Temporal Discretization
dt = 0.01      # Time step [s]
total_simulation_time = 1200  # Total simulation time [s]
num_time_steps = int(total_simulation_time / dt)

# ==========================================
# 5. STABILITY CRITERIA CHECK
# ==========================================
# Fourier number coefficients for all node types (internal, borders, corners)
C = np.zeros((15)) 

# Coefficients calculation (Heat equation stability terms)
C[0] = 1 - 4 * (dt * k_ss) / (rho_ss * cp_ss * dx**2)
C[1] = 1 - 4 * (dt * k_ni) / (rho_ni * cp_ni * dx**2)
C[2] = 1 - 4 * (dt * k_cu) / (rho_cu * cp_cu * dx**2)
C[3] = 1 - 4 * (dt * k_ss) / (rho_ss * cp_ss * dx**2) - 2 * dt * h / (rho_ss * cp_ss * dx)
C[4] = 1 - 4 * (dt * k_ni) / (rho_ni * cp_ni * dx**2) - 2 * dt * h / (rho_ni * cp_ni * dx)
C[5] = 1 - 4 * (dt * k_cu) / (rho_cu * cp_cu * dx**2) - 2 * dt * h / (rho_cu * cp_cu * dx)
C[6] = 1 - 4 * (k_ss + k_ni + h * dx) * (dt / (dx**2 * (rho_ss * cp_ss + rho_ni * cp_ni)))
C[7] = 1 - 4 * (k_ni + k_cu + h * dx) * (dt / (dx**2 * (rho_ni * cp_ni + rho_cu * cp_cu)))
C[8] = 1 - 4 * (k_cu + k_ss + h * dx) * (dt / (dx**2 * (rho_cu * cp_cu + rho_ss * cp_ss)))
C[9] = 1 - 4 * (3 * k_ss + k_ni) * (dt / (dx**2 * (3 * (rho_ss * cp_ss) + rho_ni * cp_ni)))
C[10] = 1 - 4 * (k_ni + k_cu + k_ss) * (dt / (dx**2 * (rho_ni * cp_ni + rho_cu * cp_cu + rho_ss * cp_ss)))
C[11] = 1 - 4 * (k_ni + k_ss) * (dt / (dx**2 * (rho_ni * cp_ni + rho_ss * cp_ss)))
C[12] = 1 - 4 * (k_cu + k_ni) * (dt / (dx**2 * (rho_cu * cp_cu + rho_ni * cp_ni)))
C[13] = 1 - 4 * (k_cu + k_ss) * (dt / (dx**2 * (rho_cu * cp_cu + rho_ss * cp_ss)))
C[14] = 1 - 4 * (3 * k_ss + k_cu) * (dt / (dx**2 * (3 * (rho_ss * cp_ss) + rho_cu * cp_cu)))

# Check convergence
if np.any(C < 0):
    print(f"WARNING: Fourier stability criterion NOT met!\nMin Value: {C.min()} < 0")
else:
    print(f"Stability criterion met!\nMin Value: {C.min()} > 0")

# ==========================================
# 6. MESH & MATERIAL MAPPING
# ==========================================

# Initialize property maps with Stainless Steel (Base material)
k_map = np.full((N + 1, M + 1), k_ss)
rho_cp_map = np.full((N + 1, M + 1), rho_ss * cp_ss)

# Define indices for material boundaries based on geometry
i_ni_start = int(round((15/1000) / dx))
i_cu_start = int(round((15 + 30)/1000 / dx))
i_cu_end = int(round((15 + 30 + 30)/1000 / dx))
j_block_start = int(round((15/1000) / dy))

# Apply Nickel properties to its region
k_map[j_block_start:, i_ni_start:i_cu_start] = k_ni
rho_cp_map[j_block_start:, i_ni_start:i_cu_start] = rho_ni * cp_ni

# Apply Copper properties to its region
k_map[j_block_start:, i_cu_start:i_cu_end] = k_cu
rho_cp_map[j_block_start:, i_cu_start:i_cu_end] = rho_cu * cp_cu

# ==========================================
# 7. INITIALIZE SIMULATION VARIABLES
# ==========================================

# Temperature Matrix (T)
T = np.full((N + 1, M + 1), T_initial, dtype=float)

# Apply Boundary Condition: Bottom Surface (Prescribed Temp)
T[0, :] = T_bottom_surf

# Define Coordinates for Points of Interest (A-G)
j_top = N
i_A = 0
i_B = i_ni_start
i_C = i_cu_start
i_D = i_cu_end
i_E = M
# Center points for materials
i_F = int(round((15 + 15)/1000 / dx))
j_F = int(round((15 + 15)/1000 / dy))
i_G = int(round((15 + 30 + 15)/1000 / dx))
j_G = int(round((15 + 15)/1000 / dy))

# History Lists for Plotting
time_vector = []
T_hist_A = []
T_hist_B = []
T_hist_C = []
T_hist_D = []
T_hist_E = []
T_hist_F = []
T_hist_G = []

# Lists for Animation (Snapshots)
# We calculate how many steps represent 1 second

steps_per_frame = int(1 / dt) 
T_snapshots = []      # Will store the full matrix T at specific intervals
snapshot_times = []   # Will store the time corresponding to each snapshot

# ==========================================
# 8. MAIN CALCULATION LOOP (Explicit Method)
# ==========================================

print(f"Starting simulation with M={M}, N={N}, dt={dt:.4f}s...")

for step in range(num_time_steps):
    T_old = T.copy() # Store temperature from previous time step
    current_time = step * dt

    # ---  Capture Snapshot for Animation ---
    # Save a frame every 1 second of simulation time to save memory
    if step % steps_per_frame == 0:
        T_snapshots.append(T_old.copy())
        snapshot_times.append(current_time)

    # Iterate through the grid
    for j in range(1, N + 1):   # Rows (y), bottom to top
        for i in range(M + 1):  # Columns (x), left to right
            k = k_map[j, i]
            rho_cp = rho_cp_map[j, i]
            tau = (k * dt) / (rho_cp * l**2) # Fourier number for this node

            # Nodal Equations based on position (Internal vs Boundary)
            
            if 1 <= j < N and 1 <= i < M: # Internal Node
                T[j,i] = T_old[j,i]*(1-4*tau) + \
                        tau*(T_old[j,i-1] + T_old[j,i+1] + T_old[j-1,i] + T_old[j+1,i])
            
            elif j == N and 1 <= i < M: # Top Boundary (Convection)
                T[j,i] = T_old[j,i] * (1 - 4*tau - (2*h*l*tau)/k) + \
                        tau * (T_old[j,i-1] + T_old[j,i+1] + 2*T_old[j-1,i] + (2*h*l*T_ambient)/k)
            
            elif i == 0 and 1 <= j < N: # Left Boundary (Insulated)
                T[j,i] = T_old[j,i]*(1-4*tau) + \
                        tau*(2*T_old[j,i+1] + T_old[j-1,i] + T_old[j+1,i])
            
            elif i == M and 1 <= j < N: # Right Boundary (Insulated)
                T[j,i] = T_old[j,i]*(1-4*tau) + \
                        tau*(2*T_old[j,i-1] + T_old[j-1,i] + T_old[j+1,i])
            
            elif i == 0 and j == N: # Top-Left Corner (Insulated + Convection)
                T[j,i] = T_old[j,i] * (1 - 4*tau - (2*h*l*tau)/k) + \
                        tau * (2*T_old[j,i+1] + 2*T_old[j-1,i] + (2*h*l*T_ambient)/k)
            
            elif i == M and j == N: # Top-Right Corner (Insulated + Convection)
                T[j,i] = T_old[j,i] * (1 - 4*tau - (2*h*l*tau)/k) + \
                        tau * (2*T_old[j,i-1] + 2*T_old[j-1,i] + (2*h*l*T_ambient)/k)

    # Store results for Line Graph
    if step % 100 == 0: 
        time_vector.append(current_time)
        T_hist_A.append(T[j_top, i_A])
        T_hist_B.append(T[j_top, i_B])
        T_hist_C.append(T[j_top, i_C])
        T_hist_D.append(T[j_top, i_D])
        T_hist_E.append(T[j_top, i_E])
        T_hist_F.append(T[j_F, i_F])
        T_hist_G.append(T[j_G, i_G])

        # Status Update
        print(f"Time: {current_time:.2f}s | T_A = {T[j_top, i_A]:.2f}°C | T_F (Ni) = {T[j_F, i_F]:.2f}°C | T_G (Cu) = {T[j_G, i_G]:.2f}°C")

# Finalize timer
end_time = time.time()
print(f"Simulation completed in {end_time - start_time:.2f} seconds.")

# ==========================================
# 9. GIF GENERATION
# ==========================================
print("\n--- Starting GIF Generation ---")

# We will take 1 frame every 10 simulation seconds (approx 120 frames total)
# T_snapshots and snapshot_times were filled in the loop (1 per sec)

step_skip = 1  # Skip step (e.g., 0s, 10s, 20s...)
selected_frames = T_snapshots[::step_skip]
selected_times = snapshot_times[::step_skip]

if not selected_frames:
    print("ERROR: No frames captured. Check the loop logic.")
else:
    print(f"Processing {len(selected_frames)} frames for the GIF...")

    # Figure Setup
    fig_anim, ax_anim = plt.subplots(figsize=(10, 5))

    # Initial Plot
    # vmin/vmax are fixed so colors don't flicker
    im_anim = ax_anim.imshow(selected_frames[0], cmap='inferno', interpolation='bilinear', 
                            origin='lower', vmin=T_initial, vmax=T_bottom_surf)

    plt.colorbar(im_anim, ax=ax_anim, label='Temperature (°C)')
    ax_anim.set_title('Temperature Evolution')
    ax_anim.set_xlabel('Nodes (x)')
    ax_anim.set_ylabel('Nodes (y)')

    # Update Function
    def update(frame_idx):
        im_anim.set_data(selected_frames[frame_idx])
        ax_anim.set_title(f'Temperature Evolution\nTime: {selected_times[frame_idx]:.2f} s')
        return [im_anim]

    # Create Animation
    # interval=100ms means 10 frames per second
    anim = animation.FuncAnimation(fig_anim, update, frames=len(selected_frames), interval=100, blit=False)

    # Full path to save
    filename = 'temperature_evolution.gif'
    full_path = os.path.abspath(filename)

    print("Saving file... (Please wait, this might take a few seconds)")
    
    try:
        # Explicitly use PillowWriter
        writer = animation.PillowWriter(fps=10) 
        anim.save(full_path, writer=writer)
        print(f"✅ SUCCESS! GIF saved at:\n--> {full_path}")
    except Exception as e:
        print(f"❌ ERROR saving GIF: {e}")
        print("Tip: Ensure pillow is installed (pip install pillow)")

    # Correct cleanup to avoid 'AttributeError'
    anim.event_source = None
    anim = None
    plt.close(fig_anim)

# ==========================================
# 10. RESULTS & VISUALIZATION (Static Plots)
# ==========================================

# --- PLOT 1: Contour Map (Heatmap) ---
plt.figure(figsize=(12, 7))

x_coords = np.linspace(0, total_width * 1000, M + 1) 
y_coords = np.linspace(0, total_height * 1000, N + 1)
X, Y = np.meshgrid(x_coords, y_coords)

contour = plt.contourf(X, Y, T, levels=50, cmap='inferno')
plt.colorbar(contour, label='Temperature (°C)')
plt.title(f'Final Temperature Distribution (t = {total_simulation_time}s)', y=1.05)
plt.xlabel('Width (mm)')
plt.ylabel('Height (mm)')

style_contour = {'color': 'white', 'linestyle': '--', 'linewidth': 0.5}
plt.plot([0, 90, 90, 0, 0], [0, 0, 45, 45, 0], **style_contour)
plt.plot([15, 75], [15, 15], **style_contour)
plt.plot([15, 15], [15, 45], **style_contour)
plt.plot([45, 45], [15, 45], **style_contour)
plt.plot([75, 75], [15, 45], **style_contour)

points_def = [
    (i_A, j_top, 'Point A'), (i_B, j_top, 'Point B'), (i_C, j_top, 'Point C'),
    (i_D, j_top, 'Point D'), (i_E, j_top, 'Point E'), (i_F, j_F, 'Point F'),
    (i_G, j_G, 'Point G')
]

for i_idx, j_idx, label_txt in points_def:
    plt.plot(x_coords[i_idx], y_coords[j_idx], 'o', markersize=4, color='white', markeredgecolor='black', label=label_txt)

offset_x = (total_width * 1000) * 0.02
offset_y = (total_height * 1000) * 0.05
outline_effect = [path_effects.withStroke(linewidth=3, foreground='black')]

plt.text(x_coords[i_A], y_coords[j_top] + offset_y, 'A', color='white', ha='center', va='bottom', fontsize=10, fontweight='bold', path_effects=outline_effect)
plt.text(x_coords[i_B], y_coords[j_top] + offset_y, 'B', color='white', ha='center', va='bottom', fontsize=10, fontweight='bold', path_effects=outline_effect)
plt.text(x_coords[i_C], y_coords[j_top] + offset_y, 'C', color='white', ha='center', va='bottom', fontsize=10, fontweight='bold', path_effects=outline_effect)
plt.text(x_coords[i_D], y_coords[j_top] + offset_y, 'D', color='white', ha='center', va='bottom', fontsize=10, fontweight='bold', path_effects=outline_effect)
plt.text(x_coords[i_E], y_coords[j_top] + offset_y, 'E', color='white', ha='center', va='bottom', fontsize=10, fontweight='bold', path_effects=outline_effect)
plt.text(x_coords[i_F] + offset_x, y_coords[j_F], 'F', color='white', ha='left', va='center', fontsize=10, fontweight='bold', path_effects=outline_effect)
plt.text(x_coords[i_G] + offset_x, y_coords[j_G], 'G', color='white', ha='left', va='center', fontsize=10, fontweight='bold', path_effects=outline_effect)

x_min, x_max = 0, total_width * 1000
y_min, y_max = 0, total_height * 1000
margin_x, margin_y = (x_max - x_min) * 0.05, (y_max - y_min) * 0.05
plt.xlim(x_min - margin_x, x_max + margin_x)
plt.ylim(y_min, y_max + margin_y)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=7, fontsize='small')
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout(pad=1.5)
plt.show()

# --- PLOT 2: Temperature vs. Time History ---
plt.figure(figsize=(12, 7))
plt.plot(time_vector, T_hist_A, label='Point A', color='cyan', linestyle='--')
plt.plot(time_vector, T_hist_B, label='Point B', color='lime', linestyle='--')
plt.plot(time_vector, T_hist_C, label='Point C', color='magenta', linestyle='--')
plt.plot(time_vector, T_hist_D, label='Point D', color='yellow', linestyle='--')
plt.plot(time_vector, T_hist_E, label='Point E', color='red', linestyle='--')
plt.plot(time_vector, T_hist_F, label='Point F', color='blue', linewidth=2)
plt.plot(time_vector, T_hist_G, label='Point G', color='orangered', linewidth=2)

if time_vector:
    full_history = np.array([T_hist_A, T_hist_B, T_hist_C, T_hist_D, T_hist_E, T_hist_F, T_hist_G])
    min_temp_curve = np.min(full_history, axis=0)
    indices_750 = np.where(min_temp_curve >= 750)[0]
    time_750 = None
    if len(indices_750) > 0:
        time_750 = time_vector[indices_750[0]]

    time_stable = None
    stable_duration = 300 
    tol_temp = 0.1 

    if time_vector[-1] > stable_duration:
        for i in range(len(time_vector)):
            current_t = time_vector[i]
            future_t = current_t + stable_duration
            try:
                future_idx = next(j for j, t in enumerate(time_vector) if t >= future_t)
            except StopIteration:
                break
            
            max_delta = 0
            for p in range(full_history.shape[0]):
                delta = abs(full_history[p, future_idx] - full_history[p, i])
                if delta > max_delta:
                    max_delta = delta
            
            if max_delta <= tol_temp:
                time_stable = current_t
                break

    style_ref = {'color': 'gray', 'linestyle': '--', 'linewidth': 0.5}
    plt.axhline(y=750, **style_ref)

    if time_750 is not None:
        plt.axvline(x=time_750, **style_ref)
        plt.text(time_750 + 5, 700, f'750ºC - {time_750:.0f} s', color='gray', fontsize=11)

    if time_stable is not None:
        style_stable = {'color': 'green', 'linestyle': '-.', 'linewidth': 1.5}
        plt.axvline(x=time_stable, **style_stable)
        plt.text(time_stable + 5, 500, f'Steady State\nt = {time_stable:.0f} s', color='green', fontsize=11)

plt.title('Temperature Evolution at Interest Points')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (°C)')
plt.grid(True)
plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
plt.tight_layout()
plt.show()
