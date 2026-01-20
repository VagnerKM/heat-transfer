# Transient Heat Transfer Simulation: Metal Sintering Study

This project implements a numerical solution for a transient heat transfer problem involving the sintering of metallic powders (Nickel and Copper) within a Stainless Steel mold. It was developed as part of the **EMC5417 - Heat Transfer** course at **UFSC** (Federal University of Santa Catarina).

## 📋 Problem Overview

Sintering is a metallurgical process used to manufacture metal parts by heating powders in a mold to a temperature slightly below their melting point.

In this specific study, Nickel and Copper blocks are placed inside a U-shaped Stainless Steel mold. The system is subjected to a heat source at the base, and the goal is to determine if the thermal conditions required for proper bonding are met without melting the components.

### 📐 Geometry and Schema

The domain consists of a 2D cross-section with the following layout:

* **Left/Right Walls:** Stainless Steel (15mm width)
* **Bottom Base:** Stainless Steel (15mm height)
* **Center Left:** Nickel Block (30mm x 30mm)
* **Center Right:** Copper Block (30mm x 30mm)

**Visual Representation:**

<img src="https://github.com/user-attachments/assets/3826334a-e768-401a-aea4-3e6dd5bd4e65" alt="Geometry Schema" width="600">

*Note: Points F and G represent the geometric center of the Nickel and Copper materials, respectively.*

---

## ⚙️ Physical Parameters & Boundary Conditions

### 1. Initial Conditions
The entire system starts at an initial temperature of **20°C**.

### 2. Boundary Conditions
* **Bottom Surface:** Prescribed constant temperature of **1000°C**.
* **Side Surfaces (Left/Right):** Thermally insulated (Adiabatic).
* **Top Surface:** Exposed to convection with $h = 125 W/m^2K$ and $T_{\infty} = 15^\circ C$.

### 3. Material Properties
The simulation treats the sintered material as a continuous body (approx. 99% solid density).

| Material | Density ($\rho$) $[kg/m^3]$ | Specific Heat ($c_p$) $[J/kgK]$ | Thermal Conductivity ($k$) $[W/mK]$ |
| :--- | :---: | :---: | :---: |
| **Stainless Steel** | 8000 | 500 | 16 |
| **Copper** | 8960 | 385 | 401 |
| **Nickel** | 8900 | 440 | 91 |

*(Source: Properties Table)*

---

## 🎯 The Prerogative (Sintering Constraints)

The core objective of this code is to verify if the thermal cycle meets the strict metallurgical requirements for sintering. Since Copper has a lower melting point than Nickel, it is the critical component.

**Critical Constraints:**
* **Melting Limit:** The temperature must **NOT** reach **1085°C** (Copper melting point).
* **Sintering Condition:** To ensure proper bonding, **ALL regions** of the Copper and Nickel must maintain a temperature of at least **750°C** for a minimum duration of **30 minutes**.

---

## 🚀 Objectives & Features

This repository contains the code and analysis to solve the following tasks required by the assignment:

* **Verification:** Check if the bottom boundary condition of $1000^\circ C$ satisfies the sintering constraints.
* **Correction:** If the condition is not met, determine the new required base temperature.
* **Steady State:** Estimate the time required for the system to reach thermal equilibrium.
* **Data Visualization:**
    * Temperature evolution plots for points A, B, C, D, E, F, and G.
    * Color map (heatmap) or GIF animation of the thermal distribution.

---

## 🛠️ Numerical Method

* **Discretization:** Finite Difference Method (FDM) applied to spatial and temporal domains.
* **Energy Balance:** Nodal equations derived for internal nodes, boundaries, and material interfaces.
* **Implementation:** Python (using NumPy for matrix operations and Matplotlib for plotting).

---

## 📦 Usage

To run the simulation:

```bash
# Clone the repository
git clone https://github.com/your-username/sintering-simulation.git

# Install dependencies
pip install numpy matplotlib

# Run the main script
python sintering_simulation.py

```

## 💻 Requirements

To run this simulation, ensure you have **Python 3.10+** installed along with the following scientific libraries:

* **NumPy:** Used for matrix operations and the Finite Difference Method (FDM) calculations.
* **Matplotlib:** Used for rendering the 2D heatmap and line graphs.
* **Pillow:** Required to save the thermal evolution animation as a `.gif` file.

You can install all dependencies at once using `pip`:

```bash
pip install numpy matplotlib pillow
```

## 📊 Results

Upon execution, the simulation produces visual outputs to assess the thermal behavior of the system.

### 1. Thermal Evolution Animation
The animation below demonstrates the transient heat transfer from the bottom plate (1000°C) through the Stainless Steel mold, Nickel, and Copper powders.

<p align="center">
  <img src="https://github.com/user-attachments/assets/ee50f8ab-00b2-4c2a-8d83-2e98320341dc" alt="temperature_evolution">
  <br>
  <em>Figure 1: Time-lapse of temperature distribution up to steady state.</em>
</p>

### 2. Static Analysis Plots
The script generates specific plots to verify if the sintering conditions (T ≥ 750°C for 30 min) were met.

| Final Temperature Distribution | Temperature History (Points A-G) |
| :---: | :---: |
| <img src="https://github.com/user-attachments/assets/334ca353-1e30-4f46-9e0d-b3983bf2c1d9" alt="Final_temp_distribution" width="450"> | <img src="https://github.com/user-attachments/assets/faf137b6-7362-4c5f-94bf-eda459d86500" alt="Temp_evolution" width="450"> |
| *Figure 2: Contour map at steady state.* | *Figure 3: Thermal cycles for critical points.* |

### Conclusion
The analysis of the results allows us to determine if the prescribed temperature of 1000°C is sufficient or if adjustments to the process parameters are necessary to guarantee the quality of the sintered piece.

---
