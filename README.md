# 2D Heat Transfer Simulation: Sintering Process

## Overview

This project simulates the transient heat transfer in a metallurgical sintering process. It models a **Stainless Steel mold** containing **Nickel** and **Copper** powders. The goal is to determine if the copper reaches the sintering temperature (750°C) for at least 30 minutes without melting ($T_{melt} = 1085^\circ C$).

The simulation solves the **2D Heat Conduction Equation** using the **Finite Difference Method (Explicit)**.

## Key Features

- **Multi-Material Mesh:** Automatically assigns thermal properties ($k, \rho, c_p$) based on geometry.
- **Stability Control:** Checks the Fourier number criterion before execution.
- **Visualization:**
  - **Heatmap:** Final temperature distribution.
  - **History Plot:** Temperature vs. Time for critical points.
  - **Animation (GIF):** A time-lapse of the heating process.

## Requirements

- Python 3.x
- Libraries: `numpy`, `matplotlib`
- For GIF generation: `pillow` (usually included with matplotlib) or `imagemagick`.

## Usage

Run the script directly:

```bash
python sintering_simulation.py
```
