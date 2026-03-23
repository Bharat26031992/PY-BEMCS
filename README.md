#  Interactive MCC Charge Exchange Simulator (MATLAB)

A highly optimized, fully vectorized Monte Carlo Collision (MCC) simulator built in MATLAB. This application models the dynamics of **Charge Exchange (CEX)** collisions in a plasma thruster plume (e.g., Ion or Hall thrusters) and features a fully interactive GUI with real-time analytics, directional flux probes, and GIF recording.

*Adapted and heavily upgraded from the original Java-based PIC/MCC demonstration by [Particle In Cell Consulting](http://www.particleincell.com/2011/mcc/).*


## ✨ Key Features

* **⚡ Vectorized Physics Engine:** Replaces traditional loops with matrix operations for massive performance gains in MATLAB.
* **📊 Live 2D & 1D Analytics:** Instantly visualizes the plume using a 2D spatial heatmap (`histcounts2`) and a 1D axial density profile.
* **🎯 Directional Flux Probes:** Point-and-click to place physical probes in the simulation domain. 
    * Adjustable **collection radius** and **surface angle**.
    * Probes act as physical barriers, absorbing particles and creating realistic downstream "shadows."
* **📈 Real-Time Time Series:** Tracks the particle flux (hits per time-step) across all placed probes on a live-scrolling plot. Includes instant toggles for **Log X** and **Log Y** scales.
* **💾 Data Export:** Pause the simulation at any time and export the raw probe flux histories directly to a nicely formatted `.csv` file for post-processing.
* **🎥 Built-in GIF Recorder:** Silently buffers frames into memory while running, allowing you to export high-quality `.gif` animations of the plume dynamics with a single click.

---

## 🛠️ Installation & Usage

1. **Prerequisites:** You need **MATLAB R2015b or newer** (requires `histcounts2`, `gobjects`, and standard UI components). No extra toolboxes are required.
2. **Download:** Clone this repository or download the `MCC_App.m` file.
3. **Run:** Open MATLAB, navigate to the folder containing the file, and type the following in the Command Window:
   ```matlab
   MCC_App

   
