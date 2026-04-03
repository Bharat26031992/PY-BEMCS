import numpy as np
import matplotlib
matplotlib.use('Agg') # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from PIL import Image
from physics_engine import DigitalTwinSimulator

def run_rf_sweep():
    # Vary RF potential from 50 to 900V. Using step size of 50V
    rf_amplitudes = np.arange(50, 950, 50) 
    
    for rf_amp in rf_amplitudes:
        print(f"Starting simulation for RF Amplitude: {rf_amp}V")
        sim = DigitalTwinSimulator()
        
        params = {
            'sim_mode': 'Both',
            # Default value of screen grid and second grid dimensions
            'ts': 1.0, 'ta': 1.0, 'rs': 1.0, 'ra': 0.6,
            'cham_s': 0, 'cham_a': 0,
            # Keep the screen grid at 1000V
            'Vs': 1000, 'Va': 0,
            'gap': 1.0,
            # Plasma density to be around 5e17m3
            'n0_plasma': 2e17,
            'Te_up': 3.0, 'Ti': 0.1,
            # Default neutralizer current to be zero
            'neut_rate': 0, 
            # Neutral density to be around 1e14
            'n0': 1e14,
            # RF applied only to the second grid
            'rf_enable': True,
            'rf_grid_idx': 1,
            'rf_freq': 13.56,
            'rf_amp': rf_amp
        }
        
        # Build grid config expected by the physics engine
        params['grids'] = [
            {'V': params['Vs'], 't': params['ts'], 'gap': params['gap'], 'r': params['rs'], 'cham': params['cham_s']},
            {'V': params['Va'], 't': params['ta'], 'gap': 1.0, 'r': params['ra'], 'cham': params['cham_a']}
        ]
        
        sim.build_domain(params)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        canvas = FigureCanvas(fig)
        frames = []
        
        gy, gx = np.where(sim.isBound)
        
        # Run script for 1000 timesteps
        total_steps = 500
        for step_idx in range(1, total_steps + 1):
            sim.step(params)
            
            # Save frame every 1 steps to create the animation without memory bloat
            if step_idx % 1 == 0:
                ax.clear()
                ax.set_title(f'RF Amplitude: {rf_amp}V - Step: {step_idx}')
                
                # Plot potential field
                ax.contourf(sim.X, sim.Y, sim.V, 20, cmap='viridis', alpha=0.4)
                
                # Plot solid grids
                ax.scatter(gx * sim.dx, gy * sim.dy, s=12, c='k', alpha=0.8)
                
                # Separate particles
                prim_mask = ~sim.p_isCEX
                cex_mask = sim.p_isCEX
                
                # Plot primary ions mapped by energy
                if np.any(prim_mask):
                    v_sq_prim = sim.p_vx[prim_mask]**2 + sim.p_vy[prim_mask]**2
                    energy = (0.5 * sim.m_XE * v_sq_prim) / sim.q
                    ax.scatter(sim.p_x[prim_mask], sim.p_y[prim_mask], s=2, c=energy, cmap='turbo', vmin=0, vmax=1050, alpha=0.8)
                
                # Plot charge-exchange ions
                if np.any(cex_mask):
                    ax.scatter(sim.p_x[cex_mask], sim.p_y[cex_mask], s=7, c='red', alpha=1.0)
                
                # Plot electrons
                if len(sim.e_x) > 0:
                    ax.scatter(sim.e_x, sim.e_y, s=1, c='#00FF00', alpha=0.5)
                    
                ax.set_xlim(0, sim.Lx)
                ax.set_ylim(0, sim.Ly)
                ax.set_xlabel('Axial Distance Z (mm)')
                ax.set_ylabel('Radial Distance R (mm)')
                
                canvas.draw()
                rgba = np.asarray(canvas.buffer_rgba())
                frames.append(rgba)
                
                if step_idx % 250 == 0:
                    print(f"  Step {step_idx}/{total_steps} completed")
        
        plt.close(fig)
        
        # Save each gif animation
        if frames:
            pil_images = [Image.fromarray(f) for f in frames]
            filename = f"rf_sweep_amp_{rf_amp}V.gif"
            pil_images[0].save(filename, save_all=True, append_images=pil_images[1:], optimize=False, duration=50, loop=0)
            print(f"Saved animation to {filename}\n")

if __name__ == '__main__':
    run_rf_sweep()
