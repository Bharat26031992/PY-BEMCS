import numpy as np
import matplotlib
# matplotlib.use('TkAgg') # Uncomment if window stays white
import matplotlib.pyplot as plt
import multiprocessing as mp
import queue
from physics_engine import DigitalTwinSimulator 

# --- WORKER PROCESS FUNCTION ---
def worker_sweep(gap, n0_sweep, result_queue):
    q = 1.602e-19
    m_XE = 131.293 * 1.6605e-27
    
    for n0 in n0_sweep:
        sim = DigitalTwinSimulator()
        params = {
            'sim_mode': 'Both',
            'ts': 1, 'ta': 1, 'rs': 1.0, 'ra': 0.6,
            'cham_s': 0, 'cham_a': 0,
            'Vs': 1650, 'Va': -350, 
            'gap': gap,
            'n0_plasma': n0,
            'Te_up': 3.0, 'Ti': 0.1,
            'neut_rate': 0, 'n0': 0 
        }
        
        sim.build_domain(params)
        steady_state_steps = 600
        impingement_history = []
        
        for step_idx in range(1, steady_state_steps + 1):
            # We assume your sim.step returns or updates a way to count lost particles
            # Here we extract it via the 'I_accel' telemetry if your engine provides it, 
            # or calculate it from the internal particle tracker.
            _, _, _, I_accel_step, _ = sim.step(params)
            
            # We only start averaging after the initial plasma startup (e.g., step 500)
            if step_idx > 500:
                impingement_history.append(I_accel_step)
            
            if step_idx % 50 == 0 or step_idx == 1:
                print(f"[Core Gap {gap}mm | n0: {n0:.1e}] Step {step_idx}/{steady_state_steps} -> I_imp: {I_accel_step:.2e} A", flush=True)

        # Physics Math for Perveance X-Axis
        v_bohm = np.sqrt(q * params['Te_up'] / m_XE)
        A_inject = (params['rs'] - 0.05) * 1e-3 * 1.0
        I_ion = q * 0.61 * n0 * v_bohm * A_inject
        V_total = params['Vs'] - params['Va']
        perveance = I_ion / (V_total**1.5)
        
        # Mean impingement current at steady state
        final_I_imp = np.mean(impingement_history) if impingement_history else 0.0
        
        result_queue.put(('data', gap, perveance, final_I_imp))
        
    result_queue.put(('done', gap))

# --- MAIN GUI THREAD ---
def run_impingement_benchmark():
    grid_gaps = [0.5, 0.75, 1.0] 
    n0_sweep = np.linspace(1e16, 2.5e17, 15) 
    
    plt.ion()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_title('Impingement Benchmark: Accelerator Grid Current vs. Perveance')
    ax.set_xlabel('Perveance (A/V^{1.5})')
    ax.set_ylabel('Impingement Current (Amperes)')
    ax.set_yscale('log') # Better for seeing small crossover currents
    ax.grid(True, which="both", ls="-", alpha=0.5)
    
    lines_sim = {}
    data_store = {gap: {'perv': [], 'curr': []} for gap in grid_gaps}
    colors = ['tab:blue', 'tab:orange', 'tab:green']
    
    for i, gap in enumerate(grid_gaps):
        lines_sim[gap], = ax.plot([], [], marker='o', markersize=4, color=colors[i], label=f'Gap = {gap} mm')
        
    ax.legend()
    fig.canvas.draw()
    
    manager = mp.Manager()
    result_queue = manager.Queue()
    workers = []
    
    for gap in grid_gaps:
        p = mp.Process(target=worker_sweep, args=(gap, n0_sweep, result_queue))
        p.start()
        workers.append(p)
        
    completed_workers = 0
    while completed_workers < len(grid_gaps):
        try:
            msg = result_queue.get(timeout=0.1)
            if msg[0] == 'done':
                completed_workers += 1
            elif msg[0] == 'data':
                _, gap, perv, curr = msg
                data_store[gap]['perv'].append(perv)
                data_store[gap]['curr'].append(curr)
                
                sort_idx = np.argsort(data_store[gap]['perv'])
                lines_sim[gap].set_data(np.array(data_store[gap]['perv'])[sort_idx], 
                                        np.array(data_store[gap]['curr'])[sort_idx])
                ax.relim()
                ax.autoscale_view()
        except queue.Empty:
            pass
        fig.canvas.draw_idle()
        fig.canvas.flush_events()

    for p in workers: p.join()
    plt.ioff() 
    plt.savefig('impingement_current_benchmark.png', dpi=300)
    print("Benchmark complete. Plot saved as 'impingement_current_benchmark.png'.")
    plt.show() 

if __name__ == "__main__":
    run_impingement_benchmark()