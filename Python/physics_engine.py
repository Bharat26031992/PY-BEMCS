import numpy as np
from scipy.interpolate import RegularGridInterpolator

class DigitalTwinSimulator:
    def __init__(self):
        # Fundamental Constants
        self.dt = 2e-9
        self.q = 1.602e-19
        self.m_XE = 131.293 * 1.6605e-27
        self.kB = 1.380649e-23
        self.eps0 = 8.854e-12
        
        # ARTIFICIAL ELECTRON MASS: 100x lighter than Xenon
        self.m_e = self.m_XE / 100.0 
        
        # --- THERMAL MULTIPHYSICS CONSTANTS ---
        self.macro_weight = 1e5  
        self.alpha_moly = 4.8e-6 
        self.sb_sigma = 5.67e-8  
        self.emissivity = 0.8    
        self.thermal_accel = 1e7 
        
        self.C_cell = 10280 * (0.05e-3)**2 * 1e-3 * 250 
        self.A_cell = 2 * (0.05e-3) * 1e-3 
        
        self.T_screen = 300.0 
        self.T_accel = 300.0 
        
        # Mesh parameters
        self.Lx = 8.0
        self.Ly = 4.0
        self.dx = 0.05
        self.dy = 0.05
        self.nx = int(self.Lx / self.dx) + 1
        self.ny = int(self.Ly / self.dy) + 1
        
        self.x_pts = np.linspace(0, self.Lx, self.nx)
        self.y_pts = np.linspace(0, self.Ly, self.ny)
        self.X, self.Y = np.meshgrid(self.x_pts, self.y_pts)
        
        self.iteration = 0
        self.T_map = np.full((self.ny, self.nx), 300.0)
        self.reset_arrays()

    def reset_arrays(self):
        self.p_x, self.p_y = np.array([]), np.array([])
        self.p_vx, self.p_vy = np.array([]), np.array([])
        self.p_isCEX = np.array([], dtype=bool)
        
        self.e_x, self.e_y = np.array([]), np.array([])
        self.e_vx, self.e_vy = np.array([]), np.array([])
        
        self.V = np.zeros((self.ny, self.nx))
        self.rho = np.zeros((self.ny, self.nx)) 
        self.isBound = np.zeros((self.ny, self.nx), dtype=bool)
        self.V_fixed = np.zeros((self.ny, self.nx))
        self.damage_map = np.zeros((self.ny, self.nx))
        self.Ex, self.Ey = np.zeros((self.ny, self.nx)), np.zeros((self.ny, self.nx))
        self.interp_Ex, self.interp_Ey = None, None

    def build_domain(self, params):
        self.reset_arrays()
        self.iteration = 0

        screen_start = 1.0
        screen_end = screen_start + params['ts']
        accel_start = screen_end + params['gap']
        accel_end = accel_start + params['ta']

        in_screen = (self.X >= screen_start) & (self.X <= screen_end)
        in_accel = (self.X >= accel_start) & (self.X <= accel_end)

        R_screen = params['rs'] + np.maximum(0, self.X - screen_start) * np.tan(np.radians(params['cham_s']))
        R_accel = params['ra'] + np.maximum(0, self.X - accel_start) * np.tan(np.radians(params['cham_a']))

        self.mask_s = in_screen & (self.Y >= R_screen)
        self.isBound[self.mask_s] = True
        self.V_fixed[self.mask_s] = params['Vs']

        self.mask_a = in_accel & (self.Y >= R_accel)
        self.isBound[self.mask_a] = True
        self.V_fixed[self.mask_a] = params['Va']

        self.isBound[:, 0] = True
        self.V_fixed[:, 0] = params['Vs'] + 50
        
        self.T_map[~self.isBound] = 300.0
        self.recalc_poisson(iterations=500)

    def recalc_poisson(self, iterations=500):
        self.V[self.isBound] = self.V_fixed[self.isBound]
        dx_m2 = (self.dx * 1e-3)**2 
        coeff = dx_m2 / self.eps0

        for _ in range(iterations):
            self.V[1:-1, 1:-1] = 0.25 * (
                self.V[2:, 1:-1] + self.V[:-2, 1:-1] + 
                self.V[1:-1, 2:] + self.V[1:-1, :-2] + 
                coeff * self.rho[1:-1, 1:-1]
            )
            self.V[self.isBound] = self.V_fixed[self.isBound]
            self.V[0, :] = self.V[1, :]
            self.V[-1, :] = self.V[-2, :]
            self.V[:, -1] = self.V[:, -2]

        self.Ey, self.Ex = np.gradient(-self.V, self.dy * 1e-3, self.dx * 1e-3)
        self.interp_Ex = RegularGridInterpolator((self.y_pts, self.x_pts), self.Ex, bounds_error=False, fill_value=0)
        self.interp_Ey = RegularGridInterpolator((self.y_pts, self.x_pts), self.Ey, bounds_error=False, fill_value=0)

    def step(self, params):
        if not np.any(self.Ex): return False, np.nan, np.nan, self.T_screen, self.T_accel
        sim_mode = params.get('sim_mode', 'Both')
        self.iteration += 1
        
        # --- A. INJECT PARTICLES ---
        num_inject = 30
        new_y = np.linspace(0.05, params['rs'] - 0.1, num_inject) + (np.random.rand(num_inject) - 0.5) * 0.05
        new_x = np.full(num_inject, 0.1)
        v_bohm = np.sqrt(2 * self.q * 50 / self.m_XE)
        v_spread = np.sqrt(self.q * params['Ti'] / self.m_XE)

        self.p_x = np.concatenate((self.p_x, new_x))
        self.p_y = np.concatenate((self.p_y, new_y))
        self.p_vx = np.concatenate((self.p_vx, np.full(num_inject, v_bohm) + np.random.randn(num_inject) * v_spread))
        self.p_vy = np.concatenate((self.p_vy, np.random.randn(num_inject) * v_spread))
        self.p_isCEX = np.concatenate((self.p_isCEX, np.zeros(num_inject, dtype=bool)))

        # --- DYNAMIC NEUTRALIZER CONTROL ---
        num_e = int(params.get('neut_rate', 30))
        Te_eV = params.get('Te', 5.0)

        if num_e > 0:
            new_ey = np.random.rand(num_e) * self.Ly
            new_ex = np.full(num_e, self.Lx - 0.1)
            v_e_th = np.sqrt(2 * self.q * Te_eV / self.m_e)
            
            self.e_x = np.concatenate((self.e_x, new_ex))
            self.e_y = np.concatenate((self.e_y, new_ey))
            self.e_vx = np.concatenate((self.e_vx, -np.abs(np.random.randn(num_e)) * v_e_th - v_e_th*0.2)) 
            self.e_vy = np.concatenate((self.e_vy, np.random.randn(num_e) * v_e_th))

        # --- B. POISSON & SPACE CHARGE ---
        self.rho.fill(0.0)
        cell_vol = (self.dx * 1e-3) * (self.dy * 1e-3) * 1e-3 
        charge_per_particle = self.q * self.macro_weight
        
        if len(self.p_x) > 0:
            ix_rho = np.clip(np.round(self.p_x / self.dx).astype(int), 1, self.nx - 2)
            iy_rho = np.clip(np.round(self.p_y / self.dy).astype(int), 1, self.ny - 2)
            np.add.at(self.rho, (iy_rho, ix_rho), charge_per_particle / cell_vol)
            
        if len(self.e_x) > 0:
            ix_rho_e = np.clip(np.round(self.e_x / self.dx).astype(int), 1, self.nx - 2)
            iy_rho_e = np.clip(np.round(self.e_y / self.dy).astype(int), 1, self.ny - 2)
            np.add.at(self.rho, (iy_rho_e, ix_rho_e), -charge_per_particle / cell_vol)

        if self.iteration % 2 == 0:
            self.recalc_poisson(iterations=20)

        # --- C. PUSH PARTICLES ---
        pts = np.column_stack((self.p_y, self.p_x))
        Ex_p, Ey_p = self.interp_Ex(pts), self.interp_Ey(pts)
        self.p_vx += (self.q / self.m_XE) * Ex_p * self.dt
        self.p_vy += (self.q / self.m_XE) * Ey_p * self.dt
        self.p_x += self.p_vx * self.dt * 1000
        self.p_y += self.p_vy * self.dt * 1000

        if len(self.e_x) > 0:
            pts_e = np.column_stack((self.e_y, self.e_x))
            Ex_e, Ey_e = self.interp_Ex(pts_e), self.interp_Ey(pts_e)
            self.e_vx += (-self.q / self.m_e) * Ex_e * self.dt
            self.e_vy += (-self.q / self.m_e) * Ey_e * self.dt
            self.e_x += self.e_vx * self.dt * 1000
            self.e_y += self.e_vy * self.dt * 1000

        max_grid_x = 1.0 + params['ts'] + params['gap'] + params['ta']
        post_grid = (~self.p_isCEX) & (self.p_x > max_grid_x)
        current_div = np.percentile(np.abs(np.arctan2(self.p_vy[post_grid], self.p_vx[post_grid])) * 180 / np.pi, 95) if np.sum(post_grid) > 5 else np.nan
        min_pot = np.min(self.V[0, :])

        # --- D. 2D THERMAL MULTIPHYSICS & HIT DETECTION ---
        ix = np.clip(np.round(self.p_x / self.dx).astype(int), 0, self.nx - 1)
        iy = np.clip(np.round(self.p_y / self.dy).astype(int), 0, self.ny - 1)
        hit_grid = self.isBound[iy, ix]
        out_of_bounds = (self.p_x < 0) | (self.p_x > self.Lx) | (self.p_y < 0) | (self.p_y > self.Ly) | np.isnan(self.p_x)

        remeshed = False
        
        if sim_mode in ['Thermal', 'Both'] and np.any(hit_grid):
            v_mag_sq = self.p_vx[hit_grid]**2 + self.p_vy[hit_grid]**2
            E_joules = 0.5 * self.m_XE * v_mag_sq * self.macro_weight
            dT_heat = (E_joules / self.C_cell) * self.thermal_accel
            np.add.at(self.T_map, (iy[hit_grid], ix[hit_grid]), dT_heat)

        if sim_mode in ['Thermal', 'Both']:
            T_bound = self.T_map[self.isBound]
            cooling_factor = (self.emissivity * self.sb_sigma * self.A_cell * self.dt * self.thermal_accel) / self.C_cell
            dT_cool = cooling_factor * (T_bound**4 - 300.0**4)
            
            self.T_map[self.isBound] -= dT_cool
            self.T_map[self.isBound] = np.maximum(self.T_map[self.isBound], 300.0) 
            
            if hasattr(self, 'mask_s') and np.any(self.mask_s): 
                self.T_screen = np.mean(self.T_map[self.mask_s])
            if hasattr(self, 'mask_a') and np.any(self.mask_a): 
                self.T_accel = np.mean(self.T_map[self.mask_a])

            delta_gap = self.alpha_moly * params['gap'] * (self.T_accel - 300) 
            if abs(delta_gap) > 0.05:
                params['gap'] -= delta_gap 
                self.build_domain(params) 
                remeshed = True

        is_erosion_hit = hit_grid & self.p_isCEX
        if sim_mode in ['Erosion', 'Both'] and np.any(is_erosion_hit):
            E_eV = (0.5 * self.m_XE * (self.p_vx[is_erosion_hit]**2 + self.p_vy[is_erosion_hit]**2)) / self.q
            Y_yield = np.where(E_eV > 30, 1.05e-4 * (E_eV - 30)**1.5, 0)
            np.add.at(self.damage_map, (iy[is_erosion_hit], ix[is_erosion_hit]), Y_yield * params['Accel'])

            broken_cells = (self.damage_map > params['Thresh']) & self.isBound
            if np.any(broken_cells):
                self.isBound[broken_cells] = False
                self.damage_map[broken_cells] = 0
                self.recalc_poisson()
                remeshed = True

        # Purge dead Ions
        dead_mask = hit_grid | out_of_bounds
        self.p_x, self.p_y = self.p_x[~dead_mask], self.p_y[~dead_mask]
        self.p_vx, self.p_vy = self.p_vx[~dead_mask], self.p_vy[~dead_mask]
        self.p_isCEX = self.p_isCEX[~dead_mask]
        
        # Purge dead Electrons
        if len(self.e_x) > 0:
            ix_e = np.clip(np.round(self.e_x / self.dx).astype(int), 0, self.nx - 1)
            iy_e = np.clip(np.round(self.e_y / self.dy).astype(int), 0, self.ny - 1)
            hit_grid_e = self.isBound[iy_e, ix_e]
            out_e = (self.e_x < 0) | (self.e_x > self.Lx) | (self.e_y < 0) | (self.e_y > self.Ly) | np.isnan(self.e_x)
            
            dead_e = hit_grid_e | out_e
            self.e_x, self.e_y = self.e_x[~dead_e], self.e_y[~dead_e]
            self.e_vx, self.e_vy = self.e_vx[~dead_e], self.e_vy[~dead_e]

        # --- E. CEX COLLISIONS ---
        primary_mask = (~self.p_isCEX) & (self.p_x > 1.0)
        if np.any(primary_mask):
            v_mag = np.sqrt(self.p_vx[primary_mask]**2 + self.p_vy[primary_mask]**2)
            g = np.maximum(v_mag, 1)
            sigma = ((-0.8821 * np.log(g) + 15.1262)**2) * 1e-20
            prob = 1 - np.exp(-params['n0'] * sigma * g * self.dt)
            collided = np.random.rand(np.sum(primary_mask)) < prob
            
            if np.any(collided):
                c_idx = np.where(primary_mask)[0][collided]
                neut_vth = np.sqrt(2 * self.kB * params['Tn'] / self.m_XE)
                self.p_vx[c_idx] = np.random.randn(len(c_idx)) * neut_vth
                self.p_vy[c_idx] = np.random.randn(len(c_idx)) * neut_vth
                self.p_isCEX[c_idx] = True

        return remeshed, min_pot, current_div, self.T_screen, self.T_accel