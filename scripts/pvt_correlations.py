#!/usr/bin/env python3
"""
PVT Correlations for Eclipse
Implements industry-standard black oil correlations (Standing, Glaso, Vazquez-Beggs)
"""

import numpy as np
import pandas as pd
from typing import Tuple, Optional, Dict, Any


class BlackOilPVT:
    """Black oil PVT correlations for Eclipse"""
    
    def __init__(self):
        # Default fluid properties
        self.api_gravity = 35.0        # API gravity
        self.gas_sg = 0.65            # Gas specific gravity (air=1)
        self.temp_res = 80.0          # Reservoir temperature (°C)
        self.temp_std = 15.6          # Standard temperature (°C)
        self.press_std = 1.01325      # Standard pressure (bar)
        
        # Pressure range for tables
        self.p_min = 1.0
        self.p_max = 400.0
        self.n_points = 20
        
        # Water properties
        self.water_comp = 4.5e-5      # Water compressibility (1/bar)
        self.water_visc = 0.5         # Water viscosity (cP)
        self.water_fvf_ref = 1.03     # Water FVF at reference pressure
        
    def set_fluid_properties(self, api_gravity: float, gas_sg: float, 
                           temp_res: float, **kwargs):
        """Set fluid properties"""
        self.api_gravity = api_gravity
        self.gas_sg = gas_sg
        self.temp_res = temp_res
        
        # Update optional properties
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
    
    def oil_density_sg(self) -> float:
        """Calculate oil specific gravity from API gravity"""
        return 141.5 / (131.5 + self.api_gravity)
    
    def standing_rs_correlation(self, pressure: np.ndarray) -> np.ndarray:
        """Standing correlation for solution gas-oil ratio (Rs)"""
        # Convert temperature to Fahrenheit for correlation
        temp_f = self.temp_res * 9/5 + 32
        
        # Standing correlation: Rs = yg * [(p/18.2 + 1.4) * 10^(0.0125*API - 0.00091*T)]^1.2048
        factor = (pressure / 18.2 + 1.4) * (10**(0.0125 * self.api_gravity - 0.00091 * temp_f))
        rs = self.gas_sg * (factor**1.2048)
        
        return np.clip(rs, 0, 2000)  # Reasonable limits
    
    def standing_bo_correlation(self, pressure: np.ndarray, rs: np.ndarray) -> np.ndarray:
        """Standing correlation for oil formation volume factor (Bo)"""
        # Convert temperature to Fahrenheit
        temp_f = self.temp_res * 9/5 + 32
        
        # Standing correlation
        term1 = rs * (self.gas_sg / self.oil_density_sg())**0.5
        term2 = 1.25 * temp_f
        
        bo = 0.9759 + 0.000120 * (term1 + term2)**1.2
        
        return np.clip(bo, 1.0, 3.0)  # Physical limits
    
    def beggs_robinson_oil_viscosity(self, pressure: np.ndarray, 
                                   rs: np.ndarray) -> np.ndarray:
        """Beggs-Robinson correlation for oil viscosity"""
        # Convert temperature to Fahrenheit
        temp_f = self.temp_res * 9/5 + 32
        
        # Dead oil viscosity
        z = 3.0324 - 0.02023 * self.api_gravity
        y = 10**z
        x = y * (temp_f**(-1.163))
        mu_od = 10**x - 1
        
        # Live oil viscosity
        a = 10.715 * (rs + 100)**(-0.515)
        b = 5.44 * (rs + 150)**(-0.338)
        mu_o = a * (mu_od**b)
        
        return np.clip(mu_o, 0.1, 100.0)  # Reasonable viscosity range
    
    def standing_bg_correlation(self, pressure: np.ndarray) -> np.ndarray:
        """Gas formation volume factor using real gas law"""
        # Convert temperature to Rankine
        temp_r = (self.temp_res + 273.15) * 9/5
        
        # Z-factor approximation (simplified)
        z_factor = 0.95 - 0.1 * (pressure / 100.0)**0.5
        z_factor = np.clip(z_factor, 0.7, 1.1)
        
        # Bg = 0.02829 * Z * T / P (reservoir bbl/scf)
        bg = 0.02829 * z_factor * temp_r / pressure
        
        return bg
    
    def lee_gas_viscosity(self, pressure: np.ndarray) -> np.ndarray:
        """Lee correlation for gas viscosity"""
        # Convert temperature to Rankine
        temp_r = (self.temp_res + 273.15) * 9/5
        
        # Gas density
        mg = 28.97 * self.gas_sg  # Molecular weight
        rho_g = pressure * mg / (10.732 * temp_r)  # lb/ft3
        
        # Lee correlation
        k = ((9.4 + 0.02 * mg) * (temp_r**1.5)) / (209 + 19 * mg + temp_r)
        x = 3.5 + 986/temp_r + 0.01 * mg
        y = 2.4 - 0.2 * x
        
        mu_g = k * np.exp(x * (rho_g/62.4)**y) / 10000  # Convert to cP
        
        return np.clip(mu_g, 0.01, 0.1)
    
    def generate_pvdo_table(self) -> str:
        """Generate PVDO table for Eclipse"""
        pressures = np.linspace(self.p_min, self.p_max, self.n_points)
        
        # Calculate PVT properties
        rs = self.standing_rs_correlation(pressures)
        bo = self.standing_bo_correlation(pressures, rs)
        mu_o = self.beggs_robinson_oil_viscosity(pressures, rs)
        
        # Create PVDO table
        pvdo_lines = ["PVDO"]
        for i in range(len(pressures)):
            line = f"{pressures[i]:.1f} {bo[i]:.4f} {mu_o[i]:.4f}"
            pvdo_lines.append(line)
        pvdo_lines.append("/")
        
        return "\n".join(pvdo_lines)
    
    def generate_pvdg_table(self) -> str:
        """Generate PVDG table for Eclipse"""
        pressures = np.linspace(self.p_min, self.p_max, self.n_points)
        
        # Calculate gas properties
        bg = self.standing_bg_correlation(pressures)
        mu_g = self.lee_gas_viscosity(pressures)
        
        # Create PVDG table
        pvdg_lines = ["PVDG"]
        for i in range(len(pressures)):
            line = f"{pressures[i]:.1f} {bg[i]:.6f} {mu_g[i]:.5f}"
            pvdg_lines.append(line)
        pvdg_lines.append("/")
        
        return "\n".join(pvdg_lines)
    
    def generate_pvtw_table(self, ref_pressure: float = 250.0) -> str:
        """Generate PVTW table for Eclipse"""
        # PVTW: Pref Bw Cw mu_w Cv
        pvtw = f"PVTW\n{ref_pressure:.1f} {self.water_fvf_ref:.3f} {self.water_comp:.2e} {self.water_visc:.2f} 0.0 /"
        
        return pvtw
    
    def generate_density_table(self) -> str:
        """Generate DENSITY table for Eclipse"""
        # Calculate densities at standard conditions
        oil_density = self.oil_density_sg() * 1000  # kg/m3
        water_density = 1025.0  # kg/m3
        gas_density = self.gas_sg * 1.225  # kg/m3 (air at SC = 1.225)
        
        density = f"DENSITY\n{oil_density:.1f} {water_density:.1f} {gas_density:.3f} /"
        
        return density
    
    def generate_rock_table(self, ref_pressure: float = 250.0, 
                          rock_comp: float = 4.5e-5) -> str:
        """Generate ROCK compressibility table"""
        rock = f"ROCK\n{ref_pressure:.1f} {rock_comp:.2e} /"
        
        return rock


class CompositionalPVT:
    """Advanced compositional PVT for complex fluids"""
    
    def __init__(self):
        self.components = {
            'N2': {'mw': 28.014, 'tc': -146.9, 'pc': 33.9},
            'CO2': {'mw': 44.010, 'tc': 31.0, 'pc': 73.8},
            'C1': {'mw': 16.043, 'tc': -82.6, 'pc': 46.0},
            'C2': {'mw': 30.070, 'tc': 32.3, 'pc': 48.8},
            'C3': {'mw': 44.097, 'tc': 96.7, 'pc': 42.5},
            'C4': {'mw': 58.123, 'tc': 152.0, 'pc': 38.0},
            'C5': {'mw': 72.150, 'tc': 196.6, 'pc': 33.7},
            'C6': {'mw': 86.177, 'tc': 234.5, 'pc': 30.2},
            'C7+': {'mw': 142.0, 'tc': 350.0, 'pc': 24.0}
        }
    
    def peng_robinson_eos(self, pressure: float, temperature: float, 
                         composition: Dict[str, float]) -> Dict[str, float]:
        """Peng-Robinson equation of state (simplified implementation)"""
        # This would be a full EOS implementation
        # For now, return placeholder
        return {
            'z_factor': 0.85,
            'density': 0.7,
            'fugacity': pressure * 0.9
        }


def generate_comprehensive_pvt(fluid_config: Dict[str, Any], 
                             output_dir: str = "../include/PROPS") -> Dict[str, str]:
    """Generate comprehensive PVT tables for Eclipse"""
    
    # Initialize PVT calculator
    pvt = BlackOilPVT()
    pvt.set_fluid_properties(**fluid_config)
    
    # Generate all PVT tables
    pvt_tables = {
        'PVDO': pvt.generate_pvdo_table(),
        'PVDG': pvt.generate_pvdg_table(), 
        'PVTW': pvt.generate_pvtw_table(),
        'DENSITY': pvt.generate_density_table(),
        'ROCK': pvt.generate_rock_table()
    }
    
    # Save to files
    for table_name, table_content in pvt_tables.items():
        filename = f"{output_dir}/{table_name}.INC"
        with open(filename, 'w') as f:
            f.write(f"-- {table_name} table generated by PVT correlations\n")
            f.write(f"-- API gravity: {pvt.api_gravity}°\n")
            f.write(f"-- Gas SG: {pvt.gas_sg}\n")
            f.write(f"-- Temperature: {pvt.temp_res}°C\n\n")
            f.write(table_content)
            f.write("\n\n")
    
    print(f"Generated PVT tables in {output_dir}/:")
    for table_name in pvt_tables.keys():
        print(f"  - {table_name}.INC")
    
    return pvt_tables


def main():
    """Generate example PVT tables"""
    
    # Fluid configuration
    fluid_config = {
        'api_gravity': 32.0,      # Medium crude
        'gas_sg': 0.68,           # Typical gas
        'temp_res': 85.0,         # Reservoir temperature
        'p_max': 350.0,           # Maximum pressure
        'n_points': 25            # Table resolution
    }
    
    # Generate PVT tables
    pvt_tables = generate_comprehensive_pvt(fluid_config)
    
    print("\nSample PVDO table (first 5 lines):")
    pvdo_lines = pvt_tables['PVDO'].split('\n')
    for line in pvdo_lines[:6]:
        print(line)


if __name__ == "__main__":
    main()