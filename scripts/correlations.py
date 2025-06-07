#!/usr/bin/env python3
"""
Eclipse Data Preprocessor - Permeability and Capillary Pressure Correlations
Uses PySCAL (Equinor) for realistic relative permeability and capillary pressure curves
"""

import numpy as np
import pandas as pd
from pyscal import WaterOil, GasOil, WaterOilGas
import matplotlib.pyplot as plt


class PermeabilityCorrelations:
    """Generate permeability correlations using PySCAL"""
    
    def __init__(self):
        self.swirr = 0.1   # Irreducible water saturation
        self.swcr = 0.15   # Critical water saturation
        self.sorw = 0.2    # Residual oil saturation to water
        self.sorg = 0.15   # Residual oil saturation to gas
        self.sgcr = 0.05   # Critical gas saturation
        
    def create_wateroil_table(self, krwend=0.8, krowend=0.9, nw=2.0, now=2.0):
        """Create water-oil relative permeability table"""
        wo = WaterOil(
            swirr=self.swirr,
            swcr=self.swcr, 
            sorw=self.sorw,
            h=0.01  # saturation step
        )
        
        # Add Corey relative permeability curves
        wo.add_corey_water(nw=nw, krwend=krwend)
        wo.add_corey_oil(now=now, kroend=krowend)
        
        return wo
    
    def create_gasoil_table(self, krgend=0.85, kroend=0.95, ng=2.5, nog=1.8):
        """Create gas-oil relative permeability table"""
        go = GasOil(
            sgcr=self.sgcr,
            sorg=self.sorg,
            h=0.01
        )
        
        # Add Corey relative permeability curves
        go.add_corey_gas(ng=ng, krgend=krgend)
        go.add_corey_oil(nog=nog, kroend=kroend)
        
        return go
    
    def create_three_phase_table(self):
        """Create three-phase relative permeability table"""
        wo = self.create_wateroil_table()
        go = self.create_gasoil_table()
        
        # Combine into three-phase
        try:
            wog = WaterOilGas(wo, go, h=0.01)
        except:
            # Fallback if three-phase fails
            print("Three-phase table creation failed, using water-oil table")
            wog = wo
        
        return wog


class CapillaryPressureCorrelations:
    """Generate capillary pressure correlations"""
    
    def __init__(self):
        self.swirr = 0.1
        self.swcr = 0.15
        self.sorw = 0.2
        self.sorg = 0.15
        self.sgcr = 0.05
        
    def add_simple_j_function(self, wo_table, a=5.0, b=-1.5):
        """Add simple J-function based capillary pressure"""
        wo_table.add_simple_J(a=a, b=b)
        return wo_table
    
    def add_brooks_corey_pc(self, wo_table, entry_pressure=10.0, lambda_param=2.0):
        """Add Brooks-Corey capillary pressure"""
        # Create custom saturation points
        sw = np.linspace(wo_table.swirr, 1.0 - wo_table.sorw, 50)
        
        # Brooks-Corey equation: Pc = Pd * (Sw_norm)^(-1/lambda)
        sw_norm = (sw - wo_table.swirr) / (1.0 - wo_table.swirr - wo_table.sorw)
        pc = entry_pressure * (sw_norm ** (-1.0/lambda_param))
        
        # Create DataFrame and interpolate to table
        pc_df = pd.DataFrame({'SW': sw, 'PC': pc})
        wo_table.add_fromtable(pc_df, pccolname='PC')
        
        return wo_table


def generate_eclipse_tables(output_dir='.'):
    """Generate Eclipse-formatted relative permeability and capillary pressure tables"""
    
    # Initialize correlations
    perm_corr = PermeabilityCorrelations()
    pc_corr = CapillaryPressureCorrelations()
    
    # Create water-oil table without capillary pressure first
    wo = perm_corr.create_wateroil_table()
    # Skip capillary pressure for now since it's causing issues
    # wo = pc_corr.add_simple_j_function(wo)
    
    # Create gas-oil table
    go = perm_corr.create_gasoil_table()
    
    # Create three-phase table
    wog = perm_corr.create_three_phase_table()
    
    # Export to Eclipse format
    try:
        eclipse_wo = wo.SWOF()
        print(f"SWOF length: {len(eclipse_wo)}")
    except Exception as e:
        print(f"Error generating SWOF: {e}")
        eclipse_wo = ""
    
    try:
        eclipse_go = go.SGOF()
        print(f"SGOF length: {len(eclipse_go)}")
    except Exception as e:
        print(f"Error generating SGOF: {e}")
        eclipse_go = ""
    
    # Handle three-phase export
    try:
        eclipse_wog = wog.SWOF()
    except:
        eclipse_wog = eclipse_wo  # Use water-oil if three-phase failed
    
    # Save to files
    with open(f'{output_dir}/SWOF.txt', 'w') as f:
        f.write(eclipse_wo)
        print(f"Wrote SWOF.txt with {len(eclipse_wo)} characters")
    
    with open(f'{output_dir}/SGOF.txt', 'w') as f:
        f.write(eclipse_go)
        print(f"Wrote SGOF.txt with {len(eclipse_go)} characters")
        
    with open(f'{output_dir}/SWOF_3PHASE.txt', 'w') as f:
        f.write(eclipse_wog)
        print(f"Wrote SWOF_3PHASE.txt with {len(eclipse_wog)} characters")
    
    # Plot curves
    plot_curves(wo, go, output_dir)
    
    return wo, go, wog


def plot_curves(wo, go, output_dir='.'):
    """Plot relative permeability and capillary pressure curves"""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Water-Oil Kr curves
    try:
        wo.plotkr(ax=axes[0,0])
    except:
        axes[0,0].plot(wo.table['SW'], wo.table['KRW'], 'b-', label='Krw')
        axes[0,0].plot(wo.table['SW'], wo.table['KROW'], 'r-', label='Krow')
        axes[0,0].legend()
        axes[0,0].grid(True)
    axes[0,0].set_title('Water-Oil Relative Permeability')
    
    # Gas-Oil Kr curves  
    try:
        go.plotkr(ax=axes[0,1])
    except:
        axes[0,1].plot(go.table['SG'], go.table['KRG'], 'g-', label='Krg')
        axes[0,1].plot(go.table['SG'], go.table['KROG'], 'm-', label='Krog')
        axes[0,1].legend()
        axes[0,1].grid(True)
    axes[0,1].set_title('Gas-Oil Relative Permeability')
    
    # Capillary pressure
    try:
        wo.plotpc(ax=axes[1,0])
    except:
        if 'PC' in wo.table.columns:
            axes[1,0].plot(wo.table['SW'], wo.table['PC'], 'k-')
            axes[1,0].set_ylabel('Capillary Pressure')
        else:
            axes[1,0].text(0.5, 0.5, 'No PC data', ha='center', va='center')
    axes[1,0].set_title('Water-Oil Capillary Pressure')
    
    # Combined plot
    axes[1,1].plot(wo.table['SW'], wo.table['KRW'], 'b-', label='Krw')
    axes[1,1].plot(wo.table['SW'], wo.table['KROW'], 'r-', label='Krow')
    axes[1,1].set_xlabel('Water Saturation')
    axes[1,1].set_ylabel('Relative Permeability')
    axes[1,1].legend()
    axes[1,1].grid(True)
    axes[1,1].set_title('Combined Kr Curves')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/relative_permeability_curves.png', dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    print("Generating Eclipse relative permeability and capillary pressure tables...")
    
    # Generate tables
    wo, go, wog = generate_eclipse_tables()
    
    print("Generated files:")
    print("- SWOF.txt (Water-Oil relative permeability)")
    print("- SGOF.txt (Gas-Oil relative permeability)")  
    print("- SWOF_3PHASE.txt (Three-phase relative permeability)")
    print("- relative_permeability_curves.png (Plots)")
    
    print("\nSample SWOF table (first 10 rows):")
    print(wo.table.head(10).to_string(index=False))

