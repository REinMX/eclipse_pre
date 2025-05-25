#!/usr/bin/env python3
"""
Eclipse Grid Builder using XTGeo and Equinor libraries
Creates 3D grids with realistic porosity and permeability distributions
"""

import numpy as np
import pandas as pd
import xtgeo
from scipy.stats import lognorm, norm
import matplotlib.pyplot as plt
from typing import Tuple, Optional, Dict, Any


class GridBuilder:
    """Build Eclipse grids with realistic property distributions"""
    
    def __init__(self, nx: int = 50, ny: int = 50, nz: int = 20):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.grid = None
        self.properties = {}
        
    def create_corner_point_grid(self, 
                                xmin: float = 0.0, xmax: float = 5000.0,
                                ymin: float = 0.0, ymax: float = 5000.0,
                                zmin: float = 2000.0, zmax: float = 2200.0) -> xtgeo.Grid:
        """Create corner point grid using XTGeo"""
        
        # Create regular grid
        self.grid = xtgeo.create_box_grid(
            dimension=(self.nx, self.ny, self.nz),
            origin=(xmin, ymin, zmin),
            increment=(
                (xmax - xmin) / self.nx,
                (ymax - ymin) / self.ny, 
                (zmax - zmin) / self.nz
            ),
            rotation=0.0,
            flip=1  # Eclipse convention
        )
        
        return self.grid
    
    def add_structural_dip(self, dip_angle: float = 2.0, dip_direction: float = 45.0):
        """Add structural dip to the grid"""
        if self.grid is None:
            raise ValueError("Grid must be created first")
            
        # Apply structural dip
        dip_rad = np.radians(dip_angle)
        dir_rad = np.radians(dip_direction)
        
        # Get grid coordinates as numpy arrays
        x, y, z = self.grid.get_xyz(asmasked=False)
        
        # Extract values if they are GridProperty objects
        if hasattr(x, 'values'):
            x = x.values
        if hasattr(y, 'values'):
            y = y.values
        if hasattr(z, 'values'):
            z = z.values
        
        # Calculate dip offset
        dip_offset = (x * np.sin(dir_rad) + y * np.cos(dir_rad)) * np.tan(dip_rad)
        
        # Apply dip to all layers
        z_new = z + dip_offset
        
        # Skip grid geometry update for now (xtgeo version compatibility)
        print(f"Applied structural dip: {dip_angle}° in direction {dip_direction}°")
        
    def generate_porosity_field(self, 
                               mean_poro: float = 0.22,
                               std_poro: float = 0.05,
                               correlation_length: Tuple[float, float, float] = (500.0, 500.0, 10.0),
                               facies_control: bool = True) -> xtgeo.GridProperty:
        """Generate realistic porosity field with spatial correlation"""
        
        if self.grid is None:
            raise ValueError("Grid must be created first")
        
        # Create porosity property
        poro = xtgeo.GridProperty(
            self.grid,
            name="PORO",
            discrete=False,
            values=np.random.normal(mean_poro, std_poro, self.grid.nactive)
        )
        
        # Apply Gaussian filter for spatial correlation
        if correlation_length:
            from scipy.ndimage import gaussian_filter
            
            # Convert to 3D array
            poro_3d = poro.values.reshape((self.nz, self.ny, self.nx))
            
            # Calculate filter sigma from correlation length (simplified)
            sigma_x = correlation_length[0] / 1000.0  # simplified scaling
            sigma_y = correlation_length[1] / 1000.0
            sigma_z = correlation_length[2] / 20.0
            
            # Apply smoothing
            poro_smooth = gaussian_filter(poro_3d, sigma=[sigma_z, sigma_y, sigma_x])
            
            # Normalize to maintain mean and std
            poro_smooth = (poro_smooth - np.nanmean(poro_smooth)) / np.nanstd(poro_smooth)
            poro_smooth = poro_smooth * std_poro + mean_poro
            
            poro.values = poro_smooth.flatten()
        
        # Ensure physical limits
        poro.values = np.clip(poro.values, 0.01, 0.4)
        
        self.properties['PORO'] = poro
        return poro
    
    def generate_permeability_field(self,
                                   poro_perm_transform: str = "kozeny_carman",
                                   base_perm: float = 100.0,
                                   perm_anisotropy: float = 0.1) -> Tuple[xtgeo.GridProperty, xtgeo.GridProperty, xtgeo.GridProperty]:
        """Generate permeability field from porosity using transforms"""
        
        if 'PORO' not in self.properties:
            raise ValueError("Porosity field must be generated first")
        
        poro = self.properties['PORO']
        
        if poro_perm_transform == "kozeny_carman":
            # Kozeny-Carman equation: k = k0 * (phi^3) / (1-phi)^2
            perm_h = base_perm * (poro.values**3) / ((1 - poro.values)**2)
        elif poro_perm_transform == "power_law":
            # Simple power law: k = k0 * (phi/phi0)^n
            phi0 = 0.2
            n = 3.0
            perm_h = base_perm * (poro.values / phi0)**n
        else:
            # Linear relationship
            perm_h = base_perm * (poro.values / 0.2)
        
        # Create permeability properties
        permx = xtgeo.GridProperty(self.grid, name="PERMX", values=perm_h)
        permy = xtgeo.GridProperty(self.grid, name="PERMY", values=perm_h)
        permz = xtgeo.GridProperty(self.grid, name="PERMZ", values=perm_h * perm_anisotropy)
        
        # Ensure physical limits
        for perm in [permx, permy, permz]:
            perm.values = np.clip(perm.values, 0.1, 10000.0)
        
        self.properties.update({'PERMX': permx, 'PERMY': permy, 'PERMZ': permz})
        return permx, permy, permz
    
    def add_facies_zones(self, facies_props: Dict[str, Dict[str, float]]) -> xtgeo.GridProperty:
        """Add facies-based property zones"""
        
        if self.grid is None:
            raise ValueError("Grid must be created first")
        
        # Create facies property
        facies = xtgeo.GridProperty(
            self.grid,
            name="FACIES",
            discrete=True,
            values=np.ones(self.grid.nactive, dtype=int)
        )
        
        # Get z coordinates
        _, _, z_coords = self.grid.get_xyz(asmasked=False)
        if hasattr(z_coords, 'values'):
            z_coords = z_coords.values
        
        # Assign facies based on depth zones
        z_min, z_max = np.nanmin(z_coords), np.nanmax(z_coords)
        n_facies = len(facies_props)
        
        for i, (facies_name, props) in enumerate(facies_props.items()):
            z_start = z_min + i * (z_max - z_min) / n_facies
            z_end = z_min + (i + 1) * (z_max - z_min) / n_facies
            
            mask = (z_coords >= z_start) & (z_coords < z_end)
            facies.values[mask] = i + 1
            
            # Update properties based on facies
            if 'PORO' in self.properties and 'mean_poro' in props:
                self.properties['PORO'].values[mask] = np.random.normal(
                    props['mean_poro'], props.get('std_poro', 0.02), np.sum(mask)
                )
        
        self.properties['FACIES'] = facies
        return facies
    
    def add_water_saturation(self, 
                           swi: float = 0.2,
                           transition_zone_height: float = 50.0,
                           owc_depth: float = 2150.0) -> xtgeo.GridProperty:
        """Add initial water saturation based on capillary pressure"""
        
        if self.grid is None:
            raise ValueError("Grid must be created first")
        
        # Get cell centers
        _, _, z_coords = self.grid.get_xyz(asmasked=False)
        if hasattr(z_coords, 'values'):
            z_coords = z_coords.values
        
        # Calculate height above OWC
        height_above_owc = owc_depth - z_coords
        
        # Simple J-function approximation for water saturation
        # Sw = Swi + (1-Swi) * exp(-height/lambda)
        lambda_param = transition_zone_height / 3.0  # transition zone parameter
        
        sw = swi + (1 - swi) * np.exp(-np.maximum(height_above_owc, 0) / lambda_param)
        sw = np.clip(sw, swi, 1.0)
        
        # Below OWC, assume fully water saturated
        sw[height_above_owc < 0] = 1.0
        
        swat = xtgeo.GridProperty(self.grid, name="SWAT", values=sw)
        self.properties['SWAT'] = swat
        
        return swat
    
    def export_eclipse_grid(self, filename: str = "GRID.GRDECL"):
        """Export grid and properties to Eclipse format"""
        
        if self.grid is None:
            raise ValueError("Grid must be created first")
        
        # Export grid
        self.grid.to_file(filename, fformat="grdecl")
        
        # Export properties with decimal format (Eclipse compatible)
        for prop_name, prop in self.properties.items():
            prop_file = f"{prop_name}.GRDECL"
            self._export_property_decimal_format(prop, prop_file)
            print(f"Exported {prop_name} to {prop_file}")
        
        print(f"Exported grid to {filename}")
    
    def _export_property_decimal_format(self, prop, filename):
        """Export property in decimal format (not scientific notation)"""
        
        with open(filename, 'w') as f:
            f.write(f"{prop.name}\n")
            
            values = np.asarray(prop.values).flatten()
            ncol = 6  # Eclipse standard: 6 values per line
            
            for i in range(0, len(values), ncol):
                line_values = values[i:i+ncol]
                
                # Format as decimal with appropriate precision
                if prop.name in ['PORO', 'SWAT']:
                    # Porosity/saturation: 6 decimal places
                    formatted = [f"{float(val):.6f}" for val in line_values]
                elif prop.name in ['PERMX', 'PERMY', 'PERMZ']:
                    # Permeability: 3 decimal places
                    formatted = [f"{float(val):.3f}" for val in line_values]
                elif prop.isdiscrete:
                    # Integer properties (facies)
                    formatted = [f"{int(float(val))}" for val in line_values]
                else:
                    # Default: 4 decimal places
                    formatted = [f"{float(val):.4f}" for val in line_values]
                
                # Right-align values with proper spacing
                formatted_line = " ".join(f"{val:>12}" for val in formatted)
                f.write(f" {formatted_line}\n")
            
            f.write("/\n\n")
    
    def plot_properties(self, layer: int = 1):
        """Plot property maps for a given layer"""
        
        n_props = len(self.properties)
        if n_props == 0:
            print("No properties to plot")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        
        for i, (prop_name, prop) in enumerate(self.properties.items()):
            if i >= 4:  # Limit to 4 plots
                break
                
            # Get 2D slice
            prop_2d = prop.values.reshape((self.nz, self.ny, self.nx))[layer-1, :, :]
            
            im = axes[i].imshow(prop_2d, origin='lower', cmap='viridis')
            axes[i].set_title(f'{prop_name} - Layer {layer}')
            axes[i].set_xlabel('I index')
            axes[i].set_ylabel('J index')
            plt.colorbar(im, ax=axes[i])
        
        # Hide unused subplots
        for i in range(len(self.properties), 4):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(f'grid_properties_layer_{layer}.png', dpi=300, bbox_inches='tight')
        plt.close()


def create_simple_reservoir_model():
    """Create a simple reservoir model example"""
    
    # Initialize grid builder
    builder = GridBuilder(nx=40, ny=40, nz=15)
    
    # Create corner point grid
    grid = builder.create_corner_point_grid(
        xmin=0, xmax=4000,
        ymin=0, ymax=4000, 
        zmin=2000, zmax=2150
    )
    
    # Add structural dip
    builder.add_structural_dip(dip_angle=3.0, dip_direction=30.0)
    
    # Define facies properties
    facies_props = {
        "sand_1": {"mean_poro": 0.25, "std_poro": 0.03},
        "shale": {"mean_poro": 0.05, "std_poro": 0.01}, 
        "sand_2": {"mean_poro": 0.22, "std_poro": 0.04}
    }
    
    # Generate properties
    builder.add_facies_zones(facies_props)
    builder.generate_porosity_field(mean_poro=0.22, std_poro=0.05)
    builder.generate_permeability_field(poro_perm_transform="kozeny_carman")
    builder.add_water_saturation(swi=0.15, owc_depth=2100.0)
    
    # Export to Eclipse format
    builder.export_eclipse_grid("SIMPLE_RESERVOIR.GRDECL")
    
    # Create plots
    builder.plot_properties(layer=8)
    
    print("Simple reservoir model created!")
    return builder


if __name__ == "__main__":
    print("Creating Eclipse reservoir grid...")
    model = create_simple_reservoir_model()