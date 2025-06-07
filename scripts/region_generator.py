#!/usr/bin/env python3
"""
Region Generator for Eclipse
Creates FIPNUM, SATNUM, PVTNUM regions
"""

import numpy as np
import xtgeo
from typing import Tuple


def create_region_fields(grid: xtgeo.Grid, output_dir: str = "../include/SOLUTION"):
    """Generate all region definition fields"""
    
    nx, ny, nz = grid.dimensions
    
    # Get coordinates
    _, _, z_coords = grid.get_xyz(asmasked=False)
    if hasattr(z_coords, 'values'):
        z_coords = z_coords.values.flatten()
    else:
        z_coords = z_coords.flatten()
    
    # FIPNUM - Fluid in place regions (depth-based)
    z_min, z_max = np.nanmin(z_coords), np.nanmax(z_coords)
    fipnum = np.ones(grid.nactive, dtype=int)
    
    # Create 3 FIP regions by depth
    region_thickness = (z_max - z_min) / 3
    for i, depth in enumerate(z_coords):
        region = int((depth - z_min) / region_thickness) + 1
        fipnum[i] = min(region, 3)
    
    # SATNUM - Single saturation function region (since we only have one SWOF/SGOF table)
    satnum = np.ones(grid.nactive, dtype=int)
    
    # PVTNUM - Single PVT region
    pvtnum = np.ones(grid.nactive, dtype=int)
    
    # Export regions
    regions = {
        'FIPNUM': fipnum,
        'SATNUM': satnum, 
        'PVTNUM': pvtnum
    }
    
    for region_name, region_values in regions.items():
        filename = f"{output_dir}/{region_name}.INC"
        
        with open(filename, 'w') as f:
            f.write(f"{region_name}\n")
            
            ncol = 10  # 10 integers per line
            for i in range(0, len(region_values), ncol):
                line_values = region_values[i:i+ncol]
                formatted_line = " ".join(f"{val:3d}" for val in line_values)
                f.write(f" {formatted_line}\n")
            
            f.write("/\n\n")
    
    print(f"Generated region files in {output_dir}/:")
    for region_name in regions.keys():
        print(f"  - {region_name}.INC")
    
    return regions


if __name__ == "__main__":
    # Test with existing grid
    try:
        grid = xtgeo.grid_from_file("../include/GRID/GRID.GRDECL")
        create_region_fields(grid)
    except:
        print("No grid found - run eclipse_deck_builder.py first")

