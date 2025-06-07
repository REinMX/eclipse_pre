# EXAMPLE_FIELD - Eclipse Simulation Deck

Generated: 2025-05-25 07:10:29

## Project Structure
```
EXAMPLE_FIELD/
├── EXAMPLE_FIELD.DATA          # Main Eclipse data file
├── include/
│   ├── GRID/                # Grid geometry and properties
│   ├── PROPS/               # Fluid and rock properties  
│   ├── SOLUTION/            # Initial conditions
│   └── SCHEDULE/            # Wells and production schedule
├── scripts/                 # Python generation scripts
└── output/                  # Simulation results
```

## Model Parameters
- **Grid**: 50 x 50 x 20 cells
- **Porosity**: 0.25 ± 0.05
- **Permeability**: 150 mD (average)
- **Reservoir depth**: 2000 - 2150 m
- **Initial pressure**: 250 bar @ 2075 m
- **OWC depth**: 2100 m
- **Simulation period**: 15 years

## Wells
- **PROD1**: Producer at (12, 12)
- **PROD2**: Producer at (37, 37)  
- **INJE1**: Water injector at (25, 25)

## Usage
```bash
# Run Eclipse simulation
eclipse EXAMPLE_FIELD.DATA

# Or with other simulators
flow EXAMPLE_FIELD.DATA
```
