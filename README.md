# preprocessingLumo

## Overview
This repository contains scripts for DOT (LUMO) data processing using selected methods and user-defined parameters.

- **MATLAB (`pruningComparisonsMatlab/`)**: Processes raw data and generates output tables.
- **R (`pruningComparisonsMatlab/`)**: Runs multilevel models (MLMs) on processed data. Also includes a parameter choice script.

## Repository Structure
```
pruningComparisons/
│ 
│── pruningComparisonsMatlab/
│ 
├── prepRunAll.m #master script - change parameters here
│ 
├── +prepTools/ # modular helper functions
│   ├── functions to go here ....
│ 
│── README.md


```

## Requirements
- MATLAB (Tested on version R2022b)
- Required Matlab toolboxes: qt-nirs, Homer2, DOT-HUB_toolbox

## License
[MIT License](LICENSE)