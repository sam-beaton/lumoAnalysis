# lumoAnalysis

## Overview
This repository contains scripts for data analysis using the NeuroDOT toolbox.



## Repository Structure
```
nDotAnalysis/
│ 
├── parcelAnalysisVolumetric.m # runs image reconstruction & saves data in volumetric parcel space
├── parcelAnalysisSurface.m # runs image reconstruction & saves data in parcel space using a surface-based parcellation
├── tileAnalysisParallel.m # runs channel and tile space analysis
│ 
├── +analysisTools/ # modular helper functions
├── +plotting/ # modular plotting functions
├── auxilliary/ # pre-analysis helper functions to format data etc.
│ 
│── README.md


```

## Requirements
- MATLAB (Tested on version R2022b)
- Required Matlab toolboxes: NeuroDOT, NIRFASTer, Mesh2EEG

## License
[MIT License](LICENSE)

