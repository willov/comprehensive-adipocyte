# A comprehensive dynamic model of the adipocyte

These are the scripts that are used in the publication titled *"A comprehensive mechanistic model of adipocyte signaling with layers of confidence"*.

The models are simulated using [IQM tools](https://iqmtools.intiquan.com/), a toolbox for MATLAB. The toolbox requires a working C-compiler on the system. 

**To generate the figures used in the paper**, run the script in the root folder named `PlotFigures_combined`. This will generate the supplementary Figures 1-6, and the main Figures 3-6. To generate the phosphoproteome-wide figures, use the `PlotFigures_phosphoprotome` function.

Note that running all simulation can take quite a bit of time, therefore we have added an optional, slightly faster, version of the simulations with a bit lower resolution for the dose-response simulations.

The model equations are available in the files in the `Models` folder, and the data used is available in the `Data` folder. 

The code used in the automatic model expansion are available in `phosphoproteome` folder, starting from the `ExpandModel` script. To rerun the expansion (not needed to generate the figures from the paper), the [MEIGO64 toolbox](https://github.com/gingproc-IIM-CSIC/MEIGO64) is required (not included in this repo).