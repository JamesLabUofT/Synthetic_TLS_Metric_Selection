# Code and Data for the "Improving Selection of Forest Fuel Metrics Using Synthetic Lidar Point Clouds" paper

Source code can be found in the `src` folder:

-   Getting the parameter space from real plots: `src/ParameterSpace.R`

-   Creating Scenarios for each parameter combination: `src/ScenarioCreator.R`

-   Running simulations (generating a point cloud for each scenario, computing the metrics): `src/RunSimulations.R`

-   Running statistical analyses on the results: `src/Analysis.R`

-   Validating the methods by comparing metric outputs of real TLS scans vs Synthetic equivalents: `src/ValidateMethod.R`

-   `src/CreateSyntheticCloud.R` , `src/FuelMetrics.R` ,`src/GetObjects.R` are sourced by the above files and include helper functions.

The objects for the object library (as .laz files) are distributed in the `data/Snag`, `data/Tree`, and `data/Understory` directories.

Simulation Outputs (metric outputs of simulated stands) are found in `data/SimulationOutputsRaw`
