# Cloth-Fit Project TODO

This document outlines the current tasks for the Cloth-Fit project, focused on computing existing garment-data.

## CLI Application (`cloth_fit_cli`) - COMPLETED ✅

*   [x] Queue jobs with sqlite and oban.
*   [x] Implement `cloth-fit` command with full input validation.
*   [x] Add CLI commands for listing assets, managing queues, viewing results.
*   [x] Implement error handling and user feedback.
*   [x] Implement configuration management for CLI.

## Computing Existing Garment Data (`garment-data`)

### Existing Garments
*   [ ] Test and validate `jumpsuit_dense` garment simulation.
*   [ ] Test and validate `LCL_Skirt_DressEvening_003` garment simulation.
*   [ ] Test and validate `Puffer_dense` garment simulation.

### Existing Avatars
*   [ ] Test and validate `FoxGirl` avatar compatibility.
*   [ ] Test and validate `Goblin` avatar compatibility.
*   [ ] Test and validate `T-rex` avatar compatibility.

### Existing Simulation Configurations
*   [ ] Process and validate `foxgirl_skirt` setup configuration.
*   [ ] Process and validate `Goblin_Jacket` setup configuration.
*   [ ] Process and validate `Goblin_Jumpsuit` setup configuration.
*   [ ] Process and validate `Trex_Jacket` setup configuration.

### Data Processing Tasks
*   [ ] Verify all existing garment mesh files are properly formatted.
*   [ ] Validate all existing avatar mesh files and skeletons.
*   [ ] Ensure all setup.json files contain valid configuration data.
*   [ ] Test cross-compatibility between existing garments and avatars.
*   [ ] Generate simulation results for all valid garment-avatar combinations.
