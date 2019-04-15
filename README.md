# MAP Explorer
Use MAP to explore EHR data for individual patients. This Shiny app provides a interface to select patients and visualize their MAP probabilities.



# Running the App

 1. Get code from GitHub and install required packages.
 1. Copy required data files into a subdirectory called "data" within the source code directory. This directory will be ignored by git to prevent the data from being uploaded to GitHub.
 1. Launch Shiny app.

# Data

The app requires three data files as input:

- `MAPcutoff.Rdata`
- `MAPmanhattan.Rdata`
- `pheinfo.rda`