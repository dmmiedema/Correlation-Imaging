# Correlation-Imaging
Matlab script to calculate density, velocity and run length of molecular motors on microtubules using correlation imaging (CI).

This program calculates the velocity and run length of particles moving along a straight segment, as well as the (1D) density of particles on the segment. The method is based on spatiotemporal correlations of the intensities from the moving objects. It has been developed to measure the density, velocity and run length of fluorescently labeled molecular motors moving along microtubule segments, and is applicable also at high motor density. 

The main function is divided in 4 analysis parts:
1) Segment selection from an image sequence and obtaining the intensity
   along the segment over time. (uses sub-function 'Im2In').
2) Determine the particle density along the segment from the spatially
   average intensities using fluorescence correlation spectroscopy.
3) Calculate the spatiotemporal correlation of intensities. (uses
   sub-function 'correlation')
4) Derive particle velocity and run length from analysis of the
   evolution of the correlation peak. (uses sub-function 'motility').
