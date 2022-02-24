# Verlet-integration Speedtest

This repository contains a set of small programs to assess the speed-up of particle trajectory calculations (using Verlet-integration) in a background air flow current. I made this because I wanted to try out the speed-up with OpenMP and compare the use of Eigen library to represent the position and velocity arrays. Eigen library will make the program much more elegant to read, but the question is at what cost?

- rotVelField.m: Illustration of the flow field as a Matlab quiver plot.
- verletTestArray.cpp: Main program using array representation of the position and velocity vectors.