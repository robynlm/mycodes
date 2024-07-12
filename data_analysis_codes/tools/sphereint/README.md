# sphereint
Class to compute the numerical volume in a sphere and area of a sphere.
    
This class provides a weight for each grid position based on whether or not it is in (weight = 1), out (weight = 0), or partially in (weight in between 0 and 1) a sphere of a given radius.

A cubic cell is placed around each grid position and the volume of the cell in the sphere (assuming a flat suface in the cell) is calculated and normalised by the cell volume to obtain the weight.

This has been extended to calculate the area of the sphere with the option domain='area'.

If you use this code please reference
@article{R.L.Munoz_M.Bruni_2023,
    title     = {Structure formation and quasispherical collapse from initial curvature perturbations with numerical relativity simulations},
    author    = {Munoz, R. L. and Bruni, M.},
    journal   = {Physical Review D},
    volume    = {107},
    number    = {12},
    pages     = {123536},
    numpages  = {26},
    year      = {2023},
    month     = {6},
    doi       = {10.1103/PhysRevD.107.123536},
    archivePrefix = {arXiv},
    eprint    = {astro-ph/2302.09033}}