# HaploLimit
Alleles and linkage information to estimate range of compatible haplotypes


This is the minimal working example with the initial example from Nick sent in 2015.
`./bin/testhaplo` simply gives bounds for all compatible haplotypes, and is being launched many times by the `haploscanner.py` master script.

Master script calculates haplo-boundaries for the original input data and then moves on the scanning procedure: fix haplotype frequency within the determined complatible range and recaluclate other boundaries accordingly. 




