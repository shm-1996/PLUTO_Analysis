#!/bin/bash
#PBS -P ek9
#PBS -q express
#PBS -l walltime=02:00:00
#PBS -l ncpus=4
#PBS -l mem=8GB
#PBS -l wd
#PBS -N Compute_b
#PBS -M shyam.menon@anu.edu.au
#PBS -m bea
#PBS -l storage=scratch/ek9+gdata/ek9

#eta_0.5
python b_TimeEvolution.py -directory /g/data1b/ek9/sm5890/Driving_Parameter/Turbulence_Initial/eta_0.5/Mach_10/ -N 200 -tstart 0 -tend 140 -direct
python b_TimeEvolution.py -directory /g/data1b/ek9/sm5890/Driving_Parameter/Turbulence_Initial/eta_0.5/Mach_5/ -N 200 -tstart 0 -tend 300 -direct
python b_TimeEvolution.py -directory /g/data1b/ek9/sm5890/Driving_Parameter/Turbulence_Initial/eta_0.5/Mach_20/ -N 200 -tstart 0 -tend 70 -direct

#eta_0.75
python b_TimeEvolution.py -directory /g/data1b/ek9/sm5890/Driving_Parameter/Turbulence_Initial/eta_0.75/Mach_10/ -N 200 -tstart 0 -tend 140 -direct
python b_TimeEvolution.py -directory /g/data1b/ek9/sm5890/Driving_Parameter/Turbulence_Initial/eta_0.75/Mach_5/ -N 200 -tstart 0 -tend 300 -direct
python b_TimeEvolution.py -directory /g/data1b/ek9/sm5890/Driving_Parameter/Turbulence_Initial/eta_0.75/Mach_20/ -N 200 -tstart 0 -tend 70 -direct

#eta_0.0
python b_TimeEvolution.py -directory /g/data1b/ek9/sm5890/Driving_Parameter/Turbulence_Initial/eta_0.0/Mach_10/ -N 200 -tstart 0 -tend 140 -direct

python b_TimeEvolution.py -directory /scratch/ek9/sm5890/Pillar_Simulations/Recombination_Field/ -N 200 -tstart 140 -tend 300 -table

python b_TimeEvolution.py -directory /g/data1b/ek9/sm5890/Driving_Parameter/self_gravity/ -N 200 -tstart 140 -tend 340 -table

python ColumnDensityPlots.py -directory /g/data1b/ek9/sm5890/Driving_Parameter/self_gravity/ -N 200 -tstart 140 -tend 340 