# Stimuli-Responsive_Emulsions
The following Command is used to run the code:

python3 Code.py &arg1 &arg2 &arg3 &arg4 &arg5 &arg6 &arg6 &arg7 &arg8 &arg9 &arg10 &arg11 nohup &
&arg1: Random Seed Value
&arg2: Interface energy, ε_ext (larger than the bead-bead interaction scale)
&arg3: Effective interaction range of (one) interface potential, σ_ext (smaller than the bead size) 
&arg4: Distance of interface from z=0 axis. 
&arg5: Length of the LJ tail 
&arg6: Timesteps
&arg7: size of the simulation box in x,y direction
&arg8: Size of the simulation box in z direction
&arg9: Value of alpha parameter
&arg10: kT
&arg11: Value of R_c defined in the paper 

Sample command to run the code: 

python3 Code.py 1 5.5 0.5 18 100.0 10000000 150 150 0.0 1.0 40 nohup &

Sample Input file: "Initial_Configuration.gsd"
Sample Output Files:"Final_Configuration.gsd"
		   "Trajectory.gsd"
                    "Forces.txt": This is used to calculate the mean force between all the beads.
                     awk '{printf ($3-$6)/2; printf "\n"}' Forces.txt>>Data.txt
