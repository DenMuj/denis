# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/bin:$HOME/bin2vtk

export PATH
module load intel/2019
module load openmpi/4.0.5-intel
module load fftw3/intel/3.3.8-openmpi4
module load sprng/5.0
