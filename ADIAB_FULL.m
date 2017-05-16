% This is an Octave script for normalizing non-adiabatic coupling vectors (matrix elements)
%
% Run this script using
% octave -qf ~/bin/ADIAB_FULL.m [name of the base file]
% the names must be written WITHOUT suffixes, i.e. hp-modes-calc instead of hp-modes-calc.log
%
% To use this script, you must first:
%  + extract the adiabatic excitation energies
% This can be done using the script "Qchem44-EXTRACT_adiabat.sh"
% 

1;

% set up array of arguments supplied
arg_list= argv ();

% declare variables for the arguments
basenameNACME = arg_list{1};

% set up search strings for the various files
adiabEX = sprintf('%s.adiabex',basenameNACME);

% Name of the DCnoETF output
fullEX = sprintf('%s.adiabfull',basenameNACME);

% Form the vector of the excitations
a = dlmread(adiabEX);

% Make the diagonal excitation matrix from the vector
A = diag(a);

% Write the output
dlmwrite(fullEX, A, "precision", "% 1.12f", "delimiter", "\t");
