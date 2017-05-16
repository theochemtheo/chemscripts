% This is an Octave script for reformatting diabatic matrices
%
% Run this script using
% octave -qf REFORM-DIABMAT.m [name of mask matrix] [name of input] [basename of output]
% in the directory that contains all of the diabatic matrices you are interested in
%
% To use this script, you must first:
%  + extract the relevant diabatic matrices
% This can be done using the script "Qchem44-EXTRACT_diabats.sh"
% 

1;

pkg load io

% Name of the mask matrix
MASKNAME = argv(){1};
% Name of the input matrix
INNAME = argv(){2};
% Name of the output
OUTNAME = argv(){3};

% Central point to take mask from
MASKMAT = dlmread(MASKNAME);
% Input matrix to be masked
INMAT = dlmread(INNAME);

% Get the sign mask from 
MASK =  MASKMAT./(sqrt(MASKMAT.*MASKMAT));
% Mask the point
OUTMAT = MASK.*(sqrt(INMAT.*INMAT));
% Get the eigenvalues of the new diabatic matrix
OUTEIG = eig(OUTMAT);
% Write the output diabatic matrix
dlmwrite(sprintf('%s.diabmat',OUTNAME),OUTMAT,'\t','precision','% 1.10f');
% Write the output eigenvalues
dlmwrite(sprintf('%s.adiabex',OUTNAME),OUTEIG,'\t','precision','% 1.12f');

