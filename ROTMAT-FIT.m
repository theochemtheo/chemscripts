% This is an Octave script for fitting rotation matrices
%
% Run this script using
% octave -qf ROTMATFIT.m [value of Q] [name of output]
% in the directory that contains all of the rotation matrices you are interested in
%
% To use this script, you must first:
%  + extract the relevant rotation matrices
% This can be done using the script "Qchem44-EXTRACT_NACME.sh"
% 

1;

pkg load io

% The value for Q to be used for the output
Q = str2num(argv(){1});
% Name of the output
OUTNAME = argv(){2};

% Dimensionality of the rotation matrices?
DIM = 6;

% Names of the files to be fitted
L02 = dlmread('CMIX_T-Boys-10C6_349_02l.rotmat');
L01 = dlmread('CMIX_T-Boys-10C6_349_01l.rotmat');
Q00 = dlmread('CMIX_T-Boys-10C6_349_00q.rotmat');
R01 = dlmread('CMIX_T-Boys-10C6_349_01r.rotmat');
R02 = dlmread('CMIX_T-Boys-10C6_349_02r.rotmat');
% The corresponding matrix for X
x = [1,-.02,.0004,-.000008,.00000016;1,-.01,.0001,-.000001,.00000001;1,0,0,0,0;1,.01,.0001,.000001,.00000001;1,.02,.0004,.000008,.00000016];

%%%
% Masking
%%%
% Get the sign mask from Q00
MASK = lt(Q00,0)*-2+ones(6);
% Mask the other rotmats
ML02 = MASK.*abs(L02);
ML01 = MASK.*abs(L01);
MR01 = MASK.*abs(R01);
MR02 = MASK.*abs(R02);


% Initialize a matrix for the output
OUTPUT = cell(DIM,DIM);
% Loop through the entries
for i = 1:DIM*DIM;
	% Prepare the vector containing correct matrix elements from the rotation matrices and the different Qs
	V = [ML02(i);ML01(i);Q00(i);MR01(i);MR02(i)];
	% Fit the polynomial for this matrix element
	PLY = ((pinv(x'*x))*x'*V);
	% Prepare the vector for this Q
	K = [1,Q,Q**2,Q**3,Q**4];
	% Put the matrix element into the output
	OUTPUT{i} = sprintf('% 1.11f',K*PLY);
endfor

% write the output for this Q, note that there is a literal tab (ctrl+I) in quotes
cell2csv (OUTNAME,OUTPUT,'	');

