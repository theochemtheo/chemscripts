% This is an Octave script for checking if there has been a state-crossing between two steps of a potential energy scan and setting the energy to be relative
%
% Run this script using
% octave -qf CROSSOVER-CHECK+EREL.m [step_N-1.diabmat] [step_N.diabmat] [step_M.diabmat] [name of output] [Erel of SCF in a.u.]
%
% To use this script, you must first:
%  + extract the relevant diabatic matrices
% This can be done using the script "Qchem44-EXTRACT_diabats.sh"
%

1;

% Read in matrices O, N and M
% First, check that the file for steps N and N-1 exist. If not, replace them with M if neither exist. If only N exists, use N for N-1
if exist(argv(){2},"file") == 2;
	NMAT = dlmread(argv(){2});
	if exist(argv(){1},"file") == 2;
		OMAT = dlmread(argv(){1});
	else
		OMAT = dlmread(argv(){2});
	endif
else
	NMAT = dlmread(argv(){3});
	OMAT = dlmread(argv(){3});
endif
MMAT = dlmread(argv(){3});
SCFEREL = str2num(argv(){5});

% Dimensionality of the matrices
DIM = columns(NMAT);
% Prepare the [c, x] matrix for steps O and N for later linear fit
x = [1,-1;1,0];
% Prepare the vector for Q for later linear fit
Q = [1,1];
% Add the relative energy of the SCF to the diagonal values of MMAT
MMAT = MMAT+eye(DIM)*SCFEREL;

% Loop through the states k of M
for k = 1:DIM;
	% Get the energy of state k of M -> Mk
	Ek = MMAT(k,k);
	% Form the vector of the off-diagonal elements related to state Mk
	Mk = vertcat(MMAT(k,[1:(k-1)])',MMAT([(k+1):DIM],k));
	% Loop through the states l of N -> Nl
	for l = 1:DIM;
		% Form the vector of the off-diagonal elements related to state Nl
		Nl = vertcat(NMAT(l,[1:(l-1)])',NMAT([(l+1):DIM],l));
		% Prepare the vector of the El at steps O and N
		V = [OMAT(l,l);NMAT(l,l)];
		% Fit the values of m and c for El = mx + c
		PLY = ((pinv(x'*x))*x'*V);
		% Expectated value of El at M given the polynomial fitted through O and N
		PEl = Q*PLY;
		% Project Mk onto Nl apart from their shared element and divide by the square of the difference  energies -> Store in a vector
		Skl(l,1) = ((dot(Mk,Nl)/(norm(Mk)*norm(Nl)))**2)/((Ek-PEl)**4);
	endfor
	% Find the maximal value of the similarity vector -> store into a vector containing that value and its index
	[MAX,iMAX] = max(abs(Skl));
	% If iMAX is not equal to k
	if iMAX != k
		% Swap rows iMAX and k
		MMAT([iMAX k],:) = MMAT([k iMAX],:);
		% Swap columns iMAX and k
		MMAT(:,[iMAX k]) = MMAT(:,[k iMAX]);
	endif
endfor

% Write the permuted matrix
dlmwrite(argv(){4},MMAT,'\t','precision','% 1.10f');
