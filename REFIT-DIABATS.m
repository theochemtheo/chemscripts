% This is an Octave script for fitting polynomials to sets of diabatic state matrices

1;

% Load the optim package to use nonline curve fitting
pkg load optim

% Order to which the fitting will occur
ORD = str2num(argv(){1});
% Make an array containing the list of the files to be fitted
LIST = importdata('fitting.lst');

if (ORD >= length(LIST));
	printf ("The order of the desired polynomial must be lower than the number of points being fit!");
	exit
endif

% Initialize the matrix X
X = zeros(length(LIST),(ORD+1));
X(:,1) = 1;
% Fill in the matrix for x by looping through the elements of the list and set up the tensor containing the matrices at each point
for k = 1:length(LIST);
	% Get the magnitude of Q from the name of the file
	Qk = (str2num((regexp(LIST{k},'.{2}(?=[l,r,q])','match tokens')){:})/100);
	% If k is on the left, Q is negative
	if (strncmp((regexp(LIST{k},'[l,r,q]','match tokens')){:},"l",1) == 1)
	Qk = -Qk;
	endif
	for l = 1:ORD;
		X(k,(l+1)) = Qk**l;
	endfor
	V(:,:,k) = dlmread(LIST{k});
endfor

% Equation for the Morse potential
% q is the position along Q, p(1)=De is the depth of the well, p(2)=a is related to the force constant by a = sqrt(k/2De), where k is the harmonic force constant and p(3)=qe is the eqm value of Q.
MORSE= @ (p, q) p(1)*((1-e.^(-p(2)*(q-p(3)))).^2);
% Initial conditions for fitting a Morse oscillator
%      note that we use the 2nd derivative of the harmonic fit to the values to get our initial guess for a
init = [3;sqrt(polyfit(X(:,2),vec(V(1,1,:)),2)(3));0];
% Define the independent variable and the observable
indep = X(:,2);
obs =  vec(V(1,1,:));

vec(V(1,1,:))
X(:,2)
polyfit(X(:,2),vec(V(1,1,:)),2)
%(pinv(X'*X))*X'*vec(V(1,1,:))

%% Fit the Morse
%[p,model_values] = nonlin_curvefit(MORSE, init, indep, obs);
%p
