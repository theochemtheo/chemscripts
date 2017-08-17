% This is an Octave script for projecting vibrational modes onto non-adiabatic coupling vectors (matrix elements)
%
% Run this script using
% octave -qf ~/bin/PROJ_FREQ-to-NACME.m [name of frequencies file] [name of NACME file]
% the names must be written WITHOUT suffixes, i.e. hp-modes-calc instead of hp-modes-calc.log
%
% To use this script, you must first:
%  + extract the high-precision, (Frobenius) normalized, mass-weighted cartesian displacement matrices for the normal modes.
%  + extract the non-adiabatic coupling vectors.
% This can be done using the scripts "G09-EXTRACT_hpmodes.sh" and "Qchem44-EXTRACT_NACME.sh"
%
% The normal modes must be in the following format:
% x-displacement [tab] y-displacement [tab] z-displacement
% and must be named [name].m[number].dat, where [number] is the number for the mode. The matrices should be Frobenius normalized.
% The NACMEs must be in the following format:
% x-force [tab] y-force [tab] z-force
% and must be named [name].[MI]-to-[MJ].[type].mat, where [MI] and [MJ] are the multiplicity and index of states I and J respecitively, and [type] is the kind of NACME (e.g. gradient difference or derivative coupling with or without ETF)
%

1;

% set up array of arguments supplied
arg_list= argv ();

% declare variables for the arguments
basenameFREQ = arg_list{1};
basenameNACME = arg_list{2};

% set up search strings for the various files
searchstrFREQ = sprintf('%s*dat',basenameFREQ);
searchstrDCnoETF = sprintf('%s*DCnoETF.mat',basenameNACME);
searchstrDCwithETF = sprintf('%s*DCwithETF.mat',basenameNACME);
searchstrNACV = sprintf('%s*NACV.mat',basenameNACME);

% construct sorted lists of the different types of files
listFREQ=sort_nat(glob(searchstrFREQ));
listDCnoETF=sort_nat(glob(searchstrDCnoETF));
listNACV=sort_nat(glob(searchstrNACV));
listDCwithETF=sort_nat(glob(searchstrDCwithETF));

% start looping through the DCnoETF files
for i = 1:length(listDCnoETF);
	% Name of the file containing this vector
	nameNACME = listDCnoETF{i,1};
	% Trim the name of .mat
	trimNACME = nameNACME(1:end-4);
	% Name of the output
	output = sprintf('%s.dist',trimNACME);
	% Load in the vector
	NACME = dlmread(nameNACME);
	% Frobenius norm the vector
	normNACME = NACME./norm(NACME,"fro");
	% Prepare a matrix of the correct dimensions with column for mode number and zeros elsewhere
	distNACME = [colon(1,length(listFREQ))',zeros(1,length(listFREQ))'];
	% loop through the frequencies
	for j = 1:length(listFREQ);
		% put the similarity of NACME to frequency j into the jth column
		distNACME(j,2) = sqrt(abs(trace(normNACME'*dlmread(listFREQ{j,1}))));
	endfor
	% write the output as a tab separated file with 6sf precision
	dlmwrite(output,distNACME,'\t','precision','% 1.6e');
endfor

% start looping through the DCnoETF files
for i = 1:length(listNACV);
	% Name of the file containing this vector
	nameNACME = listNACV{i,1};
	% Trim the name of .mat
	trimNACME = nameNACME(1:end-4);
	% Name of the output
	output = sprintf('%s.dist',trimNACME);
	% Load in the vector
	NACME = dlmread(nameNACME);
	% Frobenius norm the vector
	normNACME = NACME./norm(NACME,"fro");
	% Prepare a matrix of the correct dimensions with column for mode number and zeros elsewhere
	distNACME = [colon(1,length(listFREQ))',zeros(1,length(listFREQ))'];
	% loop through the frequencies
	for j = 1:length(listFREQ);
		% put the similarity of NACME to frequency j into the jth column
		distNACME(j,2) = sqrt(abs(trace(normNACME'*dlmread(listFREQ{j,1}))));
	endfor
	% write the output as a tab separated file with 6sf precision
	dlmwrite(output,distNACME,'\t','precision','% 1.6e');
endfor

% start looping through the DCnoETF files
for i = 1:length(listDCwithETF);
	% Name of the file containing this vector
	nameNACME = listDCwithETF{i,1};
	% Trim the name of .mat
	trimNACME = nameNACME(1:end-4);
	% Name of the output
	output = sprintf('%s.dist',trimNACME);
	% Load in the vector
	NACME = dlmread(nameNACME);
	% Frobenius norm the vector
	normNACME = NACME./norm(NACME,"fro");
	% Prepare a matrix of the correct dimensions with column for mode number and zeros elsewhere
	distNACME = [colon(1,length(listFREQ))',zeros(1,length(listFREQ))'];
	% loop through the frequencies
	for j = 1:length(listFREQ);
		% put the similarity of NACME to frequency j into the jth column
		distNACME(j,2) = sqrt(abs(trace(normNACME'*dlmread(listFREQ{j,1}))));
	endfor
	% write the output as a tab separated file with 6sf precision
	dlmwrite(output,distNACME,'\t','precision','% 1.6e');
endfor

