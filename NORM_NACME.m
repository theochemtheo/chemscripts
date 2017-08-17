% This is an Octave script for normalizing non-adiabatic coupling vectors (matrix elements)
%
% Run this script using
% octave -qf ~/bin/NORM_NACME.m [name of NACME file]
% the names must be written WITHOUT suffixes, i.e. hp-modes-calc instead of hp-modes-calc.log
%
% To use this script, you must first:
%  + extract the non-adiabatic coupling vectors.
% This can be done using the script "Qchem44-EXTRACT_NACME.sh"
%
% The NACMEs must be in the following format:
% x-force [tab] y-force [tab] z-force
% and must be named [name].[MI]-to-[MJ].[type].mat, where [MI] and [MJ] are the multiplicity and index of states I and J respecitively, and [type] is the kind of NACME (e.g. gradient difference or derivative coupling with or without ETF)
%

1;

pkg load io

% set up array of arguments supplied
arg_list= argv ();

% declare variables for the arguments
basenameNACME = arg_list{1};

% set up search strings for the various files
searchstrDCnoETF = sprintf('%s*DCnoETF.mat',basenameNACME);
searchstrDCwithETF = sprintf('%s*DCwithETF.mat',basenameNACME);
searchstrNACV = sprintf('%s*NACV.mat',basenameNACME);

% construct sorted lists of the different types of files
listDCnoETF=sort_nat(glob(searchstrDCnoETF));
listNACV=sort_nat(glob(searchstrNACV));
listDCwithETF=sort_nat(glob(searchstrDCwithETF));

% Name of the DCnoETF output
DCnoETFout = sprintf('%s.DCnoETF.norms',basenameNACME);
% Name of the NACV output
NACVout = sprintf('%s.NACV.norms',basenameNACME);
% Name of the DCwithETF output
DCwithETFout = sprintf('%s.DCwithETF.norms',basenameNACME);
% Prepare a matrix of the correct dimensions with column for NACME identity and zeros elsewhere
magNACME = cell(length(listDCnoETF),2);

% start looping through the DCnoETF files
for i = 1:length(listDCnoETF);
	% Name of the file containing this vector
	nameNACME = listDCnoETF{i,1};
	% Trim to get the pair of states being looked at
	trimNACME = nameNACME(length(basenameNACME)+2:end-12);
	% Load in the vector
	NACME = dlmread(nameNACME);
	% Put the correct values into column 1
	magNACME(i,1) = trimNACME;
	% Put the correct values into column 2
	magNACME(i,2) = norm(NACME,"fro");
endfor
% write the output as a tab separated file with 6sf precision
cell2csv (DCnoETFout,magNACME)

% Initialize magNACME
magNACME = cell(length(listNACV),2);
% start looping through the DCnoETF files
for i = 1:length(listNACV);
	% Name of the file containing this vector
	nameNACME = listNACV{i,1};
	% Trim to get the pair of states being looked at
	trimNACME = nameNACME(length(basenameNACME)+2:end-7);
	% Load in the vector
	NACME = dlmread(nameNACME);
	% Put the correct values into column 1
	magNACME(i,1) = trimNACME;
	% Put the correct values into column 2
	magNACME(i,2) = norm(NACME,"fro");
endfor
% write the output as a tab separated file with 6sf precision
cell2csv (NACVout,magNACME)

% Initialize magNACME
magNACME = cell(length(listNACV),2);
% start looping through the DCnoETF files
for i = 1:length(listDCwithETF);
	% Name of the file containing this vector
	nameNACME = listDCwithETF{i,1};
	% Trim to get the pair of states being looked at
	trimNACME = nameNACME(length(basenameNACME)+2:end-14);
	% Load in the vector
	NACME = dlmread(nameNACME);
	% Put the correct values into column 1
	magNACME(i,1) = trimNACME;
	% Put the correct values into column 2
	magNACME(i,2) = norm(NACME,"fro");
endfor
% write the output as a tab separated file with 6sf precision
cell2csv (DCwithETFout,magNACME)

