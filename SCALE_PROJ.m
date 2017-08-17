% This is an Octave script for scaling the vibrational mode to NACME projection by wavenumber
%
% You must have the script "sort_nat.m" from the MathWorks File Exchange placed in /usr/share/octave/[version-number]/m/ for this to work.
%
% Run this script using
% octave -qf ~/bin/SCALE_PROJ.m [list of frequencies from hpmodes file.freqlist] [name of NACME file]
% the name of NACME file must be written WITHOUT suffixes, i.e. nacmes-calculation instead of nacmes-calculation.out, and the frequencies should be extracted with G09-EXTRACT_hpmodes_freqlist.sh
%
% To use this script, you must first:
%  + extract the frequencies of your normal modes.
%  + extract the non-adiabatic coupling vectors.
% This can be done using the scripts "G09-EXTRACT_hpmodes_freqlist.sh" and "Qchem44-EXTRACT_NACME.sh"
%
% The frequencies must be in the following format:
% number of mode [tab] frequency in wavenumbers
% The NACMEs must be in the following format:
% x-force [tab] y-force [tab] z-force
% and must be named [name].[MI]-to-[MJ].[type].mat, where [MI] and [MJ] are the multiplicity and index of states I and J respecitively, and [type] is the kind of NACME (e.g. gradient difference or derivative coupling with or without ETF)
%

1;

% load the misc package which is required for the normc command
pkg load miscellaneous

% set up array of arguments supplied
arg_list= argv ();

% declare variables for the arguments
freqlist = arg_list{1};
basenameNACME = arg_list{2};

% set up search strings for the various files
searchstrDCnoETF = sprintf('%s*DCnoETF.dist',basenameNACME);
searchstrDCwithETF = sprintf('%s*DCwithETF.dist',basenameNACME);
searchstrNACV = sprintf('%s*NACV.dist',basenameNACME);

% construct sorted lists of the different types of files
listDCnoETF=sort_nat(glob(searchstrDCnoETF));
listNACV=sort_nat(glob(searchstrNACV));
listDCwithETF=sort_nat(glob(searchstrDCwithETF));

FREQS = dlmread(freqlist);

% start looping through the DCnoETF files
for i = 1:length(listDCnoETF);
	% Name of the file containing this vector
	nameNACME = listDCnoETF{i,1};
	% Trim the name of .mat
	trimNACME = nameNACME(1:end-5);
	% Name of the two outputs
	output = sprintf('%s.scaled.dist',trimNACME);
	% Load in the vector
	distNACME = dlmread(nameNACME);
	% Prepare a matrix of the correct dimensions
	scaleddistNACME = zeros(length(FREQS),3);
	% Put the frequencies in column 1
	scaleddistNACME(:,1) = FREQS(:,2);
	% Convert the frequencies to atomic units
	auFREQ = (FREQS(:,2).*(4.55534E-6));
	% Scale the distances by frequency, e.g. projection/frequency
	scaleddistNACME(:,2) = (distNACME(:,2)./auFREQ);
	% Norm the second column so that the sum is 1 and put this into the 3rd column
	scaleddistNACME(:,3) = (scaleddistNACME(:,2)./(sum(scaleddistNACME(:,2))));
	% Output the normed scaled distance
	dlmwrite(output,scaleddistNACME,'\t','precision','% 1.7e');
endfor

% start looping through the NACV files
for i = 1:length(listNACV);
	% Name of the file containing this vector
	nameNACME = listNACV{i,1};
	% Trim the name of .mat
	trimNACME = nameNACME(1:end-5);
	% Name of the two outputs
	output = sprintf('%s.scaled.dist',trimNACME);
	% Load in the vector
	distNACME = dlmread(nameNACME);
	% Prepare a matrix of the correct dimensions
	scaleddistNACME = zeros(length(FREQS),3);
	% Put the frequencies in column 1
	scaleddistNACME(:,1) = FREQS(:,2);
	% Convert the frequencies to atomic units
	auFREQ = (FREQS(:,2).*(4.55534E-6));
	% Scale the distances by frequency, e.g. projection/frequency
	scaleddistNACME(:,2) = (distNACME(:,2)./auFREQ);
	% Norm the second column so that the sum is 1 and put this into the 3rd column
	scaleddistNACME(:,3) = (scaleddistNACME(:,2)./(sum(scaleddistNACME(:,2))));
	% Output the normed scaled distance
	dlmwrite(output,scaleddistNACME,'\t','precision','% 1.7e');
endfor

% start looping through the DCwithETF files
for i = 1:length(listDCwithETF);
	% Name of the file containing this vector
	nameNACME = listDCwithETF{i,1};
	% Trim the name of .mat
	trimNACME = nameNACME(1:end-5);
	% Name of the two outputs
	output = sprintf('%s.scaled.dist',trimNACME);
	% Load in the vector
	distNACME = dlmread(nameNACME);
	% Prepare a matrix of the correct dimensions
	scaleddistNACME = zeros(length(FREQS),3);
	% Put the frequencies in column 1
	scaleddistNACME(:,1) = FREQS(:,2);
	% Convert the frequencies to atomic units
	auFREQ = (FREQS(:,2).*(4.55534E-6));
	% Scale the distances by frequency, e.g. projection/frequency
	scaleddistNACME(:,2) = (distNACME(:,2)./auFREQ);
	% Norm the second column so that the sum is 1 and put this into the 3rd column
	scaleddistNACME(:,3) = (scaleddistNACME(:,2)./(sum(scaleddistNACME(:,2))));
	% Output the normed scaled distance
	dlmwrite(output,scaleddistNACME,'\t','precision','% 1.7e');
endfor

