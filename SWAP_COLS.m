% This is an Octave script for swapping columns of a matrix in a file
%
% Run this script using octave -qf path/to/SWAP_COLS.m <column a> <column b> name_of_input.file name_of_output.file
%

1;

% set up array of arguments supplied
arg_list= argv ();

% read in the input matrix
swapmat = dlmread(arg_list{3});

% swap column a for column b
swapmat(:,[arg_list{2} arg_list{1}]) = swapmat(:,[arg_list{1} arg_list{2}]);

% write the output
dlmwrite(arg_list{4},swapmat,'\t','precision','% 1.5f');
