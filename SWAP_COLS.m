% This is an Octave script for swapping columns of a matrix in a file
%
% Run this script using octave -qf path/to/SWAP_COLS.m <column a> <column b> name_of_input.file name_of_output.file
%

1;

% get the arguments out as the right kind of type
cola = str2num(argv(){1});
colb = str2num(argv(){2});
inmatname = argv(){3};
outmatname = argv(){4};

% read in the input matrix
swapmat = dlmread(inmatname);

% swap column a for column b
swapmat(:,[colb cola]) = swapmat(:,[cola colb]);

% write the output
dlmwrite(outmatname,swapmat,'\t','precision','% 1.5f');
