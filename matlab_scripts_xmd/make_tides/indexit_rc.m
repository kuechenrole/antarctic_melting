function [r,c] = indexit_rc ( data_size, matlab_index )
% INDEXIT_RC:  converts column major index number to row and column
%
% USAGE:  [r,c] = indexit_rc ( data_size, matlab_index )      

switch ( length(data_size) )
case 2
	num_rows = data_size(1);
	num_cols = data_size(2);
end

r = mod(matlab_index-1,num_rows)+1;
c = (matlab_index-r)/num_rows + 1;

