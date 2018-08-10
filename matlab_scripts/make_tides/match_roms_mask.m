function output_data = match_roms_mask ( lon, lat, mask, data )
% MATCH ROMS MASK:  make the otps data match the roms land mask.
%
% The data is filled in used the mode of the nearby data rather
% than the mean.  This way the phase doesn't get screwed up as bad.
%
% USAGE:  output_data = match_roms_mask ( lon, lat, mask, data );
%

filled_in_data = data;

%
% Find the ROMS water points (mask == 1) which the OTPS grid thought 
% was water (input_constituent == nan).
nan_inds = find ( isnan(data) & mask==1);


%
% Get the row and column numbers of these points.
[nan_r, nan_c] = indexit_rc ( size(lon), nan_inds );
[num_rows,num_cols] = size(lon);


%
% loop thru each point, trying to fill it in with the nearest non-nan points
num_points = length(nan_inds);
for j = 1:num_points

	%
	% Construct an index mask around this point
	nan_index = nan_inds(j);
	mask_thickness = 1;
	while 1

		nbrows = [nan_r(j)-mask_thickness:nan_r(j)+mask_thickness];
		ind = find((nbrows<1) | (nbrows>num_rows));
		nbrows(ind) = [];

		nbcols = [nan_c(j)-mask_thickness:nan_c(j)+mask_thickness];
		ind = find((nbcols<1) | (nbcols>num_cols));
		nbcols(ind) = [];

		[C,R] = meshgrid(nbcols,nbrows);
		neighbor_inds = (C-1)*num_rows + R;

		%
		% check that none of these neighbor points are land points
		valid_ind = find(mask(neighbor_inds) == 1);
		neighbor_inds = neighbor_inds(valid_ind);

		finite_ind = find(isfinite(data(neighbor_inds)));
		if any(finite_ind)
			break;
		end

		mask_thickness = mask_thickness + 1;
		if ( mask_thickness > 35 )
		%	error ( 'how can it be this bad?' );
            break;
		end
    end

    if (mask_thickness <=35 )
	    filled_in_data(nan_index) = my_mode(data(neighbor_inds(finite_ind)));
	%filled_in_data(nan_index) = mean(data(neighbor_inds(finite_ind)));
    else
        filled_in_data(nan_index) = 0.0;
    end

end

output_data = filled_in_data;

return

function the_mode = my_mode ( input_data )
mode_data = sort ( input_data );
num_mode_points = length(input_data);
if floor(num_mode_points/2) ~= num_mode_points/2
	mode_index = ceil(num_mode_points/2);
	the_mode = mode_data(mode_index);
else
	mode_index = [num_mode_points/2 num_mode_points/2+1];
	the_mode = mean ( mode_data(mode_index) );
end

