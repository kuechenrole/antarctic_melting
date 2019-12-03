function output_data = zero_out_land ( input_data, land )

nt = size(input_data,1);
for j = 1:nt
	slice = squeeze(input_data(j,:,:));
	slice(land) = 0;
	input_data(j,:,:) = slice;
end

output_data = input_data;


