%% Make domain
Xloc = [xmin xmax]; %lon
Yloc = [ymin ymax]; %lat



%% Initialise data matrix for ECCO2 data.
ssh = nan(12*(MaxYear-MinYear+1),xmax-xmin+1,ymax-ymin+1);

ssh(:,:,:) = permute(ncread([external_dir,'/ecco2/SSH.nc/SSH.1440x720.1996_2016_monmean.nc'],'SSH',[Xloc(1) Yloc(1) 1],[Xloc(2)-Xloc(1)+1 Yloc(2)-Yloc(1)+1 Inf]),[3 1 2]);
save(fullfile(interim_dir,['cube92_iaf_ssh_',RunName,'.mat']),'ssh','-v7.3')
