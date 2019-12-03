%% Make domain
Xloc = [xmin xmax]; %lon
Yloc = [ymin ymax]; %lat



%% Initialise data matrix for ECCO2 data.
salt = nan(12*(MaxYear-MinYear+1),50,xmax-xmin+1,ymax-ymin+1);
    
%Cycle through correct naming convention for cube92 data
    %for i = [1+12(year_ind-MinYear)]:[12(year_ind-MinYear+1)]

    salt(:,:,:,:) = permute(ncread([external_dir,'/ecco2/SALT.nc/SALT.2013_2014_monmean.nc'],'SALT',[Xloc(1) Yloc(1) 1 1],[Xloc(2)-Xloc(1)+1 Yloc(2)-Yloc(1)+1 Inf Inf]),[4 3 1 2]);
    
  %  disp([num2str(TimeInd) ' of ' num2str(12*(MaxYear-MinYear+1)) ' month done.'])


save(fullfile(interim_dir,['cube92_iaf_salt_xmd_test_',RunName,'.mat']),'salt','-v7.3')


