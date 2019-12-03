%% Make domain
Xloc = [xmin xmax]; %lon
Yloc = [ymin ymax]; %lat



%% Initialise data matrix for ECCO2 data.
uvel = nan(12*(MaxYear-MinYear+1),50,xmax-xmin+1,ymax-ymin+1);
    
%Cycle through correct naming convention for cube92
year_ind = MinYear:MaxYear;
mon_ind = 1:12;
loop_ind=1;
clear monstr
for yy=MinYear:MaxYear; %number of years
 for mm=1:12 %number of months
  if mm<10
  monstr(loop_ind)=str2num([num2str(yy),'0',num2str(mm)]);
  elseif mm>=10
  monstr(loop_ind)=str2num([num2str(yy),num2str(mm)]);
  end
  loop_ind=loop_ind+1;
  end
end


for TimeInd = 1:12*(MaxYear-MinYear+1)

uvel(TimeInd,:,:,:) = permute(ncread([external_dir,'/ecco2/UVEL.nc/UVEL.1440x720x50.' num2str(monstr(TimeInd)) '.nc'],'UVEL',[Xloc(1) Yloc(1) 1 1],[Xloc(2)-Xloc(1)+1 Yloc(2)-Yloc(1)+1 Inf Inf]),[3 1 2]);

disp([num2str(TimeInd) ' of ' num2str(12*(MaxYear-MinYear+1)) ' month done.'])
end

save(fullfile(interim_dir,['cube92_iaf_uvel_',RunName,'.mat']),'uvel','-v7.3')
