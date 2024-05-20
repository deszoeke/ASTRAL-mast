function ship_day = read_ship_day(path_working_ddd,ddd,yyyy,PosLims,SCS_adj,zp_ship,zq_ship)
%{
Reads daily ship files at 1-sec interval for one day
returns structure with concatenated columns of each field
should be length of time = 86400 = 24 x 3600. 


Nans will be present when data are not. 

input parameters: path_working_ddd = path to daily data folder
                  ddd = julian date
                  yyyy = year string
                  PosLims = lat/lon limits for filtering position data
%}

fclose all;
  
    thedate = doy2date(ddd,yyyy);
    datestamp = datestr(thedate,'mmDD');
%     MET_Revelle_ASTRAL_2023_0609
    dfl = fullfile(path_working_ddd,['MET_Revelle_ASTRAL_2023_',datestamp,'.mat']);
    [ship_day] = read_ship_astral(dfl,ddd,yyyy,PosLims,SCS_adj,zp_ship,zq_ship);

end  %% end function
