clear all;
close all;

graphdevice = '-dpng'; % select graphic device
graphformat = '.png';  % select graphics format

% Parse the high resolution met data for one day

%%% do this manually:
% disp('-----------------------------------------------------');
% disp('do this manually first:');
% disp('cp *.MET ~/DATA/PISTON_2019/Sally_Ride/ship_day/data/');
% disp('-----------------------------------------------------');

savedir = '/Users/eliz/DATA/ASTRAL_2023/Revelle/ship/';
datadir = savedir;
files = dir([datadir '*.MET']);
nfiles = length(files);

halfway = round(length(files)/2);
lastone = length(files)-1;

for i = 1:nfiles           %%% 9/2 through 9/9    ... already done on 9/15
% for i = 9:15            %%% 9/10 (pt 1 and pt 2) through 9/15  ...
% for i = nfiles     %%% 9/16 through 9/25  ... numbered sequentially conveninently

    clear bigdata ship_day fname a;

    fname = [datadir files(i).name];
    disp(['** reading in file: ' fname]);
    try
    eval(['!grep -v "#" ' fname ' > myfile'])
    a = load('myfile');
    
    %%% because file headers changed mid cruise, don't save any data past
    %%% column 76. Otherwise program will crash because column numbers will
    %%% be different. 
    bigdata = a(:,1:101);
    catch
        disp('cant do it')
    end

delete myfile

%%% comment out stuff you don't want

ship_day.ta         = bigdata(:,2); % C
ship_day.pa         = bigdata(:,3); % mb
ship_day.psealevel  = bigdata(:,4); % mb
ship_day.paccum     = bigdata(:,5); % mm
ship_day.rh         = bigdata(:,6); % %
ship_day.td         = bigdata(:,8); % dewpoint, C
ship_day.rh2        = bigdata(:,9); % 2nd RH sensor? %
ship_day.lw_dome_t  = bigdata(:,12); % lw dome temp, K
ship_day.lw_case_t  = bigdata(:,13); % lw case temp, K
ship_day.lw_therm   = bigdata(:,14); % thermopile volts
ship_day.lw_dn      = bigdata(:,15); % longwave radiation, W/m2
ship_day.sw_dn      = bigdata(:,16); % solar radiation, W/m2
ship_day.rspd       = bigdata(:,18); % m/s
ship_day.rdir       = bigdata(:,19); % relative, from, deg
ship_day.wspd       = bigdata(:,20); % m/s
ship_day.wdir       = bigdata(:,21); % true, from, deg
ship_day.lat        = bigdata(:,67); % deg
ship_day.lon        = bigdata(:,68); % deg
ship_day.time       = bigdata(:,1); % ???
ship_day.cog        = bigdata(:,70); % deg
ship_day.sog        = bigdata(:,71); % kt
ship_day.gps_time   = bigdata(:,72); % GMT sec since 1/1/1970
ship_day.gps_timeofday  = bigdata(:,69); % GMT sec 0-86400
ship_day.tsea1      = bigdata(:,26); % C 
ship_day.tsea2      = bigdata(:,34); % C
% ship_day.tsea3      = bigdata(:,94); % C
% ship_day.tsea4      = bigdata(:,97); % C
ship_day.csea1      = bigdata(:,27); % mS/cm 
ship_day.csea2      = bigdata(:,35); % mS/cm 
% ship_day.csea3      = bigdata(:,95); % mS/cm 
% ship_day.csea4      = bigdata(:,95); % mS/cm 
ship_day.ssea1      = bigdata(:,28); % psu 
ship_day.ssea2      = bigdata(:,36); % psu 
% ship_day.ssea3      = bigdata(:,96); % psu 
% ship_day.ssea4      = bigdata(:,99); % psu 
ship_day.sigt1      = bigdata(:,29); % kg/m^3
ship_day.sigt2      = bigdata(:,37); % kg/m^3 
% ship_day.sigt3      = bigdata(:,97); % kg/m^3 
% ship_day.sigt4      = bigdata(:,97); % kg/m^3 
ship_day.flowsea1   = bigdata(:,31); % L/min 
ship_day.flowsea2   = bigdata(:,39); % L/min 
% ship_day.flowsea3   = bigdata(:,98); % L/min 
% ship_day.flowsea4   = bigdata(:,98); % L/min 
ship_day.pitch      = bigdata(:,47); % VRU deg
ship_day.roll       = bigdata(:,48); % VRU deg
ship_day.heave      = bigdata(:,49); % VRU m
ship_day.hed_ash    = bigdata(:,55); % Ashtech deg
ship_day.pitch_ash  = bigdata(:,56); % Ashtech deg
ship_day.roll_ash   = bigdata(:,57); % Ashtech deg
ship_day.hed_gyr    = bigdata(:,58); % Ashtech deg
ship_day.spdlog_v   = bigdata(:,100); % speed log longitudinal deg
ship_day.spdlog_u   = bigdata(:,101); % speed log transverse deg


%date in matlab datenum and jd year day format
ship_day.t = ship_day.gps_time/86400 + datenum(1970,1,1,0,0,0);
ship_day.jd = ship_day.gps_time/86400 + datenum(1970,1,1,0,0,0)- datenum(2023,0,0,0,0,0);


%% plots
datest = datestr(min(ship_day.t),'mmDD');
plot_checks = 0;
if plot_checks == 1
pathplot = '/Users/eliz/DATA/ASTRAL_2023/Revelle/flux/Raw_Images/ship_day/';

figure;
counter = 1;
thefields = fields(orderfields(ship_day));
nc = length(thefields);
for i = 1:4:nc
    clf;
    for j = 1:4
       if (i+j-1) <= nc
           subplot(2,2,j); hold on;
           eval(['thevar = ship_day.' thefields{i+j-1} ';']);
           plot(ship_day.t, thevar,'o');
           grid on;
           var_name = {strrep(thefields{i+j-1},'_',' ')};
           title([datest ' ' var_name]);
           xlim([min(ship_day.t) max(ship_day.t)]);
           datetick('x','keeplimits');
           grid on;
       end
    end
   print(graphdevice,[pathplot 'ship_ASTRAL_2023_' datest '_' sprintf('%i',counter) graphformat]);
   counter = counter + 1;
end  % for all vars
end % if plotting

% part_st2 = datestr(max(ship_day.t),'mmDD');
save([savedir 'MET_Revelle_ASTRAL_2023_' datest], 'ship_day');

end  % loop through files