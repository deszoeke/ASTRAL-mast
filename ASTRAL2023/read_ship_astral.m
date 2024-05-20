function [ship_day] = read_ship_astral(dfl,ddd,yyyy,PosLims,SCS_adj,zpim, zq_ship)

%{
reads ship_day file specified by dfl and returns array with the
following columns:

test with: ship_day = read_ship(dfl,ddd,hhh,PosLims,SCS_adj,zpim)
% path_working_ddd = '/Users/eliz/DATA/PISTON_2019/Sally/flux/Raw/19248';
% dfl = '/Users/eliz/DATA/PISTON_2019/Sally/flux/Raw/19248/scs01924804_raw.txt';
% ddd = 9;  % yearday
% hhh = 5;  % hour
% PosLims = [119 139 4 27];   % limits for map
% SCS_adj = [0 0];
% zpim = 157;

% testing of:
    %%% input parameters: dfl = string path to ship_day file
                  ddd = julian date
                  hhh = hour
                  PosLims = lat/lon limits for filtering position data

    %%% output: currently 38 parameters but could change each cruise.


    %%% input/output data 
%     t             (calculated here) matlab date/time
%     jd            (calculated here) julian date, need to print out past 5th decimal place to see the precision to 1 sec    
%     lat           Lat DDDMM.MM, primary gps
%     lon           Lon DDDMM.MM, primary gps
%     cog           COG deg, primary gps
%     sog           SOG kts, primary gps * 0.514 for m/s
%     Ttsg          T SBE38, C
%     Stsg          S SBE45, psu
%     rh            RH, %
%     Ta            air T, C
%     rs            solar, W/m2
%     rl            IR, W/m2
%     hed           heading POSMV, deg
%     wspd          true wind speed, kt, foremast 2D sonic 
%     wdir          true wind dir, deg, foremast 2D sonic anemometer
%     psealevel_reported  pressured corrected for MSL by ship_day, mb
%     
%     pitch         pitch POSMV, deg
%     roll          roll POSMV, deg
%     heave         heave POSMV, m
%     rspd          relative wind speed, kt, foremast 2D sonic
%     rdir          relative wind direction, deg, foremast 2D sonic
%  
%     Td            dew point T, C
%     pa            measured barometric pressure, mb
%     paccum        accumulated precipitation, mm

%%% starboard bow TSG
%     Ftsg1          fluorometer dry value micrograms/L
%     Ctsg1          C SBE45 (S/m)
%     Stsg1          S SBE45 (psu)
%     Ttsg1_int      T of internal sensor SBE45 (C)

%%% aft porside TSG
%     Ftsg2          fluorometer dry value micrograms/L
%     Ctsg2          C SBE45 (S/m)
%     Stsg2          S SBE45 (psu)
%     Ttsg2_int      T of internal sensor SBE45 (C)

    %%% calculated here:
%     psealevel           pressure corrected to MSPL by NOAA equations, mb
%     qa            specific humidity of air, g/kg
%     sogE          eastward sog, kt (aka Egps)
%     sogN          northward sog, kt (aka Ngps)

%}

    ship_fields = ...
    {'ta'          ;...
    'pa'           ;...
    'psealevel'    ;...
    'paccum'       ;...
    'rh'           ;...
    'td'           ;...
    'rh2'          ;...
    'lw_dome_t'    ;...
    'lw_case_t'    ;...
    'lw_therm'     ;...
    'lw_dn'        ;...
    'sw_dn'        ;...
    'rspd'         ;...
    'rdir'         ;...
    'wspd'         ;...
    'wdir'         ;...
    'lat'          ;...
    'lon'          ;...
    'time'         ;...
    'cog'          ;...
    'sog'          ;...
    'gps_time'     ;...
    'gps_timeofday';...
    'tsea1'        ;...
    'tsea2'        ;...
    'csea1'        ;...
    'csea2'        ;...
    'ssea1'        ;...
    'ssea2'        ;...
    'sigt1'        ;...
    'sigt2'        ;...
    'flowsea1'     ;...
    'flowsea2'     ;...
    'pitch'        ;...
    'roll'         ;...
    'heave'        ;...
    'hed_ash'      ;...
    'pitch_ash'    ;...
    'roll_ash'     ;...
    'hed_gyr'      ;...
    'spdlog_v'     ;...
    'spdlog_u'     ;...
    'psealevel2';'psealevel_reported';'pa_at_zq';'qa';'sogN';'sogE'};


%% apply constants and offsets
% sig_sb = 5.67e-8;    % Stefan Boltzmann constant
% C2K = 273.15;        % temp conversion constant
Td_adj = SCS_adj(1);
Tc_adj = SCS_adj(2);

%% reference timestamp for this hour of data at 1-sec interval
% start = (ddd);
% delta = double(1.0/86400);
% last = start + 24*(3600*delta);
% reference time stamp in julian date
% jd_s = (start:delta:last)';

% make a matlab version of the full time array with 3600 sec in the hour
[the_mo, the_day] = yd2md(yyyy, ddd);

% fill time arrays

% 86400 sec in a day
sec_array = repmat((0:59)',24*60,1);

%%% 3600 min in one hour
for i = 1:60
    if i == 1
        min_array = repmat(i-1, 60, 1);
    else
        min_array = [min_array; repmat(i-1, 60, 1)];
    end
end

% 24 hours in the day
hr_array = [];
min_array_full = [];
for k = 0:23    
    hr_array = [hr_array; repmat(k,60*60,1)];
    min_array_full = [min_array_full; min_array];
end

% full matlab date and time array. 3600 members of 1-sec interval for 1 hr 86400 sec in 1 day
t_s = datenum(yyyy, the_mo, the_day, hr_array, min_array_full, sec_array);
jd_s = datenum(0, the_mo, the_day, hr_array, min_array_full, sec_array);

jd0 = md2yd_vn(yyyy,the_mo,the_day);
dayfraction = days(duration(hr_array, min_array_full, sec_array));
jd_s = jd0 + dayfraction;


% jd_s = t_s - datenum(yyyy,0,0,0,0,0);

% make a new number not susceptible to rounding, to check later for
% missing data with time_reported
time_ref = round(jd_s.*86400);


if exist(dfl) == 2
    load(dfl);


%% fields in files

% after date/time in first 3 fields, the scs files had these fields


% the_fields = {'lat';'lon'; 'cog';'sog';'hed';'Ta';'rh';'wdir';'wspd';...
%     'rdir';'rspd';'Ttsg1';'Stsg1';'pa';'sw_dn';'lw_dn';'lw_T_case';'lw_T_dome';'paccum'};
% fs = length(the_fields);
% 
% % these fields were then calculated
% full_fields = vertcat(the_fields,{'qa';'psealevel';'psealevel_reported';'sogN';'sogE'});



    %% derived data... cruise dependent

    % convert from kt to m/s if not already done so
    ship_day.sog     = ship_day.sog./1.944;
%     wspd_f  = wspd./1.944;
%     rspd_f  = rspd./1.944;

    %%% how barometric pressure corrections are handled by ship_day
    % zpim is elevation of sensor in m

    %%% how barometric pressure corrections are handled by ship_day
    % zpim is elevation of sensor in m
    pcorrection_ship    = 1013.25 *( 1 - ( 1 - zpim/44307.69231 ) ^5.253283 );
    pcorrection_noaa    = 0.125*zpim;
    % difference is 0.0743;
    ship_day.psealevel2 = ship_day.pa + pcorrection_noaa;
    ship_day.psealevel_reported = ship_day.pa + pcorrection_ship;
    ship_day.pa_at_zq = ship_day.psealevel - (0.125*zq_ship);
    % compute specific humidity
    ship_day.qa = qair_p(ship_day.ta, ship_day.rh, ship_day.pa_at_zq); 
    
    %% clean up lat/lon if necessary
%     lon(lon<PosLims(1)) = NaN;
%     lon(lon>PosLims(2)) = NaN;
%     lat(lat<PosLims(3)) = NaN;
%     lat(lat>PosLims(4)) = NaN;
    ship_day.lat = despike2(ship_day.lat); % despike & replace NaNs
    ship_day.lon = despike2(ship_day.lon);

    %% check for unreasonable values
%     imrs(imrs>2500) = NaN;
%     imrs(imrs<-20) = NaN;
%     imrl(imrl>1500) = NaN;
%     imrl(imrl<-20) = NaN;
%     paccum(paccum < 0) = NaN;    
%     ship_day.ssea1(ssea1>40) = NaN;
%     ssea1(ssea1<25) = NaN;
%     Ttsg1(Ttsg1>50) = NaN;
%     Ttsg1(Ttsg1<-2) = NaN;
%     Stsg2(Stsg2>40) = NaN;
%     Stsg2(Stsg2<25) = NaN;
%     Ttsg2(Ttsg2>50) = NaN;
%     Ttsg2(Ttsg2<-2) = NaN;
%     imts(imts>50) = NaN;
%     imts(imts<-2) = NaN;
    ship_day.cog(ship_day.cog>360) = NaN;
    ship_day.cog(ship_day.cog<0) = NaN;
%     wdir_f(wdir_f>360) = NaN;
%     wdir_f(wdir_f<0) = NaN;
%     wspd_f(wspd_f<0) = NaN;
%     wspd_f(wspd_f>50) = NaN;
    ship_day.sog(ship_day.sog>20) = NaN;
    ship_day.ta(ship_day.ta>50) = NaN;
%     prcum(prcum>100) = NaN;
%     qa(qa>25) = NaN;
    ship_day.qa(ship_day.qa<0.5) = NaN;
    ship_day.rh(ship_day.rh>100) = 100;
    ship_day.rh(ship_day.rh<0) = NaN;

    %% despike and replace NaNs with nearest neighbor
    ship_day.cog = unwrap(ship_day.cog*pi/180); % unwrap and convert to radians
    [ship_day.cog,~] = despike2(ship_day.cog);

    ship_day.hed = unwrap(ship_day.hed_gyr*pi/180); % unwrap and convert to radians
    [ship_day.hed,~] = despike2(ship_day.hed);
    
    %%% 2D sonic on mast
    [ship_day.wspd,~] = despike2(ship_day.wspd);
    [ship_day.rspd,~] = despike2(ship_day.rspd);
    ship_day.wdir = unwrap(ship_day.wdir*pi/180); % unwrap and convert to radians
    [ship_day.wdir,~] = despike2(ship_day.wdir);
    ship_day.rdir = unwrap(ship_day.rdir*pi/180); % unwrap and convert to radians
    [ship_day.rdir,~] = despike2(ship_day.rdir);

%     %%% just despike the rest... why do this on some of the vars and not
%     %%% others?
%     [psealevel,~] = despike2(psealevel);
%     [qa,~] = despike2(qa);

    %%% ship_day speed components
    ship_day.sogN = ship_day.sog.*cos(ship_day.cog);
    ship_day.sogE = ship_day.sog.*sin(ship_day.cog);
    
    %%% convert heading, cog, and true wind directions back to deg: 0 to 360
    ship_day.hed = (ship_day.hed*180/pi);
    ship_day.hed = mod(ship_day.hed + 360, 360);
    ship_day.cog = (ship_day.cog*180/pi);
    ship_day.cog = mod(ship_day.cog + 360, 360);
    ship_day.wdir = (ship_day.wdir*180/pi);
    ship_day.wdir = mod(ship_day.wdir + 360, 360);
    
    %%% convert rdir to degrees: +/- 180
    ship_day.rdir = (ship_day.rdir*180/pi);
    ship_day.rdir(ship_day.rdir>180) = ship_day.rdir(ship_day.rdir>180)-360; 

%     %% format output structure called ship_day with all reported data, derived data, and matlab & jd time fields
%   
%     %%% some time entries are repeated. Make sure they don't count.
%     [time_unique_0, ID_unique_0] = unique(ship_day.gps_time);
%     
%     %%% throw out data that exceeds or comes before t_s; 
%     %%% sometimes this happens because the time stamp / logger slips forward
%     overtime = find(time_unique_0 > time_ref(end));
%     undertime = find(time_unique_0 < time_ref(1));
%     badtimes = [undertime overtime];
%     time_unique = time_unique_0;
%     ID_unique = ID_unique_0;
%     if ~isempty(badtimes)
%         time_unique(badtimes) = [];
%         ID_unique(badtimes) = [];
%     end
% 
%     %%% fill missing values for full 3600 seconds in the hour
%     not_reported = find(ismember(time_ref, time_reported) == 0);
%     reported = find(ismember(time_ref, time_reported) == 1);
%     
% %     clear ship_day;
% %     %%% dropping _s suffix for time varibales in final structure
% %     ship_day.jd = jd_s;
% %     ship_day.t = t_s;
% 
%     if length(reported) ~= length(ID_unique)
%      disp(['*** n(reported)= ' sprintf('%i',length(reported)) ' ... but n(ID_unique)= ' sprintf('%i',length(ID_unique))]);
%     end
%     % this is for everything but t and jd... since the vars with those
%     % prefixes already exist in this program and have been added without
%     % the prefixes to the structure already.
% %     for i = 1:length(full_fields)
% %         eval([full_fields{i} '_s = nan(3600*24,1);']);
% %         eval([full_fields{i} '_s(reported) = ' full_fields{i} '(ID_unique);']);
% %         eval(['ship_day.' full_fields{i} ' = ' full_fields{i} '_s;']);
% %     end
%     
% %     disp('size of ship_day dataset: time x fields');
% %     disp(size(ship_day));
else
    
    %%% create empty structure for this hour
    disp('file does not exist for this day');
        
    clear ship_day;
    ship_day.jd = jd_s;
    ship_day.t  = t_s;
        

    
    for i = 1:length(ship_fields)
        eval(['ship_day.' ship_fields{i} ' = nan(24*3600,1);']);
    end
    
    
end %%% end of file i/o check and read loop

end  %%% end of function

