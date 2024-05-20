function out = read_wxt(dfl3,ddd,hhh,zwxt)
%{
reads wxt0 files and returns 1 min avg array with the
following columns:

    1  jd_ref          decimal day of year
    2  rwdir_wxt       rel wind direction ,0-360 deg from bow
    3  rwspd_wxt       rel wind speed, m/s
    4  Tair_wxt        air temperature, C
    5  RH_wxt          RH, %
    6  Pmb_wxt         air pressure, mb or hPa
    7  rain_wxt        rain rate, mm/hr
    8  run_wxt          relative wind component: positive south-to-north (v from normal convention), WXT
    9  ruw_wxt          relative wind component: positive east-to-west (-u from normal convention), WXT

    input parameters: dfl3 = string path to wxt0 file
                       ddd = julian date
                        hh = hour

    wxt0 file format:
    0000814,  227,  1.7,  28.9,  87.0,  1004.3,  0.00,  0,  0.0,  0.0,  0,  0.0
    1 2 3     4     5     6      7      8        9     10   11    12   13   14
    1:3 PSD timecode: MMSSsss
    4   relative wind direction, degrees
    5   relative wind speed, m/s
    6   air temperature, C
    7   RH, %
    8   pressure, mb or hPa
    9   cumulative rain, mm
    10  rain event duration, sec
    11  rain intensity, mm/hr
    12  cumulative hail, hits/cm^2
    13  hail event duration, s
    14  hail intensity, hits/cm^2/hr

%}

%% initialize

g = 0.0098;          % adiabatic lapse rate
sig_sb = 5.67e-8;    % Stefan Boltzmann constant
C2K = 273.15;        % temp conversion constant

nr3 = 14;     % number of variables in wxt0
raw = [];
temp = {};

%% read wxt0 file
disp(['Reading wxt0 file for hour ',int2str(hhh)]);
flist = fopen(dfl3,'r');

temp = textscan(flist,['%2f%2f%3f ',repmat('%f ',1,11),'%*[^\n]'],'delimiter',',');
raw = [raw cell2mat(temp)];

[wxt_run,wxt_ruw] = sd_to_uv(raw(:,5),raw(:,4));  % compute relative wind components ru and rv
raw = [raw,wxt_run,wxt_ruw];                      % add columns 15,16 to array

jd_pc_3 = ddd+(hhh+(raw(:,1)+(raw(:,2)+raw(:,3)/1000)/60)/60)/24;   % timestamp from PSD system
fclose(flist);

%% average wxt to 1 min timestamp
start = (ddd+hhh/24);
delta = double(1.0/1440);
last = start + 59*delta;
jd_bin = start:delta:last+delta;

% average all columns
wxt = interval_avg(jd_pc_3,raw(:,4:16),jd_bin');
[wspd,wdir] = uv_to_sd(wxt(:,13),wxt(:,14));   % compute wspd, wdir from averaged rUn, rUw components of measured rspd, rdir
wxt(:,2) = wdir;
wxt(:,3) = wspd;
wxt(:,6) = wxt(:,6) + 0.125*zwxt;            % sea level P
out = [wxt(:,1:6),wxt(:,9),wxt(:,13:14)];
end

