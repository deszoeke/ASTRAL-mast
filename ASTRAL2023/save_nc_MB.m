%%% program to save the final netcdf file for PISTON 2018. 

% Run these first:
% PISTON_2018_save_mat.m ... which takes advantage of the most important part of fix_met_sea.m
% the latter of which examined whether corrections were needed for the sea snake, 
% TSG T, or TSG S. Since only the TSG S needed adjusting, I copied what was needed to this
% program and so it performs the fix directly. 

% check the output of this file, prior files, and that for 2019 as well in
% PISTON_check_nc.m ... which is in Documents/MATLAB/

% this program calls assign_nc_info.m which is very minimally modified from
% experiment to experiment, ideally. The only things that are still hard
% coded are the corrections applied to seawater data. The rest of the
% instrument, height, cruise info are automatically loaded where needed in
% setup_cruise.m

% this program calls coare because we had to recreate some fields to
% converge on the final product between Elizabeth and Ludovic's files. 

% August 2022 EJT


%% Initialize run parameters
close all;
fclose all;
clear all;
warning ('off','MATLAB:MKDIR:DirectoryExists');
setup_cruise_MB;

% system specific path defs
sysType = computer;
username=char(java.lang.System.getProperty('user.name'));
if strncmp(sysType,'MACI64',7)     % set Mac paths
    data_drive = '/Users/eliz/DATA/';
    path_prog = fullfile(data_drive,'PISTON_2018',ship,'Scientific_analysis','programs');
elseif strncmp(sysType,'PCWIN64',7)  % set PSL DAC paths
    data_drive = 'D:\DATA\';
    path_prog = fullfile(data_drive,'PISTON_2018',ship,'Scientific_analysis','programs');
end

% matlab script path
restoredefaultpath
cd(fullfile(path_prog,'save_nc'));
% cd('/Users/eliz/DATA/PISTON_2018/TGT/scientific_analysis/programs/save_nc/subs/');
addpath(genpath(fullfile(path_prog,'save_nc','subs')));
addpath('/Users/eliz/DATA/PISTON_2019/Sally_Ride/scientific_analysis/programs/flux/read_map_fix_print_plot/');
addpath(genpath(fullfile(path_prog,'save_nc')));
rehash toolboxcache;


% where both input and output data live
indir = '/Users/eliz/DATA/MISOBOB/';
outdir = '/Users/eliz/DATA/MISOBOB/';
plotdir = '/Users/eliz/DATA/MISOBOB/nccheck/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 min, 10 min, hourly input mat files:
% vin = 'r2'; 
% f1  = [cruise  '_1min_' vin '.mat'];
% f10 = [cruise '_10min_' vin '.mat'];
% f60 = [cruise '_60min_' vin '.mat'];
% 
% load([indir f1]);
% load([indir f10]);
% load([indir f60]);

load('/Users/eliz/DATA/MISOBOB/PISTON_MISOBOB_2018_flux_10_v4_decorr_SdS_allvars.mat');
% f10 = table2struct(flux10);

% output files we are making here:
vnum = 'R1'; 
product_version_nc = '1'; % used to fill in global attributes automatically

fillval = NaN;

%% load

clear a c e60;
f10 = table2struct(flux10_all,"ToScalar",true);
f10.t = f10.datenum;
tres = '10 minutes';
nt = length(f10.t);

f10 = rmfield(f10,'sog2');

%% REMOVE EEZ and come up with resulting filenames
% EEZ data 
eez_start  = datenum(yr,6,4,23,20,0);
eez_in     = datenum(yr,6,16,0,0,0); 
eez_out    = datenum(yr,6,30,0,0,0);
eez_end    = datenum(yr,7,20,0,0,0);

eez_start_jd  = eez_start - datenum(yr,0,0,0,0,0);
eez_in_jd     = eez_in - datenum(yr,0,0,0,0,0);
eez_out_jd    = eez_out - datenum(yr,0,0,0,0,0);
eez_end_jd    = eez_end - datenum(yr,0,0,0,0,0);


% % check for uniqueness of a.t ... good news! all unique!
% [t_un, ID_un] = unique(a.t);

% take out EEZ if still present
wh_bad1 = find(f10.t < eez_start);
wh_bad2 = find(f10.t > eez_in & f10.t < eez_out); % between leg 2 and leg 3
wh_bad3 = find(f10.t > eez_end); % after leg 3
wh_bad = [wh_bad1; wh_bad2; wh_bad3];
ID_bad = zeros(nt,1);
ID_bad(wh_bad) = 1;

fmet = fields(f10);
%%% if you want to remove bad data altogether
for j = 1:length(fmet)
    x10 = f10.(fmet{j});
    x10(wh_bad) = [];
    a.(fmet{j}) = x10; 
end

a.jd = a.t - datenum(yr,0,0,0,0,0);
nt = length(a.t);

day1_str = datestr(a.t(1),'yyyymmdd');
dayn_str = datestr(a.t(end),'yyyymmdd');

ncname10=['MISOBOB-nav-met-sea-flux-10min_' ship '_' day1_str '_thru_' dayn_str '_' vnum '.nc'];
ncnamehr=['MISOBOB-nav-met-sea-flux-1hour_' ship '_' day1_str '_thru_' dayn_str '_' vnum '.nc'];

%% define new vars

% time elapsed in seconds since Jan 1 of that year
% note... don't add this to a structure. just call it time. used below in a
% special netcdf way. 

time = (a.t - datenum(yr,1,1,0,0,0))*24*3600;
ndays = a.jd(end) - a.jd(1);
zsea_ship = zsea1_ship;

[a.year, a.month, a.day, a.hour, a.minute, ~] = datevec(a.t);

a.cq = a.cql;
a.jmanuv = a.manuv;
a.co2_std_lic = a.std_co2_lic;
a.qair_std_lic = a.std_qair_lic;
% a.ssea_ship = ones(nt,1)*33.0;
% a.uvar = nan(nt,1);
% a.vvar = nan(nt,1);
% a.wvar = nan(nt,1);
% a.tvar = nan(nt,1);
a.rdir(a.rdir>180) = a.rdir(a.rdir>180)-360;    % +/- 180 deg format
a.paccum = cumsum(a.prate)/6; % 6 X 10 min segments per hour
a.ssea_ship(a.ssea_ship < 28) = NaN;

%%% fix or fill bad/missing rh values from vaisala with corrected ship data
plot_tsg = 0;
if plot_tsg == 1
    figure;
    plot(a.t, a.rhair_ship+nanmedian(a.rhair - a.rhair_ship),'o',a.t, a.rhair,'x'); 
    grid on; datetick('x','mm/dd','keeplimits');legend('ship corrected','psl');
    hold on;
end

rhair_shipC = a.rhair_ship+nanmedian(a.rhair - a.rhair_ship);
a.rhair_ship = rhair_shipC;
rh_best  = a.rhair;
bad_rh = find(isnan(a.rhair) == 1);
rh_best(bad_rh) = a.rhair_ship(bad_rh);

% plot(a.t, rh_best,'.k');
a.rhair = rh_best;


% ship manuevering and sampling the ship plume. 1 = bad; 0 = good;
jship = ones(nt, 1); % assume all are bad (1)
jship(a.jmanuv < 0.5 | a.jplume < 0.5) = 0; % data are good (0) where jmanuv < 0.5 (3 is max), jplume = 0
a.flag_bad_ship = jship;

% wind from aft.  1 = bad; 0 = good;
jbulk = ones(nt,1); % assume all are bad (1)
jbulk(abs(a.rdir) < 140) = 0; % all data with rdir within +/- 140 deg are good (0)
a.flag_bad_bulk = jbulk;

%% RERUN COARE
B = coare36vnWarm_et(a.jd, a.wspd, zu, a.tair, zt, a.rhair, zq,...
             a.psealevel, a.tsea, a.sw_down, a.lw_down, a.lat, a.lon, 100,...
             a.prate, zsnk, a.ssea_ship, nan, nan, 2, 2, 2);

cfields = {'ustar';'tau_bulk';'hs_bulk';'hl_bulk';'hb_bulk';'hb_son';'hl_webb';'tstar';'qstar';...
    'rough_u';'rough_t';'rough_q';...
    'cd';'ch';'ce';'mo_length';'zeta';'dt_skin';'dq_skin';'dz_skin';'wspd_2';'tair_2';'qair_2';...
    'rhair_2';'wspd_2N';'tair_2N';'qair_2N';'lw_net';'sw_net';'Le';'rhoair';...
    'wspd_N';'wspd_10';'wspd_10N';'cd10N';'ch10N';'ce10N';...
    'hrain';'qskin';'erate';'tair_10';'tair_10N';'qair_10';'qair_10N';'rhair_10';'pair_10';...
    'rhoair_10';'gust';'wave_whitecap_frac';'wave_edis';...
    'dt_warm';'dz_warm';'dt_warm_to_skin';'du_warm'};

for i = 1:length(cfields)
    eval(['a.' cfields{i} '   = B(:,i);']);
end

% interface
a.tskin       = a.tsea + a.dt_warm_to_skin - a.dt_skin;

% recalculate upwelling radiative fluxes
a.lw_up       = a.lw_net - a.lw_down;
a.sw_up       = a.sw_net - a.sw_down;

% net heat flux: heating into the ocean (use minus sign if fluxes haven't
% been flipped yet!)
a.hnet        = a.sw_net + a.lw_net - a.hs_bulk - a.hl_bulk - a.hrain;
[a.tau_bulk_u, a.tau_bulk_v] = sd_to_uv_met(a.tau_bulk, a.wdir);
[a.u_10N, a.v_10N] = sd_to_uv_met(a.wspd_10N, a.wdir);
[a.u, a.v] = sd_to_uv_met(a.wspd, a.wdir);

%%% option: sign changes of COARE fields... the bulk sensible, latent, and rain
%%% fluxes are defined positive by COARE, but we don't provide them to others that way
flips = {'hs_bulk';'hl_bulk';'hrain'};
for k = 1:length(flips)
    eval(['a.' flips{k} ' = - a.' flips{k} ';']);
end

a = orderfields(a);

%% produce hourly averages for ludovic's big flux database

% some are plotted, some are saved, but should save all in the future

% hourly time base
jd_h_bin = floor(a.jd(1)) : (1.0/24.0) : ceil(a.jd(end)); % hourly bin edges
jd_h_bin = jd_h_bin'; % hourly bin edges
eez_h = find(jd_h_bin > eez_in_jd & jd_h_bin < eez_out_jd);
jd_h_bin(eez_h) = [];
e60.jd = jd_h_bin(1:end-1); % just days
e60.t = datenum(yr,0,0,0,0,0) + e60.jd;

[e60.year, e60.month, e60.day, e60.hour, e60.minute, ~] = datevec(e60.t);
e60.minute = e60.minute*0.0;
e60.t = datenum(e60.year, e60.month, e60.day, e60.hour, e60.minute, 0);
e60.jd = e60.t - datenum(yr,0,0,0,0,0);
jd_h_bin = [e60.jd; e60.jd(end) + (1./24)];

var_avg = fields(a);

%%% perfrom a straight average, just not the time fields;
for i = 1:length(var_avg)
    if strcmp(var_avg{i},'jd') ~=1 && strcmp(var_avg{i},'t') ~=1 
    	eval(['e60.' var_avg{i} ' = interval_avg_var(a.jd, a.' var_avg{i} ', jd_h_bin);']);
    end
end

[e60.year, e60.month, e60.day, e60.hour, e60.minute, ~] = datevec(e60.t);

%%  special values, directions, and vectors
 
% e60.minute(isfinite(e60.minute)==1) = 0;
e60.tilt = interval_avg_dir(a.jd,a.tilt,jd_h_bin); 

% [e60.cspd, e60.cdir]            =  interval_avg_vect(a.jd, a.cspd, a.cdir, jd_h_bin);
% [e60.wspd, e60.wdir]            =  interval_avg_vect(a.jd, a.wspd, a.wdir, jd_h_bin);
% [e60.wspd_sfc, e60.wdir_sfc]    =  interval_avg_vect(a.jd, a.wspd_sfc, a.wdir_sfc, jd_h_bin);

%%% check to make sure this is resulting in +/- 180 for all
% rdir_new, rspd_new2, rspd_raw
% [e60.rspd, e60.rdir]            =  interval_avg_rvect(a.jd, a.rspd, a.rdir, jd_h_bin);
 
plot_rdir = 0;
if plot_rdir == 1
   figure;
   plot(a.t, a.rdir,'x',e60.t, e60.rdir,'o'); 
   legend('10-min','hr');    
   title('rdir check');
end

%% fix hed, cog, sog if needed
% These are equivalent to the ones below, and so are not needed. 

% hr2.sog = sqrt(e60.sogN.^2 + e60.sogE.^2);  % recompute sog from speed components
% 
% hr2.cog = atan2(e60.sogE, e60.sogN)*r2d;  % recompute cog from speed components
% hr2.cog = mod(e60.cog + 360, 360); % reorient units2
% 
% hr2.heading = atan2(e60.hedE, e60.hedN)*r2d; % recompute avg hed from avg N/E components
% hr2.heading = mod(e60.hed + 360, 360); % reorient units    

%% fix any negative angles (NOT the rdirs... which are now +/- 180 on purpose
% is this necessary? there weren't any bad angles in PISTON 2019
angles1 = {'wdir'};
angles2 = {'cog';'hed'};
for i = 1:length(angles1)
    this = e60.(angles1{i});
    bad_ang = find(this<0);
    if ~isempty(bad_ang)
        disp(['bad angles found for ' angles1{i}]);
        this(bad_ang) = this(bad_ang)+360;
        e60.(angles1{i}) = this;
    end
end


%%  special special or coare recalculations after averaging

temp = e60.prate;
temp(isnan(temp)) = 0;
e60.paccum = cumsum(temp)/1; % 1 x 60 min segment per hour... already in increments of hours

C = coare36vnWarm_et(e60.jd, e60.wspd, zu, e60.tair, zt, e60.rhair, zq,...
             e60.psealevel, e60.tsea, e60.sw_down, e60.lw_down, e60.lat, e60.lon, 600,...
             e60.prate, zsnk, e60.ssea_ship, nan, nan, 2, 2, 2);
cfields = {'ustar';'tau_bulk';'hs_bulk';'hl_bulk';'hb_bulk';'hb_son';'hl_webb';'tstar';'qstar';...
    'rough_u';'rough_t';'rough_q';...
    'cd';'ch';'ce';'mo_length';'zeta';'dt_skin';'dq_skin';'dz_skin';'wspd_2';'tair_2';'qair_2';...
    'rhair_2';'wspd_2N';'tair_2N';'qair_2N';'lw_net';'sw_net';'Le';'rhoair';...
    'wspd_N';'wspd_10';'wspd_10N';'cd10N';'ch10N';'ce10N';...
    'hrain';'qskin';'erate';'tair_10';'tair_10N';'qair_10';'qair_10N';'rhair_10';'pair_10';...
    'rhoair_10';'gust';'wave_whitecap_frac';'wave_edis';...
    'dt_warm';'dz_warm';'dt_warm_to_skin';'du_warm'};

for i = 1:length(cfields)
    if contains(cfields(i),'hb_son') ~= 1 && contains(cfields(i),'zeta') ~= 1 && ...
        contains(cfields(i),'dq_skin') ~= 1 && ...
        contains(cfields(i),'Le') ~= 1 && ... 
        contains(cfields(i),'wspd_N') ~= 1 && ...
        contains(cfields(i),'du_warm') ~= 1 && ...
        contains(cfields(i),'tair_2N') ~= 1 && ...
        contains(cfields(i),'qair_2N') ~= 1 && ...
        contains(cfields(i),'dt_warm_to_skin') ~= 1 && ...
        contains(cfields(i),'dz_warm') ~= 1 && ...
        contains(cfields(i),'dt_warm') ~= 1 && ...
        contains(cfields(i),'du_warm') ~= 1
        %% only rewrite variables if they aren't what was just written out
        eval(['e60.' cfields{i} '   = C(:,i);']);
    end
end

% interface
e60.tskin       = e60.tsea + e60.dt_warm_to_skin - e60.dt_skin;

% recalculate upwelling radiative fluxes
e60.lw_up       = e60.lw_net - e60.lw_down;
e60.sw_up       = e60.sw_net - e60.sw_down;

% net heat flux: heating into the ocean (use minus sign for this if the
% signs for hs, hl, and hrain haven't been flipped yet!)
e60.hnet        = e60.sw_net + e60.lw_net - e60.hs_bulk - e60.hl_bulk - e60.hrain;
[e60.tau_bulk_u, e60.tau_bulk_v] = sd_to_uv_met(e60.tau_bulk, e60.wdir);
[e60.u_10N, e60.v_10N] = sd_to_uv_met(e60.wspd_10N, e60.wdir);

% for some reason the warm layer looks weird at hourly levels... so just
% average 10-min output instead

e60.dt_warm_to_skin = interval_avg_var(a.jd, a.dt_warm_to_skin, jd_h_bin);
e60.dt_warm = interval_avg_var(a.jd, a.dt_warm, jd_h_bin);
e60.dz_warm = interval_avg_var(a.jd, a.dz_warm, jd_h_bin);

%% option: sign changes of COARE fields... the bulk sensible, latent, and rain
%%% fluxes are defined positive by COARE, but we don't provide them to others that way
%%% **** don't need to do id or cov fluxes again because they are just
%%% averages of the already flipped 10-min data
flips = {'hs_bulk';'hl_bulk';'hrain'};
for k = 1:length(flips)
    eval(['e60.' flips{k} ' = - e60.' flips{k} ';']);
end

%% do these separately because of wrapping issues
cog_10r = unwrap(a.cog*pi/180); % unwrap and convert to radians
cog_60r = interval_avg_var(a.jd, cog_10r, jd_h_bin);
cog_60 = (cog_60r*180/pi);
cog_60 = mod(cog_60 + 360, 360);
e60.cog = cog_60;

[hedN_10, hedE_10] = pol2cart(a.heading*d2r,1);  % N/E heading components in radians
hedN_60 = interval_avg_var(a.jd, hedN_10, jd_h_bin);
hedE_60 = interval_avg_var(a.jd, hedE_10, jd_h_bin);

hed_60 = atan2(hedE_60, hedN_60)*r2d; % recompute avg hed from avg N/E components
hed_60 = mod(hed_60 + 360, 360);
e60.heading = hed_60;

[rUn_10, rUw_10] = sd_to_uv(a.rspd, a.rdir);
rUn_60 = interval_avg_var(a.jd, rUn_10, jd_h_bin);
rUw_60 = interval_avg_var(a.jd, rUw_10, jd_h_bin);

[Un_10, Uw_10] = sd_to_uv(a.wspd, a.wdir);
Un_60 = interval_avg_var(a.jd, Un_10, jd_h_bin);
Uw_60 = interval_avg_var(a.jd, Uw_10, jd_h_bin);

% [Uns_10, Uws_10] = sd_to_uv(a.wspd_sfc, a.wdir_sfc);
% Uns_60 = interval_avg_var(a.jd, Uns_10, jd_h_bin);
% Uws_60 = interval_avg_var(a.jd, Uws_10, jd_h_bin);
% 
% [Unc_10, Uwc_10] = sd_to_uv(a.cspd, a.cdir);
% Unc_60 = interval_avg_var(a.jd, Unc_10, jd_h_bin);
% Uwc_60 = interval_avg_var(a.jd, Uwc_10, jd_h_bin);

% [rspd_60,rdir_60] = uv_to_sd(rUn_60,rUw_60);% compute true wspd, wdir, and components Un and Uw from rUn, rUw components of adjusted components of rspd, rdir, and hed + cog + sog
% e60.rdir = rdir_60;
% e60.rspd = rspd_60;
% e60.rdir(e60.rdir>180) = e60.rdir(e60.rdir>180)-360;    % +/- 180 deg format

[~,wdir_60] = uv_to_sd(Un_60,Uw_60);% compute true wspd, wdir, and components Un and Uw from rUn, rUw components of adjusted components of rspd, rdir, and hed + cog + sog
e60.wdir = wdir_60;

% [~,wdirs_60] = uv_to_sd(Uns_60,Uws_60);% compute true wspd, wdir, and components Un and Uw from rUn, rUw components of adjusted components of rspd, rdir, and hed + cog + sog
% e60.wdir_sfc = wdirs_60;
% 
% [~,cdir_60] = uv_to_sd(Unc_60,Uwc_60);% compute true wspd, wdir, and components Un and Uw from rUn, rUw components of adjusted components of rspd, rdir, and hed + cog + sog
% e60.cdir = cdir_60;

%% fix integers data flags... set to 0 if nan or decimal place after avg-ing
ints = {'flag_bad_ship';'flag_bad_bulk';'jplume'};

for i = 1:length(ints)
    this = e60.(ints{i});
    bad_int = find(isnan(this) == 1 | mod(this,1) ~= 0);
    if ~isempty(bad_int)
        % could try this to save more data
        % this(bad_int) = round(this(bad_int)); 
        this(bad_int) = 0; % just says anything not perfect is bad
        e60.(ints{i}) = this;
    end
    %%% also remove eez in case that got messed up just now

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e60 = orderfields(e60);


%% for 10-min and 1-hour files, assign vars to new structure

for k = 1:2

    
    if k == 1
        thencname = ncname10;
    elseif k == 2
        thencname = ncnamehr;
    end

%%% time, navigation
vars1 = {'lat';'lon';...
'cog';'sog';'heading';...
'flag_bad_ship';'flag_bad_bulk';...
'year';'month';'day';'hour';'minute';};

% excluded for this particular cruies:

%   allll the wxt variables... bad quality
%    wspd_wxt';'wdir_wxt';'rspd_wxt';'rdir_wxt';...
%        'Ta_wxt';'qa_wxt';'rh_wxt';'rhoa_wxt';'H2O_wxt';...
%        'prate_wxt';'paccum_wxt';'psealevel_wxt
%   ship variables... they confuse people
%       'tair_ship'; ;'psealevel_ship' 'wspd_ship' 'wdir_ship';
%       'sw_up_ship'; 'lw_up_ship' 'sw_down_ship'; 'lw_down_ship'
%       rhair_ship, prate_ship,
%   ship variables... that were actually bad on this cruise

%%% met and seawater, including warm layer and cool skin
vars2 = {'tair';'tair_2';'tair_10';'tair_10N';...
'rhair';'rhair_2';'rhair_10';...
'qair';'qair_2';'qair_10';'qair_10N';...
'rhoair';'rhoair_10';...
'prate';'paccum';'psealevel';'pair_10';...
'rdir';'rspd';'wdir';...
'wspd';'wspd_2';'wspd_2N';'wspd_10';'wspd_10N';...
'u';'v';'u_10N';'v_10N';...
'tskin';'qskin';'tsea';'qsea';'tsea_ship';'ssea_ship';...
'dt_skin';'dz_skin';'dt_warm';'dz_warm';'dt_warm_to_skin';...
'wave_edis';'wave_whitecap_frac';...
'sw_down';'sw_down_clear';'sw_up';...
'lw_down';'lw_down_clear';'lw_up';};

%%% fluxes
vars3 = {'hnet';'hs_bulk';...
'hs_cov';'hs_id';
'hl_bulk';'hl_webb';
'hl_cov';'hl_id';
'hrain';'hb_bulk';...
'tau_cov';'tau_cov_cross';'tau_id';...
'tau_bulk';'tau_bulk_u';'tau_bulk_v';'tilt';'mo_length';...
'ustar';'tstar';'qstar';'cd';'ce';'ch';'cd10N';'ce10N';'ch10N';...
'erate';'rough_u';'rough_t';'rough_q';'gust'};

%%% more stuff for Ludovic's big ship database
vars4 = {'tair_ship';'psealevel_ship';'wspd_ship';'wdir_ship';...
       'sw_down_ship';'lw_down_ship';'rhair_ship';...
       'ct';'cq';'cu';'cw';...
       'qair_lic';'qair_std_lic';'co2_lic';'co2_std_lic';...
       'jmanuv';'jplume';'uvar';'vvar';'wvar';'tvar';'qvar'};


vars = vertcat(vars1, vars2, vars3, vars4);

% make new array c for the variables we want to save from a (input) 
disp(['var = time ...  n(NaN) = ' sprintf('%i',length(find(isnan(time) == 1)) )]);

for i = 1:length(vars)
    
    if k == 1
        eval(['c.' vars{i} ' = a.' vars{i} ';']);
    elseif k == 2
        eval(['c.' vars{i} ' = e60.' vars{i} ';']);
    end
        
    %=========== reset nans to fill value, except for time, position, flags
    clear bad_data thisvar;

    eval(['thisvar = c.' vars{i} ';']);
    bad_data = find(isnan(thisvar)==1);
    disp(['var = ' vars{i} ' ... n(NaN) = ' sprintf('%i',length(bad_data)) ]);
    if strcmp(vars{i},'t') == 0 && strcmp(vars{i},'jd') == 0 && ...
        strcmp(vars{i},'lat') == 0 && strcmp(vars{i},'lon') == 0 && ...
        contains(vars{i},'flag') == 0
%         strcmp(vars{i},'year') == 0 && strcmp(vars{i},'month') == 0 && ...
%         strcmp(vars{i},'day') == 0 && strcmp(vars{i},'hour') == 0 && ...
%         strcmp(vars{i},'minute') == 0 
            thisvar(bad_data) = fillval;
            eval(['c.' vars{i} ' = thisvar;']);
    end

end

thefields = fields(c);
nf = length(thefields);

%% assign netcdf variable attributes
[units, standard_names, long_names, comments, types, heights,...
    instruments, locations, methods] = ...
    assign_nc_info_MB(thefields, nf, ins, ins_ship);
for j = 1:nf
    if isempty(units{j}) == 1
        disp(['missing unit for field ' sprintf('%i',j) ' = ' string(thefields(j))]);
    end
    if isempty(standard_names{j}) == 1
        disp(['missing standard name for field ' sprintf('%i',j) ' = ' string(thefields(j))]);
    end
    if isempty(long_names{j}) == 1
        disp(['missing longname for field ' sprintf('%i',j) ' = ' string(thefields(j))]);
    end    
    if isempty(types{j}) == 1
        disp(['missing type for field ' sprintf('%i',j) ' = ' string(thefields(j))]);
    end    
end

%%% write a table when 10-min file is being produced since that's when
%%% fluxes will be computed and listed 
T = table(thefields, units, standard_names, long_names, comments, types, heights, instruments, locations, methods);
writetable(T, [indir cruise '_data_info.csv'],'Delimiter',',','QuoteStrings',true); 

%% readmes
make_readme = 0;
if make_readme == 1

readme_data_f10     = [pathreadme '_fields_' cruise '_' version_num '.txt'];
header_data_f10     = [pathreadme '_header_' cruise '_' version_num '.txt'];
fileID_d_f10 = fopen(readme_data_f10,'w');
fileID_h_f10 = fopen(header_data_f10,'w');

for j = 1:nfa
    formatSpec = 'i(%i) = %15s %-60s\n';
    fprintf(fileID_d_f10, formatSpec, j, string(thefields(j)), string(lf(j)) );
    fprintf(fileID_h_f10,'%20s',string(thefields(j)) );
end
fclose(fileID_d_f10);
fclose(fileID_h_f10);

end % if making readme


%% netcdf file construction

if k == 2
    time = (e60.t - datenum(yr,1,1,0,0,0))*24*3600;
end

thedata = c;
vartypes = cell(nf,1);
dr = nan(nf,1);
dc = nan(nf,1);
% note: varname is an internal netcdf library reserved special name. do not use.
for i = 1:nf
    eval(['thisvar = thedata.' thefields{i} ';']);
    thevar = char(thefields(i));
    vartype = class(thisvar);
    vartypes(i) = {vartype};
    [dr(i), dc(i)] = size(thisvar);
%     disp([' saving ' thevar ' = ' vartype]);
end

dt = find(dr == length(time) & dc == 1);   % discretized in time
dn = find(dr > 2 & dr < 100 & dc == 1);          % some other short array or short readme
d1 = find(dr == 1 & dc == 1);                     % 1D fields
nd1 = length(d1);
ndn = length(dn);
ndt = length(dt);
% nd =  ndt;
disp([' out of ' sprintf('%i',nf) ' total vars, ' sprintf('%i',ndt) ' are accounted for and time discretized']);

% initialize file:
thenc = netcdf.create([outdir thencname], 'NETCDF4');

% define the dimensions

%%% for multiple trajectories
% [obs_dimlen,c]=size(jdy);
% trajectory_dimlen=1;  %number of cruise legs
% obs_dimID = netcdf.defDim(ncid,'obs',obs_dimlen); %simply length of data set
% traj_dimID= netcdf.defDim(ncid,'trajectory',trajectory_dimlen); %for single trajectory this dimension could be omitted but keep it for easiness and consistency for cruises with multiple legs
% dimIDs = [traj_dimID, obs_dimID];%used to pass the dimids of the dimensions of the NETCDF variables. All the NETCDF variables we are creating share the same dimensions.

ntime=length(time);

% create ID for (define) the dimensions
timedimID=netcdf.defDim(thenc,'time',ntime);

% time vars must be in seconds since first day of year. Must be double or float.
timeID=netcdf.defVar(thenc,'time','double',timedimID);
netcdf.putAtt(thenc,timeID,'units',['seconds since ' sprintf('%i',yr) '-01-01 00:00 UTC']);
netcdf.putAtt(thenc,timeID,'standard_name','time');
netcdf.putAtt(thenc,timeID,'long_name','time');
netcdf.putAtt(thenc,timeID,'comment','');
netcdf.putAtt(thenc,timeID,'method',...
   ['this time marks the beginning of the averaging interval, i.e. this time ' ...
    'is the leading bin edge for each bin average. values provided are averages '...
    'of all samples from this time step to the next time step. Only 1 good data '...
    'point was required for a valid average to be computed and reported. Data '...
    'were originally collected at intervals of 5 min (skin ocean temp ROSR), ' ...
    '1 min (radiation, temp, humidity, pressure, rain, sea water, extra met sensor), ' ...
    '10 Hz (wind, fast humidity, GPS, heading, pitch/roll, motion)']);
netcdf.putAtt(thenc, timeID,'coverage_content_type', 'coordinate');
netcdf.defVarFill(thenc,timeID,false,fillval);
netcdf.putAtt(thenc,timeID, 'axis', 'time'); % define that these are time (T) axis variables
netcdf.putVar(thenc,timeID,time); % actually put the variable in

theID = nan(ndt,1);
clear j;

%% create ID (define) for the rest of the time-based variables
for i = 1:ndt
    j = dt(i);
    theID(j) = netcdf.defVar(thenc,char(thefields(j)),char(vartypes(j)),timedimID);
end

% define the variable attributes
for i = 1:ndt

    netcdf.putAtt(thenc,theID(i),'units',               char(units{i}));

    if isempty(standard_names{i}) == 0  
        netcdf.putAtt(thenc,theID(i),'standard_name',   char(standard_names{i}));
    end
    
    netcdf.putAtt(thenc,theID(i),'long_name',           char(long_names{i}));
    
    if isempty(methods{i}) == 0 
        netcdf.putAtt(thenc,theID(i),'method',          char(methods{i}));
    end
    if isempty(comments{i}) == 0 
        netcdf.putAtt(thenc,theID(i),'comment',         char(comments{i}));
    end
    if isempty(heights{i}) == 0 
        netcdf.putAtt(thenc,theID(i),'height',          heights{i});
    end
    
    netcdf.putAtt(thenc,theID(i),'cdm_data_type',       char(types{i}));
    
    if isempty(instruments{i}) == 0 
        netcdf.putAtt(thenc,theID(i),'instrument',      char(instruments{i}));
    end
    if isempty(locations{i}) == 0 
        netcdf.putAtt(thenc,theID(i),'location',        char(locations{i}));
    end
    
    netcdf.defVarFill(thenc,theID(i),false,fillval);
    netcdf.putAtt(thenc,theID(i), 'axis', 'time'); % define that these are time (T) axis variables

    %%% other things that we could define for each variable but won't for now
    %     valid_range %%% ugh
    %     actual_range %%% ugh     
    %     platform %%% leaving out because it's the same for all variables in file... ? 
    %     statistic  %%% leave out because they are all means? 
    %     grid_mapping %%% only seems necessary if the dimensions are not lat/lon
    %     coordinates %%% only seems necessary if the dimensions are not lat/lon

end

%% define global attributes: Part 1
% % % TGT:   ncei_code = '3250'
% % %        call_sign = 'KTDQ'
% % %        imo_code  = '8814419'

gvarid = netcdf.getConstant('GLOBAL');

acknowledgement_nc = 'NOAA Climate Program Office Climate Variability and Predictability (NOAA CPO CVP) and Office of Naval Research Marine Meteorology (ONR)';
cdm_data_type_nc = 'Trajectory';
comment_nc = 'Corrections and Data Quality Notes not contained in global or variable attributes: Unavailable data, bad data, and data within restricted Exclusive Economic Zones were assigned _FillValue = NaN. Please use the variables named flag_bad_ship and flag_bad_bulk to further mask out questionable or non-ideal data points depending on the application for state variables and bulk fluxes respectively. Variables that are either excluded or are not continuously available that would normally be included: To convert wind speed and direction to u and v components such that zonal u is positive to E, meridional v is positive to N in meteorological from convention: v = spd .* cos((dir+180)*pi/180); u = spd .* sin((dir+180)*pi/180.';
% contributor_name_nc = 'Jay Orson Hyde and Joe (Hadrina) Fernando, Notre Dame University';
% contributor_role_nc = 'installation, maintenance, demobilization, and shipping of instruments';
conventions_nc = 'CF-1.6 ACCD-1.3';
coverage_content_type_nc = 'physicalMeasurement, qualityInformation, modelResult, coordinate';
creator_email_nc = 'elizabeth.thompson@noaa.gov';
creator_institution_nc = '(1) NOAA Physical Sciences Lab (PSL); (2) CIRES Cooperative Institute for Research in Environmental Sciences at the University of Colorado Boulder in partnership with NOAA PSL, (3) Oregon State University';
creator_name_nc = 'Elizabeth J. Thompson (1), Chris Fairall (1), Ludovic Bariteau (2), Sergio Pezoa (1), Byron Blomquist (2), Simon de Szoeke (3)';
creator_type_nc = 'group';
creator_url_nc = 'https://psl.noaa.gov/boundary-layer/';
geospatial_lat_units_nc = 'degrees_north';
geospatial_lon_units_nc = 'degrees_east';
geospatial_vertical_units_nc = 'meters';
% geospatial_vertical_positive = 'up'; % we don?t use this really because both depth and height are defined positive
history_nc = 'v0: original data, v1: first release, v2: reformatted second release';
id_nc = 'doi = not yet assigned';
institution_nc = creator_institution_nc;
% instrument_nc = 'i';
instrument_vocabulary_nc = 'GCMD Version 12.3';
keywords_library_nc = 'GCMD Version 12.3';
licence_nc = 'Please acknowledge data and MISOBOB and PISTON 2018 project according to global attribute info: acknowledgement, creator_name, creator_institution. These data may be redistributed and used without restriction.';
% metadata_link_nc = 'xxx';
naming_authority_nc = 'gov.noaa.ncei';
platform_nc = 'Research Vessel Thomas G. Thompson';
platform_vocabulary_nc = 'GCMD Version 12.3';
processing_level_nc = 'processed and quality controlled';
% product version is defined in program based on another field
program_nc = 'Funding agencies: NOAA Climate Program Office Climate Variability and Predictability, ONR Marine Meteorology ';
project_nc = 'Monsoon Intra-Seasonal Oscillations in the Tropical Indian Ocean and the Bay of Bengal (MISOBOB) in collaboration with Propagation of Intra-Seasonal Tropical Oscillations (PISTON, by ONR) and Years of the Maritime Continent (YMC, by NOAA CPO CVP) ';
% publisher_email_nc = 'elizabeth.thompson@noaa.gov';
% publisher_institution_nc = 'NOAA Physical Sciences Laboratory';
% publisher_name_nc = 'Elizabeth J. Thompson';
% publisher_type_nc = 'person';
% publisher_url_nc = 'https://psl.noaa.gov/boundary-layer/';
references_nc = ['Fairall et al. 1996a JGR https://doi.org/10.1029/95JC03190 ...'  ...
'Fairall et al. 1996b JGR https://doi.org/10.1029/95JC03205 .... ' ...
'Fairall et al. 2003 JClim https://doi.org/10.1175/1520-0442(2003)016%3C0571:BPOASF%3E2.0.CO;2 ... ' ....
'Edson et al. 2013 JPO with corrigendum: the value should be m = 0.0017, and not m = 0.017 as originally appeared https://doi.org/10.1175/JPO-D-12-0173.1'];
source_nc = 'observations from NOAA PSL sensors and ship permanent sensors (_ship named variables, which are less accurate), derivations from those observations using eddy covariance and inertial dissipation methods of estimating fluxes, model results from COARE 3.6 bulk air-sea flux algorithm. Wave parameters were not used as input to COARE since they were either unavailable or not consistently available at all times.';
standard_name_vocabulary_nc = 'CF Standard Name Table, Version 77, 19 January 2021, https://cfconventions.org/Data/cf-standard-names/77/build/cf-standard-name-table.html';
summary_nc = 'Data collected from this cruise is critical for supporting the study of physical oceanography, air-sea interaction, tropical meteorology, as well as global weather and climate variability and predictability. This includes improvement to our fundamental understanding of these processes in the far western Pacific Ocean and their influence around the globe including the Continental United States. The data will also support improvement and validation of prediction models including parameterizations.';
title_nc = 'Ship-based data of navigation, meteorology, seawater, and air-sea fluxes from R/V Thomas G. Thompson in the Philippine Sea from the 2018 PISTON / YMC experiments (Propagation of Intraseasonal Tropical Oscillations by ONR / Years of the Maritime Continent by NOAA)';
sea_name_nc = 'Philippine Sea';
ncei_template_version_nc = 'netCDF_single_trajectory_v2.0 ';


keywords_nc = [...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC PRESSURE > AIR MASS/DENSITY; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC PRESSURE > ATMOSPHERIC PRESSURE MEASUREMENTS; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC PRESSURE > SEA LEVEL PRESSURE; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC PRESSURE > SURFACE PRESSURE; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > HEAT FLUX; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > INCOMING SOLAR RADIATION; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > LONGWAVE RADIATION > DOWNWELLING LONGWAVE RADIATION; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > LONGWAVE RADIATION > UPWELLING LONGWAVE RADIATION; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > LONGWAVE RADIATION; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > NET RADIATION; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > OUTGOING LONGWAVE RADIATION; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > RADIATIVE FLUX; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > SHORTWAVE RADIATION > DOWNWELLING SHORTWAVE RADIATION; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > SHORTWAVE RADIATION; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC RADIATION > SOLAR RADIATION; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC TEMPERATURE > SURFACE TEMPERATURE > AIR TEMPERATURE; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC TEMPERATURE > SURFACE TEMPERATURE > DEW POINT TEMPERATURE; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC TEMPERATURE > SURFACE TEMPERATURE > POTENTIAL TEMPERATURE; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC TEMPERATURE > SURFACE TEMPERATURE > SKIN TEMPERATURE; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC TEMPERATURE > SURFACE TEMPERATURE > VIRTUAL TEMPERATURE; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC WATER VAPOR > WATER VAPOR INDICATORS > HUMIDITY > RELATIVE HUMIDITY; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC WATER VAPOR > WATER VAPOR INDICATORS > HUMIDITY > SATURATION SPECIFIC HUMIDITY; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC WATER VAPOR > WATER VAPOR INDICATORS > HUMIDITY > SPECIFIC HUMIDITY; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC WATER VAPOR > WATER VAPOR INDICATORS > HUMIDITY; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC WATER VAPOR > WATER VAPOR INDICATORS > WATER VAPOR; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC WATER VAPOR > WATER VAPOR PROCESSES > EVAPORATION; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC WINDS > SURFACE WINDS > U/V WIND COMPONENTS; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC WINDS > SURFACE WINDS > WIND DIRECTION; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC WINDS > SURFACE WINDS > WIND SPEED; ' ...
'EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC WINDS > SURFACE WINDS; ' ...
'EARTH SCIENCE > ATMOSPHERE > CLOUDS > CLOUD RADIATIVE TRANSFER > CLOUD RADIATIVE FORCING; ' ...
'EARTH SCIENCE > ATMOSPHERE > PRECIPITATION > LIQUID PRECIPITATION > LIQUID SURFACE PRECIPITATION RATE; ' ...
'EARTH SCIENCE > ATMOSPHERE > PRECIPITATION > LIQUID PRECIPITATION > RAIN; ' ...
'EARTH SCIENCE > ATMOSPHERE > PRECIPITATION > LIQUID PRECIPITATION; ' ...
'EARTH SCIENCE > ATMOSPHERE > PRECIPITATION > PRECIPITATION AMOUNT; ' ...
'EARTH SCIENCE > ATMOSPHERE > PRECIPITATION > PRECIPITATION RATE; ' ...
'EARTH SCIENCE > OCEANS > OCEAN CIRCULATION > OCEAN CURRENTS > SUBSURFACE CURRENTS; ' ...
'EARTH SCIENCE > OCEANS > OCEAN CIRCULATION > OCEAN CURRENTS > SURFACE CURRENTS; ' ...
'EARTH SCIENCE > OCEANS > OCEAN CIRCULATION > OCEAN CURRENTS > SURFACE SPEED; ' ...
'EARTH SCIENCE > OCEANS > OCEAN CIRCULATION > OCEAN CURRENTS; ' ...
'EARTH SCIENCE > OCEANS > OCEAN HEAT BUDGET > HEAT FLUX > LATENT HEAT FLUX; ' ...
'EARTH SCIENCE > OCEANS > OCEAN HEAT BUDGET > HEAT FLUX > SENSIBLE HEAT FLUX; ' ...
'EARTH SCIENCE > OCEANS > OCEAN HEAT BUDGET > HEAT FLUX > TURBULENT HEAT FLUX; ' ...
'EARTH SCIENCE > OCEANS > OCEAN HEAT BUDGET > HEAT FLUX; ' ...
'EARTH SCIENCE > OCEANS > OCEAN HEAT BUDGET > HEATING RATE; ' ...
'EARTH SCIENCE > OCEANS > OCEAN HEAT BUDGET > LONGWAVE RADIATION; ' ...
'EARTH SCIENCE > OCEANS > OCEAN HEAT BUDGET > REFLECTANCE; ' ...
'EARTH SCIENCE > OCEANS > OCEAN HEAT BUDGET > SHORTWAVE RADIATION; ' ...
'EARTH SCIENCE > OCEANS > OCEAN HEAT BUDGET; ' ...
'EARTH SCIENCE > OCEANS > OCEAN PRESSURE > SEA LEVEL PRESSURE; ' ...
'EARTH SCIENCE > OCEANS > OCEAN TEMPERATURE > SEA SURFACE TEMPERATURE > SEA SURFACE  SUBSKIN TEMPERATURE; ' ...
'EARTH SCIENCE > OCEANS > OCEAN TEMPERATURE > SEA SURFACE TEMPERATURE > SEA SURFACE FOUNDATION TEMPERATURE; ' ...
'EARTH SCIENCE > OCEANS > OCEAN TEMPERATURE > SEA SURFACE TEMPERATURE > SEA SURFACE SKIN TEMPERATURE; ' ...
'EARTH SCIENCE > OCEANS > OCEAN TEMPERATURE > SEA SURFACE TEMPERATURE; ' ...
'EARTH SCIENCE > OCEANS > OCEAN TEMPERATURE > WATER TEMPERATURE; ' ...
'EARTH SCIENCE > OCEANS > OCEAN TEMPERATURE; ' ...
'EARTH SCIENCE > OCEANS > OCEAN WINDS > SURFACE WINDS > WIND DIRECTION; ' ...
'EARTH SCIENCE > OCEANS > OCEAN WINDS > SURFACE WINDS > WIND SPEED; ' ...
'EARTH SCIENCE > OCEANS > OCEAN WINDS > SURFACE WINDS; ' ...
'EARTH SCIENCE > OCEANS > OCEAN WINDS; ' ...
'EARTH SCIENCE > OCEANS > SALINITY/DENSITY > CONDUCTIVITY > SURFACE CONDUCTIVITY; ' ...
'EARTH SCIENCE > OCEANS > SALINITY/DENSITY > CONDUCTIVITY; ' ...
'EARTH SCIENCE > OCEANS > SALINITY/DENSITY > DENSITY; ' ...
'EARTH SCIENCE > OCEANS > SALINITY/DENSITY > OCEAN SALINITY > OCEAN SURFACE SALINITY; ' ...
'EARTH SCIENCE > OCEANS > SALINITY/DENSITY > OCEAN SALINITY > PRACTICAL SALINITY; ' ...
'EARTH SCIENCE > OCEANS > SALINITY/DENSITY > OCEAN SALINITY; ' ...
'EARTH SCIENCE > OCEANS > SALINITY/DENSITY > SALINITY; ' ...
'EARTH SCIENCE > OCEANS > SALINITY/DENSITY ' ...
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atts_traj = {...
'acknowledgement';'cdm_data_type';'comment';
'conventions';'coverage_content_type';'creator_email';'creator_institution';...
'creator_name';'creator_type';'creator_url';'geospatial_lon_units';'geospatial_lat_units';...
'geospatial_vertical_units';'history';...
'id';'institution';'instrument_vocabulary';'keywords';'keywords_library';...
'licence';'naming_authority';'platform';'platform_vocabulary';...
'processing_level';'product_version';'program';'project';
'references';'source';'standard_name_vocabulary';'summary';
'title';'sea_name';'ncei_template_version'};

atts = atts_traj;

for h = 1:length(atts)
    eval(['netcdf.putAtt(thenc,gvarid, "' atts{h} '", ' atts{h} '_nc);']);
end

%% define global attributes: Part 1... these are easier to do with commands

netcdf.putAtt(thenc,gvarid,'date_created',datestr(now));
netcdf.putAtt(thenc,gvarid,'date_issued',datestr(now));
netcdf.putAtt(thenc,gvarid,'date_metadata_modified',datestr(now));
netcdf.putAtt(thenc,gvarid,'date_modified',datestr(now));
netcdf.putAtt(thenc,gvarid,'product_version',product_version_nc);

geospatial_lat_min_nc = sprintf('%8.3f', min(a.lat));
geospatial_lat_max_nc = sprintf('%8.3f', max(a.lat));
geospatial_lon_min_nc = sprintf('%8.3f', min(a.lon));
geospatial_lon_max_nc = sprintf('%8.3f', max(a.lon));
geospatial_bounds_nc = ['POLYGON [' geospatial_lon_min_nc ', ' geospatial_lon_max_nc ', ' geospatial_lat_min_nc ', ' geospatial_lat_max_nc ']'];
netcdf.putAtt(thenc,gvarid,'geospatial_lat_min',geospatial_lat_min_nc);
netcdf.putAtt(thenc,gvarid,'geospatial_lat_max',geospatial_lat_max_nc);
netcdf.putAtt(thenc,gvarid,'geospatial_lon_min',geospatial_lon_min_nc);
netcdf.putAtt(thenc,gvarid,'geospatial_lon_max',geospatial_lon_max_nc);
netcdf.putAtt(thenc,gvarid,'geospatial_lat_bounds',geospatial_bounds_nc);
netcdf.putAtt(thenc,gvarid,'geospatial_vertical_min',sprintf('%5.2f',-1*zsea_ship));
netcdf.putAtt(thenc,gvarid,'geospatial_vertical_max',sprintf('%5.2f',zu));

netcdf.putAtt(thenc,gvarid,'time_coverage_start', [datestr(min(a.t), 0) ' UTC']);
netcdf.putAtt(thenc,gvarid,'time_coverage_end',   [datestr(max(a.t), 0) ' UTC']);
netcdf.putAtt(thenc,gvarid,'time_coverage_duration', [sprintf('%6.3f',ndays) ' days']);

% % tell the netcdf file that you are done defining how it will work
netcdf.endDef(thenc);

%% write the data to the variables you just prepared.

for i = 1:ndt
   eval(['data = thedata.' thefields{i} ';']);
   netcdf.putVar(thenc,theID(i),data);
end


%% end process. IT IS VERY IMPORTANT TO CLOSE THE FILE!
netcdf.close(thenc)

disp(['saved netcdf file: ' thencname]);
% ncdisp([thedir thencname]);

end

%% check files;

% ncdisp([thedir thencname]);

plot_map = 0;
if plot_map == 1
    map_all_cruise(c.lon,c.lat,ptitle,PosLims)
    annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA PSL'},'FontSize',10,'FitBoxToText','off','LineStyle','none');
    print(graphdevice,[plotdir cruise '_TrackMap'  graphformat]);
end

plot_check_all_times = 1;
if plot_check_all_times == 1
    figure;
    counter = 1;
    thefields = fields(c);
    nc = length(thefields);
    for i = 1:4:nc
        clf;
        for j = 1:4
           if (i+j-1) <= nc
               subplot(2,2,j); hold on;
               plot(a.t, a.(thefields{i+j-1}),'x');
               plot(e60.t, e60.(thefields{i+j-1}),'.');
               if j == 1
                   legend('10 min','1 hr','location',[0.46,0.464, 0.05, 0.05],'box','on');
               end
               grid on;
               var_name = {strrep(thefields{i+j-1},'_',' ')};
               title(var_name);
               xlim([min(a.t) max(a.t)]);
               datetick('x','mm/dd','keeplimits');
               grid on;
           end
        end
       print(graphdevice,[plotdir cruise '_check_' sprintf('%i',counter)  graphformat]);
       counter = counter + 1;
    end  % for all vars
end %%% if making plots
