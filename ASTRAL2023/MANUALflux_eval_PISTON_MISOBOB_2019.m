%{
 PISTON / MISOBOB :: 2019-08 :: BWB
 Process daily met files, make diagnostic plots and save intermediate
 1-min and 10-min data files.
 *** PSD heading not working on MISOBOB ***
 Calls:
  read_met_day:     PSD met data from loggers 1 & 2 @ 1 min
  read_gps_day:     PSD gps course & speed @ 1 Hz
  read_hed_day:     PSD gps heading & pitch @ 1 Hz
  read_motion_day:  PSD MotionPak @ 1 Hz
  read_sonic_day:   PSD sonic wind @ 1 Hz
  read_licor_day:   PSD licor fast water vapor / CO2 @ 1 Hz
  read_scs_day:     Sally Ride ship data @ 1 Hz
  read_wxt_day      PSD WXT met data from ship's A-frame (aft) @ 1 Hz
%}

%%  initialize run parameters
clear;
close('all');
fclose('all');
warning ('off','MATLAB:MKDIR:DirectoryExists');

% set leg number manually, 1-3, as defined below
% leg 1 jd 149:174, MISOBOB1 Phuket to Phuket
% leg 2 jd 186:214, MISOBOB2 Chennai to Chennai
% leg 3 jd 240:257, PISTON1 Taiwan to Taiwan
LegNo = 3;
LegStr = sprintf('%1i',LegNo);

% set day-of-year range for this run
jdStart = 246; jdStop = 246;

plotit = true; % for doing plots.
prtit = true;  % for saving plots--plotit must also be true
graphformat = '.png';  % select graphics format files
graphdevice = '-dpng'; % select graphic device
cruise = 'PISTON_MISOBOB_2019';     % Cruise name
yyyy = '2019';                      % year string
ship = 'Sally_Ride';                % Research vessel

% system specific path defs
sysType = computer;
username=char(java.lang.System.getProperty('user.name'));
if strncmp(sysType,'MACI64',7)     % set Mac paths
	data_drive = '/Users/byronb/Raw_Data_Files.noindex';
    path_prog = fullfile(data_drive,cruise,ship,'Scientific_analysis','programs');
elseif strncmp(sysType,'PCWIN',7)  % set PSD DAC paths
    data_drive = 'D:\Data';
    path_prog = fullfile(data_drive,cruise,ship,'Scientific_analysis','programs');
end

% matlab script path
restoredefaultpath
cd(fullfile(path_prog,'flux'));
addpath(genpath(fullfile(path_prog,'flux')));
rehash toolboxcache;

% paths relative to data_drive - these should be preexisting
path_python = fullfile(path_prog,'python');
path_raw_data = fullfile(data_drive,cruise,ship,'flux','Raw');
% define folders for saving data and plots
path_proc_data = fullfile(data_drive,cruise,ship,'flux','Field_Processed');
mkdir(path_proc_data);
mkdir(fullfile(path_proc_data,['flux_',LegStr]));
mkdir(fullfile(path_proc_data,['gps_',LegStr]));
mkdir(fullfile(path_proc_data,['met_',LegStr]));
mkdir(fullfile(path_proc_data,['rad_',LegStr]));
path_raw_images = fullfile(data_drive,cruise,ship,'flux','Raw_Images');
mkdir(path_raw_images);
mkdir(fullfile(path_raw_images,['qSpectra_',LegStr]));
mkdir(fullfile(path_raw_images,['Rwspd_',LegStr]));
mkdir(fullfile(path_raw_images,['sonicQC_',LegStr]));
mkdir(fullfile(path_raw_images,['tSpectra_',LegStr]));
mkdir(fullfile(path_raw_images,['uSpectra_',LegStr]));
mkdir(fullfile(path_raw_images,['wSpectra_',LegStr]));
mkdir(fullfile(path_raw_images,['SL_pressure_',LegStr]));
mkdir(fullfile(path_raw_images,['Temps_',LegStr]));
mkdir(fullfile(path_raw_images,['RH_',LegStr]));
mkdir(fullfile(path_raw_images,['Rainrate_',LegStr]));
mkdir(fullfile(path_raw_images,['T_RH_fan_',LegStr]));
mkdir(fullfile(path_raw_images,['SST_',LegStr]));
mkdir(fullfile(path_raw_images,['Track_plot_',LegStr]));
mkdir(fullfile(path_raw_images,['COG_SOG_',LegStr]));
mkdir(fullfile(path_raw_images,['Heading_',LegStr]));
mkdir(fullfile(path_raw_images,['pitch_',LegStr]));
mkdir(fullfile(path_raw_images,['Relative_Wind_',LegStr]));
mkdir(fullfile(path_raw_images,['True_Wind_',LegStr]));
mkdir(fullfile(path_raw_images,['Motion_',LegStr]));
mkdir(fullfile(path_raw_images,['Licor_AGC_',LegStr]));
mkdir(fullfile(path_raw_images,['Net_Heat_',LegStr]));
mkdir(fullfile(path_raw_images,['Heat_Flux_Components_',LegStr]));
mkdir(fullfile(path_raw_images,['u-star_',LegStr]));
mkdir(fullfile(path_raw_images,['IR_flux_',LegStr]));
mkdir(fullfile(path_raw_images,['Solar_flux_',LegStr]));

%% other initialization settings

% select lat/lon limits and adjustments based on leg number
switch LegNo
    case 1
        Lonmin = 79; Lonmax = 92; Latmin = 5; Latmax = 17;
        td1_adj = 0; % adjustment for MISOBOB
        tc1_adj = 0; % adjustment for MISOBOB
        td2_adj = 0; % adjustment for MISOBOB
        tc2_adj = 0; % adjustment for MISOBOB
        td_scs_adj = 0;
        tc_scs_adj = 0;
    case 2
        Lonmin = 79; Lonmax = 92; Latmin = 4; Latmax = 19;
        td1_adj = 0; % adjustment for MISOBOB
        tc1_adj = 0; % adjustment for MISOBOB
        td2_adj = 0; % adjustment for MISOBOB
        tc2_adj = 0; % adjustment for MISOBOB
        td_scs_adj = 0;
        tc_scs_adj = 0;
    case 3
        Lonmin = 119; Lonmax = 139; Latmin = 4; Latmax = 27;
        td1_adj = 0; % adjustment for PISTON
        tc1_adj = 0; % adjustment for PISTON
        td2_adj = 0; % adjustment for PISTON
        tc2_adj = 0; % adjustment for PISTON
        td_scs_adj = 0;
        tc_scs_adj = 0;
end
% sea snake adjustment
tsea_adj = 0;
% air temp adjustments
ta_adj = 0;

PosLims = [Lonmin, Lonmax, Latmin, Latmax];
adj = [tsea_adj,ta_adj,td1_adj,tc1_adj,td2_adj,tc2_adj];
SCS_adj = [td_scs_adj,tc_scs_adj];

% sonic/licor options
fsonic = 10;
flicor = 10;
rotationsonic = false;
sonicmodel = 'WindMasterPro';

% Clear sky model configuration
k1 = 0.1;   % aerosol optical depth, band 1
k2 = 0.1;   % aerosol optical depth, band 2
oz = 0.2;   % column ozone
iv = 3.5;   % column water vapor (cm), if not calculated from obs.

% sensor heights above SL, PSD heights estimated w/respect to Sally Ride
% bow mast sensor height given in SAMOS metadata (15.24m).
zu = 14.75;        % wind speed measurement height (m)
zt = 13.75;        % air T measurement height (m)
zq = 13.75;        % air q measurement height (m)
zp = 10.39;            % psd pressure sensor height (m), guess only...
zwxt = 11.31;          % wxt height, guess only...
zuim = 15.24;      % ship sensor heights (m)
ztim = 15.24;
zqim = 15.24;
zpim = 15.24;
dtsg = 3.5;       % depth of SBE 21 intake SST

rdcon = pi/180;

%% loop through days
for ddd = jdStart:jdStop
    %% run python scripts on gps and WXT raw files
    [m,d] = yd2md(str2num(yyyy), ddd);
    Vdate = [str2num(yyyy),m,d];
    path_working_ddd = fullfile(path_raw_data,[yyyy(3:4),sprintf('%03i',ddd)]);
    % convert PSD gps and met3 files to gprm and wxt format with python script
    % python language must be installed on the computer.
    cd(path_working_ddd);

    files = dir('gps*.txt');
    if ~isempty(files)
        system(['python ',fullfile(path_python,'parseGpsFiles.py')]);
    end

%     files = dir('hed0*.txt');
%     if ~isempty(files)
%         system(['python ',fullfile(path_python,'parseHedFiles.py')]);
%     end
    
    files = dir('met3*.txt');
    if ~isempty(files)
        system(['python ',fullfile(path_python,'parseWXTFiles.py')]);
    end

    %% read data - this takes time...
    % read PSD 1-min avg bulk met data - metm is 1440x17 array
    metm = read_met_day_PISTON_MISOBOB_2019(path_working_ddd,ddd,yyyy,adj,zp);

    % read PSD gps - gprm is 1Hz 86400x7 array
    gprm = read_gps_day_PISTON_MISOBOB_2019(path_working_ddd,ddd,yyyy,PosLims);

    % read PSD heading-pitch - hedm is 1Hz 86400x3 array - returns NaN
    % array for MISOBOB - heading not working
    hedm = read_hed_day_PISTON_MISOBOB_2019(path_working_ddd,ddd,yyyy);

    % read PSD motionPak - motm is 1Hz 86400x13
    motm = read_motion_day_PISTON_MISOBOB_2019(path_working_ddd,ddd,yyyy);

    % read PSD sonic - sonm is 1Hz 86400x5
    sonm = read_sonic_day_PISTON_MISOBOB_2019(path_working_ddd,ddd,yyyy,sonicmodel,...
         rotationsonic,fsonic,prtit,path_raw_images,graphformat,graphdevice,LegStr);

    % read PSD Licor - licm is 1Hz 86400x8
    licm = read_licor_day_PISTON_MISOBOB_2019(path_working_ddd,ddd,yyyy,...
                flicor,prtit,path_raw_images,graphformat,graphdevice,LegStr);

    % read ship SCS data - scsm is 1Hz 86400x26
    scsm = read_scs_day_PISTON_MISOBOB_2019(path_working_ddd,ddd,yyyy,PosLims,SCS_adj,zpim);

    % read wxt data from A-frame - wxtm is 1440 x 9 array
    wxtm = read_wxt_day_PISTON_MISOBOB_2019(path_working_ddd,ddd,yyyy,zwxt);

    close all

    %% 1-min averages for some variables

    jd_min = metm(:,1); % ref 1-min timestamp
    delta = double(1.0/1440);
    last = ddd + 1440*delta;
    jd_1bin = (ddd:delta:last);	jd_1bin = jd_1bin'; % ref 1-min bin edges

    % 1-min psd gps data
    gprm1 = interval_avg(gprm(:,1), gprm(:,2:7), jd_1bin); % average everything
    Ngps_1 = gprm1(:,4); % mean ship speed toward N
    Egps_1 = gprm1(:,5); % mean ship speed toward W
    sog_1 = sqrt(Ngps_1.^2 + Egps_1.^2);  % recompute sog from speed components
    cog_1 = atan2(Egps_1, Ngps_1)/rdcon; % recompute cog from speed components
    cog_1 = mod(cog_1 + 360, 360);
    gprm1(:,2) = cog_1; % replace with updated 1-min cog/sog
    gprm1(:,3) = sog_1;

    % 1-min psd heading
    [hed_n, hed_e] = pol2cart(hedm(:,2)*rdcon,1);  % N/E heading components
    hedm = [hedm, hed_n, hed_e];                   % add to array
    hedm1 = interval_avg(hedm(:,1), hedm(:,2:5), jd_1bin); % average everything
    hed_1 = atan2(hedm1(:,5), hedm1(:,4))/rdcon; % recompute hed from avg N/E components
    hed_1 = mod(hed_1 + 360, 360);
    hedm1(:,2) = hed_1; % replace with updated 1-min hed

    [hed_n, hed_e] = pol2cart(scsm(:,12)*rdcon,1);  % N/E heading components
    scsm = [scsm, hed_n, hed_e]; % add heading components to scsm array, cols 28,29
    scsm1 = interval_avg(scsm(:,1), scsm(:,2:29), jd_1bin); % average all scs
    Ngps_1 = scsm1(:,25); % mean ship speed toward N
    Egps_1 = scsm1(:,26); % mean ship speed toward W
    sog_1 = sqrt(Ngps_1.^2 + Egps_1.^2);  % recompute sog from speed components
    cog_1 = atan2(Egps_1, Ngps_1)/rdcon;  % recompute cog from speed components
    cog_1 = mod(cog_1 + 360, 360);
    scsm1(:,4) = cog_1; % replace with updated 1-min cog/sog
    scsm1(:,5) = sog_1;
    hed_1 = atan2(scsm1(:,29), scsm1(:,28))/rdcon; % recompute hed from avg N/E components
    hed_1 = mod(hed_1 + 360, 360);
    scsm1(:,12) = hed_1; % replace with updated 1-min hed

    licm1 =  interval_avg(licm(:,1), licm(:,2:8), jd_1bin); % average all licm

    %% correct ORG offset
    org_bkgd = nanmedian(despike2(metm(:,17)));
    org_adj = 0.06484 - org_bkgd;
    org_recalc = (metm(:,17)+org_adj).^1.87*25 - 0.15;  % recompute rain rate after adjustment

    %% compute 10-min averages
    delta = double(10.0/1440);
    last = ddd + 144*delta;
    jd_10bin = (ddd:delta:last);	jd_10bin = jd_10bin'; % ref 10-min bin edges

    % 10-min SCS
    %{
    1    Decimal DOY
    2    Lat
    3    Lon
    4    COG
    5    SOG, m/s
    6    tsgm  TSG temperature
    7    tssm  TSG salinity
    8    imrh RH, @ 15.5 m
    9    imta air temp, C, @ 15.5 m
    10   imsol downwelling solar radiation, W/m2
    11   orgm rain rate (no data)
    12   lrg heading, deg
    13   imum true wind speed, m/s
    14   imdm true wind direction
    15   imir downwelling IR, W/m2
    16   imPress sealevel pressure, mb
    17   imqa specific humidity, g/kg
    18   imts SST at intake, C
    19   imrosr radiometric SST, C
    20   impir raw PIR thermopile, W/m2
    21   imtc PIR case temp, C
    22   imtd PIR dome temp, C
    23   imrwspd relative wind speed, m/s
    24   imrwdir relative wind dir, +/- 180 deg from bow
    25   Ngps, N ship speed component, m/s
    26   Egps, E ship speed component, m/s
    27   imrl2 recomputed longwave radiation, W/m2
    28   hed_n heading N component
    29   hed_e heading E component
    %}

    scsm10 = interval_avg(scsm(:,1), scsm(:,2:29), jd_10bin); % average everything
    scsm10_std = interval_std(scsm(:,1), scsm(:,2:29), jd_10bin); % std dev of everything
    jd_10min = scsm10(:,1); % reference10-min timestamp
    latm_10 = scsm10(:,2);
    lonm_10 = scsm10(:,3);
    lrg_10n = scsm10(:,28);
    lrg_10e = scsm10(:,29);
    lrgm_10 = atan2(lrg_10e, lrg_10n)/rdcon; % from avg components
    lrgm_10 = mod(lrgm_10 + 360, 360);
    SOG_10n = scsm10(:,25);
    SOG_10e = scsm10(:,26);
    cogm_10 = atan2(SOG_10e, SOG_10n)/rdcon; % from avg components
    cogm_10 = mod(cogm_10 + 360, 360);
    sogm_10 = sqrt(SOG_10e.^2 + SOG_10n.^2);  % from avg components
    sogm_10_std = scsm10_std(:,5);    % std dev in ship speed
    imrh_10 = scsm10(:,8);
    imta_10 = scsm10(:,9);
    imum_10 = scsm10(:,13);
    imrspd_10 = scsm10(:,23);
    % recompute wdir / rwdir
    imdm_10 = interval_avg(scsm(:,1), unwrap(scsm(:,14)*rdcon), jd_10bin);
    imdm_10 = imdm_10(:,2);
    imdm_10 = mod(imdm_10/rdcon + 360*6, 360);
    imrdir_10 = interval_avg(scsm(:,1), unwrap(scsm(:,24)*rdcon), jd_10bin);
    imrdir_10 = imrdir_10(:,2);
    imrdir_10 = mod(imrdir_10/rdcon + 360*6, 360);
    imrdir_10(imrdir_10>180) = imrdir_10(imrdir_10>180) - 360; % +/- 180 deg format

    impress_10 = scsm10(:,16); % @ sea level
    tsgm_10 = scsm10(:,6);
    tssm_10 = scsm10(:,7);
    imts_10 = scsm10(:,18);
    imrosr_10 = scsm10(:,19);
    imsol_10 = scsm10(:,10);
    imir_10 = scsm10(:,15);
    imtc_10 = scsm10(:,21);
    imtd_10 = scsm10(:,22);
    imir2_10 = scsm10(:,27);
    imqa_10 = qair_p([imta_10, imrh_10,impress_10]);

    % 10-min PSD gps
    gprm10 = interval_avg(gprm(:,1), gprm(:,2:7), jd_10bin); % average everything
    Ngps_10 = gprm10(:,4); % ship speed toward N
    Egps_10 = gprm10(:,5); % ship speed toward W
    lat_10 = gprm10(:,6);
    lon_10 = gprm10(:,7);
    sog_10 = sqrt(Ngps_10.^2 + Egps_10.^2);  % sog from speed components
    cog_10 = atan2(Egps_10, Ngps_10)/rdcon; % cog from speed components
    cog_10 = mod(cog_10 + 360, 360);

    % 10-min PSD heading
    hedm10 = interval_avg(hedm(:,1), hedm(:,2:5), jd_10bin); % average everything
    hed_10n = hedm10(:,4);
    hed_10e = hedm10(:,5);
    hed_10 = atan2(hed_10e, hed_10n)/rdcon; % recompute hed from avg N/E components
    hed_10 = mod(hed_10 + 360, 360);
    hedm10(:,2) = hed_10; % replace with updated 10-min hed
%     hed_10 = lrgm_10;
%     hedm10 = lrgm_10;  % use ship heading

    % 10-min PSD met & radiation
    %{
    1   jd_met          1 min timestamp
    2   Tvais           T vaisala, C
    3   Rhvais          RH vaisala
    4   Tsea            sea snake surface SST
    5   psp1            psp1 radiation, W/m2
    6   pir1            pir1 radiation, W/m2
    7   psp2            psp2 radiation, W/m2
    8   pir2            pir2 radiation, W/m2
    9   Tc1             pir1 case temperature
    10  Td1             pir1 dome temperature
    11  Tc2             pir2 case temperature
    12  Td2             pir2 dome temperature
    13  org             rain, mm/hr
    14  press           air pressure, mb
    15  aspir_on        T/RH diagnostic
    16  org_carrier     org diagnostic
    17  org_V           ORG raw measurement, Volts
    %}

    metm10 = interval_avg(metm(:,1), metm(:,2:17), jd_10bin);
    Tvais_10 = metm10(:,2);
    Rhvais_10 = metm10(:,3);
    Tsea_10 = metm10(:,4);
    psp1_10 = metm10(:,5);
    pir1_10 = metm10(:,6);
    psp2_10 = metm10(:,7);
    pir2_10 = metm10(:,8);
    Tc1_10 = metm10(:,9);
    Td1_10 = metm10(:,10);
    Tc2_10 = metm10(:,11);
    Td2_10 = metm10(:,12);
    org_10 = metm10(:,13);
    if all(isnan(org_10)); org_10 = 0; end
    press_10 = metm10(:,14);  % sea level P
    aspir_on_10 = metm10(:,15);
    org_carrier_10 = metm10(:,16);
    qair_10 = qair_p([Tvais_10,Rhvais_10,press_10]);
    qsea_10 = qsea_p([Tsea_10,press_10]);
    org_10_recalc = interval_avg(metm(:,1), org_recalc, jd_10bin);
    org_10_recalc = org_10_recalc(:,2);

    % compute 10-min PSD winds from sonic
    sonm10 = interval_avg(sonm(:,1), sonm(:,2:end), jd_10bin); % average everything
    % compute wind spd/dir - apply flow distortion corrections for mean
    % spd, using ship heading for corrections
    [rwspd_10,rwdir_10] = uv_to_sd(sonm10(:,2)/0.95,sonm10(:,3)/1.15);
    rwdir_10(rwdir_10>180) = rwdir_10(rwdir_10>180)-360;    % +/- 180 deg format
    [wspd_10,wdir_10,~,~] = uv_rel_to_sd_true(sonm10(:,2)/0.95,sonm10(:,3)/1.15,lrgm_10,cogm_10,sogm_10);

    % 10-min psd licor
    licm10 = interval_avg(licm(:,1), licm(:,2:end), jd_10bin);

    % 10-min wxt met
    wxt10 = interval_avg(wxtm(:,1), wxtm(:,2:9), jd_10bin);      % average everything
    [rwspd_10_wxt,rwdir_10_wxt] = uv_to_sd(wxt10(:,8),wxt10(:,9)); % recompute wspd/wdir
    rr = rwdir_10_wxt>180; rwdir_10_wxt(rr) = rwdir_10_wxt(rr)-360;
    [wspd_10_wxt,wdir_10_wxt,~,~] = uv_rel_to_sd_true(wxt10(:,8),wxt10(:,9),lrgm_10,cogm_10,sogm_10);
    Twxt_10 = wxt10(:,4);
    RHwxt_10 = wxt10(:,5);
    Pwxt_10 = wxt10(:,6);
    rainwxt_10 = wxt10(:,7);  % mm/hr

    %% Some plots...

    if plotit
        % PRESSURE
        figure; plot(jd_10min,press_10,'r-',jd_10min,impress_10,'g-',jd_10min,Pwxt_10,'k-','linewidth',1); grid;
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  Sea Level Pressures',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none')
        xlabel('Hour (UTC)'); ylabel('SL Pressure (mb)');
        axis([ddd ddd+1 median(press_10(~isnan(press_10)))-5 median(press_10(~isnan(press_10)))+5]);
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        legend('PSD','SCS','WXT','location','BestOutside');
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['SL_pressure_',LegStr],['Pressure_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % TEMPERATURES
        imta_10(imta_10<15) = NaN;  % remove bad ship temperatures
%         figure; plot(jd_10min,Tc1_10,'r',jd_10min,Td1_10,'y',jd_10min,Tc2_10,'m',jd_10min,Td2_10,'k',jd_10min,Tair_wxt_10,'b',jd_10min,Tvais_10,'g.',jd_10min,imta_10,'k.','linewidth',1);
        figure; plot(jd_10min,Tc1_10+tc1_adj,'r',jd_10min,Td1_10+td1_adj,'y',jd_10min,Tc2_10+tc2_adj,...
            'm',jd_10min,Td2_10+td2_adj,'k',jd_10min,Tvais_10,'g.',jd_10min,Twxt_10,'r.',...
            jd_10min,imta_10,'k.',jd_10min,imtc_10+tc_scs_adj,'r+',jd_10min,imtd_10+td_scs_adj,'b+','linewidth',1);
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  PSD temperatures',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        legend('Tc1','Td1','Tc2','Td2','Tvais','WXT','SCS','SCS-Tc','SCS-Td','Location','BestOutside');
        xlabel('Hour (UTC)'); ylabel('Temperature (degC)');
        xlim([ddd ddd+1]); grid;
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['Temps_',LegStr],['Temperatures_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
            ppath = fullfile(path_raw_images,['Temps_',LegStr],['Temperatures_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd)]);
            savefig(ppath);
        end

        % check PIR temp offsets from night air temps
        zz = jd_10min>ddd+12/24 & jd_10min<ddd+22/24;
        tc1_adj2 = nanmean1(Tvais_10(zz)-Tc1_10(zz)); disp(['tc1_adj = ',sprintf('%4.2f',tc1_adj2)]);
        td1_adj2 = nanmean1(Tvais_10(zz)-Td1_10(zz)); disp(['td1_adj = ',sprintf('%4.2f',td1_adj2)]);
        tc2_adj2 = nanmean1(Tvais_10(zz)-Tc2_10(zz)); disp(['tc2_adj = ',sprintf('%4.2f',tc2_adj2)]);
        td2_adj2 = nanmean1(Tvais_10(zz)-Td2_10(zz)); disp(['td2_adj = ',sprintf('%4.2f',td2_adj2)]);
        figure; plot(jd_10min(zz),Tvais_10(zz)-Tc1_10(zz),'r',jd_10min(zz),Tvais_10(zz)-Td1_10(zz),'y',...
            jd_10min(zz),Tvais_10(zz)-Tc2_10(zz),'m',jd_10min(zz),Tvais_10(zz)-Td2_10(zz),'k','linewidth',1);
        title(sprintf('%s %4.2f %s %4.2f %s %4.2f %s %4.2f','PIR Delta Temps: tc1 ',tc1_adj2,' ; td1 ',td1_adj2,...
            ' ; tc2 ',tc2_adj2,' ; td2 ',td2_adj2),'FontWeight','Bold','Interpreter','none');
        legend('Tc1','Td1','Tc2','Td2','Location','BestOutside');
        xlabel('Hour (UTC)'); ylabel('Temperature (degC)'); grid;
        set(gca(gcf),'XTick',ddd+12/24:2/24:ddd+22/24); datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='22:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['Temps_',LegStr],['TcTd_delta T_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % RELATIVE HUMIDITY
        figure; plot(jd_10min, Rhvais_10, 'r',jd_10min, imrh_10,'g',jd_10min, RHwxt_10,'k','linewidth',1); grid;
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  Relative Humidity',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        xlabel('Hour (UTC)'); ylabel('Relative Humidity (%)');
        xlim([ddd ddd+1]); legend('PSD','SCS','WXT','location','BestOutside');
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['RH_',LegStr],['RH_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % ORG background adjustment
        figure; plot(metm(:,1),metm(:,17),'r.',metm(:,1),despike2(metm(:,17)),'g.',...
            [metm(1,1),metm(end,1)],[org_bkgd,org_bkgd],'k--',[metm(1,1),metm(end,1)],[0.06484,0.06484],'k:','linewidth',2); grid;
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  ORG offset voltage ',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        xlabel('Hour (UTC)'); ylabel('Volts'); xlim([ddd ddd+1]); ylim([0,0.1]);
        text(metm(30,1),0.025,['org-adj = ',sprintf('%6.4f',org_adj),' V'],'fontsize',12);
        set(gca(gcf),'XTick',ddd:2/24:ddd+1);datetick('x','HH:MM','keepticks');
        legend({'raw','despike','median offset','target offset'},'location','south','fontsize',12);
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['Rainrate_',LegStr],['ORG_offset_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % RAIN RATE
        figure; plot(jd_10min, org_10, 'r',jd_10min, org_10_recalc, 'b',jd_10min, rainwxt_10, 'k','linewidth',1); grid;
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  Rain Rate ',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        xlabel('Hour (UTC)'); ylabel('Rain Rate (mm/hr)'); xlim([ddd ddd+1]);
        set(gca(gcf),'XTick',ddd:2/24:ddd+1);datetick('x','HH:MM','keepticks');
        legend({'ORG','ORG adj','WXT'},'location','BestOutside');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['Rainrate_',LegStr],['Rainrate_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % PSD T/RH ASPIRATOR FAN
        figure; plot(jd_10min, aspir_on_10,'b','linewidth',1);
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  PSD T/RH aspirator fan',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        xlabel('Hour (UTC)');ylabel('Aspirator Fan Speed (rpm)'); xlim([ddd,ddd+1]);
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks'); grid;
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['T_RH_fan_',LegStr],['Backflow_indicator_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % SST
        Tsea_10(Tsea_10<0) = NaN; tsgm_10(tsgm_10<0) = NaN;
        figure; plot(jd_10min, Tsea_10,'r',jd_10min, tsgm_10,'b','linewidth',1); grid;
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  SST',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        xlabel('Hour (UTC)');ylabel('Sea Surface Temperature (degC)'); xlim([ddd ddd+1]);
        legend('tsnk','tsgm','location','BestOutside');
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['SST_',LegStr],['SST_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % LAT/LON TRACK PLOT - use ship lat/lon
        map_PISTON_MISOBOB_2019(lonm_10,latm_10,cruise,PosLims,ddd)
        if prtit
            ppath = fullfile(path_raw_images,['Track_plot_',LegStr],['GPS_track_LAT_LON_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % GPS SOG/COG
        figure;
        subplot(2,1,1); plot(jd_10min, sog_10, jd_10min, sogm_10,'linewidth',1);
        xlabel('Hour (UTC)');ylabel('SOG (m/s)'); legend('PSD','SCS','location','BestOutside'); grid;
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  GPS SOG/COG',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        xlim([ddd ddd+1]); ylim([0 8]); set(gca(gcf),'XTick',ddd:2/24:ddd+1);datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        subplot(2,1,2); plot(jd_10min, cog_10,'b.', jd_10min, cogm_10,'r.'); legend('PSD','SCS','location','BestOutside');
        xlabel('Hour (UTC)'); ylabel('COG (deg)'); axis([ddd ddd+1 0 360]); grid;
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        orient tall;
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['COG_SOG_',LegStr],['COG_SOG_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % HEADING
        figure; plot(jd_10min, hed_10,'b.', jd_10min, lrgm_10, 'r.');
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  Heading',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        xlabel('Hour (UTC)');ylabel('Heading (deg)'); legend('PSD','SCS','location','BestOutside');
        axis([ ddd ddd+1 0 360]); grid;
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['Heading_',LegStr],['Heading_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % pitch ANGLE
        pitch = hedm(:,3); % this is at 1Hz
        ii = find(isnan(pitch)); pitch(ii) = nanmean1(pitch);
        figure; plot(hedm(:,1),pitch,'b-'); grid;
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  Pitch Angle ',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        xlabel('Hour (UTC)'); ylabel('Pitch (deg)'); xlim([ddd ddd+1]);
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['Pitch_',LegStr],['Pitch_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % RELATIVE WIND SPEED & DIRECTION
        figure;
        subplot(2,1,1);
        plot(jd_10min, rwdir_10,'r.', jd_10min, imrdir_10, 'g.', jd_10min, rwdir_10_wxt, 'k.',...	% data points
            [ddd ddd+1],[60 60],'k:', ...
            [ddd ddd+1],[-60 -60],'k:', ...
            [ddd ddd+1],[90 90],'k-', ...
            [ddd ddd+1],[-90 -90],'k-','linewidth',1);
        xlabel('Hour (UTC)'); ylabel('Relative Wind Direction (deg)'); axis([ddd ddd+1 -180 180]); grid;
        legend('PSD','SCS','WXT','location','BestOutside');
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  Relative Wind Direction',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        subplot(2,1,2); plot(jd_10min, rwspd_10, 'r',jd_10min, imrspd_10, 'g',jd_10min, rwspd_10_wxt, 'k','linewidth',1);
        legend('PSD','SCS','WXT','location','BestOutside'); grid;
        xlabel('Hour (UTC)'); ylabel('Relative Wind Speed (m/s)'); xlim([ddd ddd+1]);
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  Relative Wind Speed',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        orient tall;
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['Relative_Wind_',LegStr],['Relative_Winds_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % TRUE WIND SPEED & DIRECTION wspd_10_wxt
        figure;
        subplot(2,1,1);
        plot(jd_10min, wdir_10,'r.', jd_10min, imdm_10, 'g.', jd_10min, wdir_10_wxt, 'k.','linewidth',1); grid;
        xlabel('Hour (UTC)'); ylabel('Wind Direction (deg)'); axis([ddd ddd+1 -5 365]);
        legend('PSD','SCS','WXT','location','Best');
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  True Wind Direction',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
        set(gca(gcf),'YTick',0:45:360);
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        subplot(2,1,2); plot(jd_10min, wspd_10, 'r', jd_10min, imum_10, 'g', jd_10min, wspd_10_wxt, 'k','linewidth',1); grid;
        xlabel('Hour (UTC)'); ylabel('Wind Speed (m/s)'); xlim([ddd ddd+1]);
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  True Wind Speed',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        orient tall;
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['True_Wind_',LegStr],['True_Winds_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % MOTION PLOT - 1 Hz
        figure;
        subplot(3,2,1);plot(motm(:,1),motm(:,2),'.','markersize',2);ylabel('accx (m.s^-^2)','FontSize',8);
        ylim([-5 5]);set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        grid;
        subplot(3,2,2);plot(motm(:,1),motm(:,3),'.','markersize',2);ylabel('accy (m.s^-^2)','FontSize',8);
        ylim([-5 5]);set(gca(gcf),'XTick',ddd:2/24:ddd+1);datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        grid;
        subplot(3,2,3);plot(motm(:,1),motm(:,4),'.','markersize',2);ylabel('accz (m.s^-^2)','FontSize',8);
        ylim([5 15]);set(gca(gcf),'XTick',ddd:2/24:ddd+1);datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        grid;
        subplot(3,2,4);plot(motm(:,1),motm(:,5),'.','markersize',2);ylabel('ratex (rad/s)','FontSize',8);
        ylim([-5 5]*.02);set(gca(gcf),'XTick',ddd:2/24:ddd+1);datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        grid;
        subplot(3,2,5);plot(motm(:,1),motm(:,6),'.','markersize',2);ylabel('ratey (rad/s)','FontSize',8);xlabel('Hour (UTC)','FontSize',8);
        ylim([-5 5]*.02);set(gca(gcf),'XTick',ddd:2/24:ddd+1);datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        grid;
        subplot(3,2,6);plot(motm(:,1),motm(:,7),'.','markersize',2);ylabel('ratez (rad/s)','FontSize',8);xlabel('Hour (UTC)','FontSize',8);
        ylim([-5 5]*.02);set(gca(gcf),'XTick',ddd:2/24:ddd+1);datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        grid;
        annotation(gcf,'textbox',[0.25 0.96 0.8 0.02462],'String',sprintf('%s (%04i-%02i-%02i, DOY%03i).  PSD motion pack',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','FitBoxToText','off','LineStyle','none','Interpreter','none');
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['Motion_',LegStr],['MotionPak_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % LICOR AGC
        figure; plot(jd_min, licm1(:,5), 'r.','linewidth',1); grid;
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  Licor AGC ',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        xlabel('Hour (UTC)'); ylabel('AGC'); xlim([ddd ddd+1]); ylim([45,100]);
        set(gca(gcf),'XTick',ddd:2/24:ddd+1);datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['Licor_AGC_',LegStr],['AGC_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end
    end

    close all;

    %% save 1-min results to text files, historical format
    % with added columns in scs, rad and met files for scs and WXT variables

    % set lower limit on org rain rate sensitivity
    org_recalc(org_recalc<0.05) = 0;

    % 1-min PSD met data
    % jd pir1 psp1 Tc1 Td1 tsnk ta rh rain org_diag aspir Pmb rwdir_wxt rwspd_wxt ta_wxt rh_wxt Pmb_wxt rain_wxt U_wxt V_wxt
    fpath = fullfile(path_proc_data,['met_',LegStr],[cruise,'_met_1min_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),'.txt']);
    fout = fopen(fpath,'w');
    jazz = zeros(12,1440)*NaN;
    jazz(1,:) = jd_min;
    jazz(2,:) = metm(:,6);      % pir1 W/m2
    jazz(3,:) = metm(:,5);      % psp1 W/m2
    jazz(4,:) = metm(:,9);      % Tc1 C
    jazz(5,:) = metm(:,10);     % Td1 C
    jazz(6,:) = metm(:,4);      % Tsnake C
    jazz(7,:) = metm(:,2);      % Tvais C
    jazz(8,:) = metm(:,3);      % RHvais
    jazz(9,:) = org_recalc;     % org mm/hr
    jazz(10,:) = metm(:,16);    % org carrier
    jazz(11,:) = metm(:,15);    % aspir on
    jazz(12,:) = metm(:,14);	% press mb at zp

    fprintf(fout,['%11.5f ',repmat('%9.2f ',1,12),' \n'],jazz);
    fclose all;

    % 1-min radiation data
    % jd pir1 pir2 psp1 psp2 Tc1 Td1 Tc2 Td2 ta rh tsnk imir imsol impir imtc imtd imir2
    fpath = fullfile(path_proc_data,['rad_',LegStr],[cruise,'_rad_1min_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),'.txt']);
    fout = fopen(fpath,'w');
    jazz = zeros(17,1440)*NaN;
    jazz(1,:) = jd_min;
    jazz(2,:) = metm(:,6);      % pir1 W/m2
    jazz(3,:) = metm(:,8);      % pir2 W/m2
    jazz(4,:) = metm(:,5);      % psp1 W/m2
    jazz(5,:) = metm(:,7);      % psp2 W/m2
    jazz(6,:) = metm(:,9);      % Tc1 C
    jazz(7,:) = metm(:,10);     % Td1 C
    jazz(8,:) = metm(:,11);     % Tc2 C
    jazz(9,:) = metm(:,12);     % Td2 C
    jazz(10,:) = metm(:,2);     % Tvais C
    jazz(11,:) = metm(:,3);     % RHvais
    jazz(12,:) = metm(:,4);     % Tsnake C
    jazz(13,:) = scsm1(:,15);   % imet pir
    jazz(14,:) = scsm1(:,10);   % imet psp
    jazz(15,:) = scsm1(:,20);   % imet thermopile mV
    jazz(16,:) = scsm1(:,21);   % imet T case
    jazz(17,:) = scsm1(:,22);   % imet T dome
    fprintf(fout,['%11.5f ',repmat('%9.3f ',1,17),' \n'],jazz);
    fclose('all');

    % 1-min gps data
    % jd lat lon sog cog hed pitch lat lon head cog sog
    fpath = fullfile(path_proc_data,['gps_',LegStr],[cruise,'_gpsnav_1min_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),'.txt']);
    fout = fopen(fpath,'w');
    jazz = zeros(12,1440)*NaN;
    jazz(1,:) = jd_min;
    jazz(2,:) = gprm1(:,6);     % PSD decimal latitude, deg
    jazz(3,:) = gprm1(:,7);     % PSD decimal longitude, deg
    jazz(4,:) = gprm1(:,3);     % PSD GPS SOG, m/s
    jazz(5,:) = gprm1(:,2);     % PSD GPS COG, deg
    jazz(6,:) = hedm1(:,2);     % PSD GPS heading, deg
    jazz(7,:) = hedm1(:,3);     % PSD GPS angle (pitch), deg
    jazz(8,:) = scsm1(:,2);     % SCS Lat
    jazz(9,:) = scsm1(:,3);     % SCS Lon
    jazz(10,:) = scsm1(:,12);   % SCS Heading
    jazz(11,:) = scsm1(:,4);    % SCS COG
    jazz(12,:) = scsm1(:,5);    % SCS SOG
    fprintf(fout,['%11.5f ',repmat('%9.3f ',1,11),' \n'],jazz);
    fclose('all');

    %% 10-min bulk flux run
    % 10-min averaging was done above
    % using vectorized COARE 3.5 with output adjusted to match historical cor30b
    % y = [hsb hlb tau zo zot zoq L usr tsr qsr dter dqer tkt RF hlwebb Cd Ch Ce Cdn_10 Chn_10 Cen_10 Urf  Trf  RHrf];
    %       1   2   3   4  5   6  7  8   9   10  11   12   13 14   15   16 17 18  19     20     21     22   23   24

    % set lower limit on org rain rate sensitivity
    org_10_recalc(org_10_recalc<0.05) = 0;

    % using PSD data, pressure at sealevel
    y = coare35vn(wspd_10,zu,Tvais_10,zt,Rhvais_10,zq,press_10,Tsea_10,psp1_10,pir1_10,lat_10,600,org_10_recalc,NaN,NaN);

    % using SCS wind and psp only, using PSD pressure, rain, pir, Tair, RH, & SST
    yim = coare35vn(imum_10,zuim,Tvais_10,zt,Rhvais_10,zq,press_10,Tsea_10,imsol_10,pir1_10,latm_10,600,org_10_recalc,NaN,NaN);

    hsb = y(:,1);
    hlb = y(:,2);
    taub = y(:,3);
    tcol = y(:,11);
    qcol = y(:,12);
    rf = y(:,14);
    usr = y(:,8);
    zo = y(:,4);
    Urf = y(:,22);
    Trf = y(:,23);
    RHrf = y(:,24);
    imUrf = yim(:,22);
    imTrf = yim(:,23);
    imRHrf = yim(:,24);

    Rns = psp1_10*0.955;
    Rnl = 0.97*(pir1_10 - 5.67e-8*(Tsea_10+273.15).^4);
    hnet = Rns + Rnl - hsb - hlb - rf;

    hm = nanmean1(hnet); sm = nanmean1(hsb); lm = nanmean1(hlb);
    Sm = nanmean1(Rns); Lm = nanmean1(Rnl); rfm = nanmean1(rf);

    %% 10-min bulk model plots

    if plotit && ~isempty(metm10(isfinite(metm10(:,2:end))))
        % NET HEAT FLUX
        figure; plot(jd_10min, hnet,'b-','linewidth',1); xlabel('Hour (UTC)'); ylabel('Heat flux (W/m^2)'); xlim([ddd ddd+1])
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  PSD Net Heat flux = %04.2f',cruise,Vdate(1),Vdate(2),Vdate(3),ddd,hm),'FontWeight','Bold','Interpreter','none');
        set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks'); grid;
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['Net_Heat_',LegStr],['Net_Heat_Flux_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % HEAT FLUX COMPONENTS - SEA
        figure; plot(jd_10min, -hsb, jd_10min, -hlb, jd_10min, Rnl, jd_10min, Rns, jd_10min, -rf,'g.','LineWidth',1,'MarkerSize',6);
        hold on; plot(jd_10min, hnet,'k','linewidth',2); xlabel('Hour (UTC)');ylabel('Heat Flux (W/m^2)');
        legend(sprintf('H_{sens}= %4.1f',-sm),sprintf('H_{lat}= %4.1f',-lm),sprintf('R_{long}= %4.1f', Lm), sprintf('R_{short}= %4.1f',Sm),...
            sprintf('Rain= %4.1f',-rfm),sprintf('Net= %4.1f',nanmean1(hnet)),'Location','BestOutside'); grid;
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  PSD Heat Flux Components (sea)',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        xlim([ddd ddd+1]); set(gca(gcf),'XTick',ddd:2/24:ddd+1);datetick('x','HH:MM','keepticks');
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['Heat_Flux_Components_',LegStr],['Heat_Flux_Components_sea_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        % BULK FRICTION VELOCITY
        figure; plot(jd_10min, usr,'b','linewidth',1); xlabel('Hour (UTC)'); ylabel('u_*(m/s)'); xlim([ddd ddd+1]);
        title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  PSD Friction Velocity from COARE algorithm 3.5',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
        set(gca(gcf),'XTick',ddd:2/24:ddd+1);datetick('x','HH:MM','keepticks'); grid;
        xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
        annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
        if prtit
            ppath = fullfile(path_raw_images,['u-star_',LegStr],['Friction_Velocity_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
            print(graphdevice,ppath);
        end

        close all;
    end

    %% save 10-min met and bulk model results to txt file - historical format with WXT data at end

    fpath = fullfile(path_proc_data,['flux_',LegStr],[cruise,'_PSD_flux_10min_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),'.txt']);
    fout = fopen(fpath,'w');
    jazz = NaN(46,144);
    jazz(1,:) = jd_10min;            % Decimal julian day at beginning of line
    jazz(2,:) = wspd_10;             % psd true wind speed at zu, m/s
    jazz(3,:) = wdir_10;             % psd true wind direction, deg
    jazz(4,:) = Tsea_10;             % psd seasnake T, C
    jazz(5,:) = tsgm_10;             % ship theromsalinograph T, C  bow
    jazz(6,:) = tssm_10;             % ship theromsalinograph salinity, psu  bow
    jazz(7,:) = Tvais_10;            % psd air T at zt, C
    jazz(8,:) = qsea_10;             % psd air specific humidity at sea surface, g/kg
    jazz(9,:) = qair_10;             % psd air specific humidity, g/kg
    jazz(10,:) = psp1_10;            % psd solar flux, w/m^2
    jazz(11,:) = pir1_10;            % psd IR flux, w/m^2
    jazz(12,:) = org_10_recalc;      % psd org precip rate, mm/hr
    jazz(13,:) = sogm_10;            % ship sog from gps, m/s
    jazz(14,:) = lrgm_10;            % ship heading from gyrocompass, deg
    jazz(15,:) = rwspd_10;           % psd rel wind speed, m/s
    jazz(16,:) = rwdir_10;           % psd rel wind direction, deg
    jazz(17,:) = lat_10;             % decimal latitude, deg
    jazz(18,:) = lon_10;             % decimal longtude, deg
    jazz(19,:) = ones(size(jd_10min))*dtsg;	% Depth of SST sensor used in heat flux calc, m
    jazz(20,:) = sogm_10_std;        % standard deviation of ship speed, m/s
    jazz(21,:) = taub;               % wind stress, coare 3.0, N/m^2
    jazz(22,:) = hsb;                % sensible heat flux, coare 3.5, w/m^2
    jazz(23,:) = hlb;                % latent heat flux, coare 3.5, w/m^2
    jazz(24,:) = rf;                 % rain heat flux, w/m^2
    jazz(25,:) = imta_10;            % IMET air temp, C
    jazz(26,:) = imqa_10;            % IMET air specific humidity, g/kg
    jazz(27,:) = imum_10;            % IMET true wind speed, m/s
    jazz(28,:) = imdm_10;            % IMET true wind direction, deg
    jazz(29,:) = imsol_10;           % IMET solar flux, w/m^2
    jazz(30,:) = imir2_10;           % IMET IR flux, w/m^2
    jazz(31,:) = impress_10;         % IMET BP, mb @ SL
    jazz(32,:) = imrh_10;            % IMET RH (%)
    jazz(33,:) = Urf;                % PSD wspd at 3m from COARE
    jazz(34,:) = Trf;                % PSD ta at 3m from COARE
    jazz(35,:) = RHrf;               % PSD rh at 3m from COARE
    jazz(36,:) = imUrf;              % SCS wspd at 3m from COARE
    jazz(37,:) = imTrf;              % SCS ta at 3m from COARE
    jazz(38,:) = imRHrf;             % SCS rh at 3m from COARE
    jazz(39,:) = rwspd_10_wxt;       % WXT relative wind speed, m/s
    jazz(40,:) = rwdir_10_wxt;       % WXT relative wind direction, deg
    jazz(41,:) = wspd_10_wxt;        % WXT true wind speed, m/s
    jazz(42,:) = wdir_10_wxt;        % WXT true wind direction, deg
    jazz(43,:) = Twxt_10;            % WXT air temp, C
    jazz(44,:) = RHwxt_10;           % WXT RH, %
    jazz(45,:) = Pwxt_10;            % WXT sea level pressure, mb
    jazz(46,:) = rainwxt_10;         % WXT rain rate, mm/hr

    fprintf(fout,['%11.5f ',repmat('%9.3f ',1,45),' \n'],jazz);
    fclose all;
    clear jazz;

    %% radiation models @ 10 min resolution
    % Downward IR
    if all(isnan(qair_10))
        rlclr_10 = (0.52+0.13/60*abs(latm_10)+(0.082-0.03/60.*abs(latm_10)).*sqrt(imqa_10)).*(5.67e-8*(imta_10+273.15).^4);
    else
        rlclr_10 = (0.52+0.13/60*abs(latm_10)+(0.082-0.03/60.*abs(latm_10)).*sqrt(qair_10)).*(5.67e-8*(Tvais_10+273.15).^4);
    end

    % Solar model
    latd = nanmean1(latm_10);
    lond = nanmean1(lonm_10);
    qa = qair_10; % humidity needed for solar model
    qrat = 4.0;
    watvap = iv;
    jdx = jd_10min + double(2.5/1440); % 5-min jd bin centers
    tutc = (jdx - floor(jdx)).*24; % decimal hour of the day, 0-24
    clear sw sz saz dirs sky;
    [sw,sz,saz,dirs,sky] = SolarRadiancex(latd,lond,jdx,tutc,watvap,nanmean1(press_10),k1,k2,oz);
    Rscl_10 = sw;

    %% more plots

    if plotit
        if ~isempty(rlclr_10(isfinite(rlclr_10)))
            ii = find(isfinite(rlclr_10));
            imir_10(imir_10<250) = NaN;  % remove bad ship PIR data
            figure; plot(jd_10min, pir1_10,'b',jd_10min, pir2_10,'r',jd_10min, imir_10,'g',jd_10min(ii),rlclr_10(ii),'k--','linewidth',1); xlim([ddd ddd+1])
            title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  PSD/SCS Infrared Flux',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
            xlabel('Hour (UTC)'); ylabel('IR Flux (W/m^2)'); grid; legend('PIR1','PIR2','SCS','Clearsky','Location','BestOutside');
            set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
            xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00'; set(gca(gcf),'XTickLabel',xax(:,1:2));
            annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
            if prtit
                ppath = fullfile(path_raw_images,['IR_flux_',LegStr],['IRflux_Comparison_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
                print(graphdevice,ppath);
            end
        end

        if ~isempty(Rscl_10(isfinite(Rscl_10)))
            ii = find(isfinite(Rscl_10));
            figure; plot(jd_10min, psp1_10,'b',jd_10min, psp2_10,'r', jd_10min, imsol_10,'g',jd_10min(ii), Rscl_10(ii),'k--','linewidth',1); xlim([ddd ddd+1]);
            title(sprintf('%s (%04i-%02i-%02i, DOY%03i).  PSD/SCS solar flux',cruise,Vdate(1),Vdate(2),Vdate(3),ddd),'FontWeight','Bold','Interpreter','none');
            xlabel('Hour (UTC)');ylabel('Solar Flux (W/m^2)'); grid; legend('PSP1','PSP2','SCS','Clearsky','Location','BestOutside');
            set(gca(gcf),'XTick',ddd:2/24:ddd+1); datetick('x','HH:MM','keepticks');
            xax = get(gca(gcf),'XTickLabel'); xax(end,:)='24:00';set(gca(gcf),'XTickLabel',xax(:,1:2));
            annotation(gcf,'textbox',[0.007154 0.01077 0.4498 0.02462],'String',{'NOAA/ESRL/PSD/Weather & Climate Physics'},'FontSize',6,'FitBoxToText','off','LineStyle','none');
            if prtit
                ppath = fullfile(path_raw_images,['Solar_Flux_',LegStr],['Solarflux_Comparison_',sprintf('%04i_%02i_%02i_%03i',Vdate(1),Vdate(2),Vdate(3),ddd),graphformat]);
                print(graphdevice,ppath);
            end
        end
        close all;
    end

    disp(['FINISHED MANUALflux_eval_',cruise,' yearday ',sprintf('%03i',ddd),'.']);

end
