function met = read_met(dfl1,dfl2,ddd,hhh,zp)
%{
reads met1, met2 files and returns 1 min avg array with the
following columns: (this is the output)

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
    13  org             rain, mm/hr - none for this cruise
    14  slp             air slp, mb
    15  aspir_on        T/RH diagnostic
    16  org_carrier     org diagnostic
    17  org_V           ORG raw measurement, Volts

input parameters: dfl1 = string path to met1 file
                  dfl2 = string path to met2 file
                  ddd = julian date
                  hh = hour
                  adj =  1x7 array of temperature adjustments
                        [tsea,tair,td1,tc1,td2,tc2]
                            *** this might still be a good idea for
                            radiometers, but I don't think this is a good
                            idea to do for sea snake and air temp. Would
                            rather correct after the fact than rerun all
                            the raw data intake. 

    Mean 1 raw file format:
    1:3 PSD timecode: MMSSsss
    4   Batt_Volt
    5   PTemp
    6   ORG_Car_v
    7   ORG_Sig_v
    8   ORG_mmhr
    9   airflow
    10  No data
    11  No data
    12  RH
    13  Tair
    14  Std Dev airflow
    15  No data
    16  No data

    Mean 2 raw file format:
                    ... and actual column number in raw files:
    1:3 PSD timecode: MMSSsss   1
    4   Batt_Volt               2
    5   PTemp                   3
    6   Case1_R                 4
    7   Case1_T                 5
    8   Dome1_R                 6
    9   Dome1_T                 7
    10  Case2_R                 8
    11  Case2_T                 9
    12  Dome2_R                 10
    13  Dome2_T                 11 
    14  PIR1_V                  12   
    15  PIR1_wm2                13    
    16  PSP1_V                  14
    17  PSP1_wm2                15
    18  PIR2_V                  16
    19  PIR2_wm2                17
    20  PSP2_V                  18
    21  PSP2_wm2                19
    22  Snake_R                 20
    23  Snake_C                 21
    24  BP_hpa                  22

%}

%% initialize

g = 0.0098;          % dry adiabatic lapse rate, C/m
sig_sb = 5.67e-8;    % Stefan Boltzmann constant
C2K = 273.15;        % temp conversion constant

% PIR1 s/n 30433F3 calibration on 06/10/2021 
% Const PIR1_Cnt = 0.00339
% 
% 'PIR2 s/n 30432F3 calibration on 06/10/2021 
% Const PIR2_Cnt = 0.00434
% 
% PSP1 s/n 30593F3 calibration on 06/29/2018
% Const PSP1_Cnt = 0.00814
% 
% PSP2 s/n 30434F3 calibration on 06/29/2018
% Const PSP2_Cnt = 0.00794
k1 = 3.63; % coefficients for PIR1 from June 2017 calibration
k2 = 3.37; % coefficients for PIR2 from June 2017 calibration

nr1 = 16;     % number of variables in met1
nr2 = 24;     % number of variables in met2
mean1 = [];   % array for raw met1 data
mean2 = [];   % array for raw met2 data
temp = {};

% reference timestamp and 1-min averaging bins for this hour
start = (ddd+hhh/24);
delta = double(1.0/1440);
last = start + 59*delta;
jd_ref = start:delta:last;
jd_bin = start:delta:last+delta;

%% read met1 file
disp(['Reading met1 file for hour ',int2str(hhh)]);
if exist(dfl1,'file')==2
    flist = fopen(dfl1,'r');
    while ~feof(flist)
        try
            temp = textscan(flist,['%2f%2f%3f ',repmat('%f ',1,13),'%*[^\n]'],'delimiter',', ','headerlines',1);
            mean1 = [mean1 cell2mat(temp)];
        catch
            % if one line is cut short, blank it out totally
            [r,c] = size(temp{1,nr1}); % check size of last cell
            if r==0 || c==0; continue; end % read next batch of lines
            for ll=1:nr1
                if length(temp{1,ll})~=r
                    temp{1,ll}(r+1)=[];
                end
            end
            for ii=1:length(temp)  % Look for NaNs and blanks in data fields
                temp{ii}(ismember(temp{ii},'NAN'))=cellstr('NaN');
                temp{ii}(ismember(temp{ii},''))=cellstr('NaN');
                mean1 = [mean1 str2num(char(temp{ii}))];
            end
        end
    end
    mean1 = mean1';
    jd_pc_1 = ddd+(hhh+(mean1(1,:)+(mean1(2,:)+mean1(3,:)/1000)/60)/60)/24;   % timestamp from PSD system
    fclose(flist);
else
    mean1 = NaN(nr1,60);
    jd_pc_1 = jd_ref;
end

%% read met2 file
disp(['Reading met2 file for hour ',int2str(hhh)]);
if exist(dfl2,'file')==2
    flist = fopen(dfl2,'r');
    while ~feof(flist)
        try
            temp = textscan(flist,['%2f%2f%3f ',repmat('%f ',1,21),'%*[^\n]'],'delimiter',', ','headerlines',1);
            mean2 = [mean2 cell2mat(temp)];
        catch
            
            %%% need to add another catch for this: see met2 day 18 hr 12
%             3756301 ?Э? &H??Э? &H??Э? &H??ь?p? 	  ?????Э? &H??Э? &H??Э? &H??ь?p? 	 ??ях??Э? &H??Э? &H??Э? &H??ь?p? 	 ??9═??Э? &H??Э? &H??Э? &H??ь?p? 	 ??ш??12.00408,27.84956,9353.08,26.63952,9356.822,26.62949,9348.03,26.65302,9349.902,26.64801,0.07131554,24.42313,-2.08768,254.9059,0.06003983,24.20961,-1.535108,194.5638,2825.017,26.35916,1017.25

            
            % if one line is cut short, blank it out totally
            [r,c] = size(temp{1,nr2});
            if r==0 || c==0; continue; end % read next batch of lines
            for ll=1:nr2
                if length(temp{1,ll})~=r
                    temp{1,ll}(r+1)=[];
                end
            end
            for ii=1:length(temp)  % Look for NaNs and blanks in data fields
                temp{ii}(ismember(temp{ii},'NAN'))=cellstr('NaN');
                temp{ii}(ismember(temp{ii},''))=cellstr('NaN');
                mean2 = [mean2 str2num(char(temp{ii}))];  % 60xnr2 array
            end
        end
    end
    mean2 = mean2';
    jd_pc_2 = ddd+(hhh+(mean2(1,:)+(mean2(2,:)+mean2(3,:)/1000)/60)/60)/24;   % timestamp from PSD system
    fclose(flist);
else
    mean2 = NaN(nr2,60);
    jd_pc_2 = jd_ref;
end

%% format variables
slp         = mean2(24,:)+0.125*zp;  % sea level pressure (mbar)
pa          = mean2(24,:);           % ambient pressure (mbar)
psp1        = mean2(17,:);           % solar radiometers, should be positive value
psp2        = mean2(21,:);           % solar radiometers, should be positive value
org_carrier = mean1(6,:);            % org diagnostic
org_V       = mean1(7,:);            % raw org voltage measurement
org         = mean1(8,:);            % rain rate
aspir_on    = mean1(9,:);
Tsea        = mean2(23,:);           % sea snake
Tc1         = mean2(7,:);            % case temperature for PIR1 
Td1         = mean2(9,:);            % dome temperature for PIR1 
Tc2         = mean2(11,:);           % case temperature for PIR2 
Td2         = mean2(13,:);           % dome temperature for PIR2
therm1      = -mean2(15,:);           % thermopile PIR1 outputs W/m2, should be negative value since they measure net IR, which is negative in tropics
therm2      = -mean2(19,:);           % thermopile PIR2 outputs W/m2, should be negative value since they measure net IR, which is negative in tropics
Rhvais      = mean1(12,:);           % Rh vaisala
Tvais       = mean1(13,:);           % Ta vaisala
% net solar should be positive heating the ocean, net ir should be negative cooling ocean. 
% pyranometer will give positive number (positive voltage) when sun shines on it. 
% if clear sky overhead, IR will measure a negative value because the net radiation
% is negative - earth is losing more IR radiation then incoming. Before,
% there would be wires hooked up manually, and you never knew which way it
% got hooked up. Just wait for sunny day. Look for voltages. If you put
% hand over IR system, it should should be positive because you radiate at
% 98.6 F and solar should become less positive because you shaded it. Dome
% T should go up rapidly when you touch it, case T should stay the same.
% Sergio has hardwired the systems so we shouldn't have to deal with sign
% ambiguities. Unless there is an issue in campbell software. For instance,
% in PISTON, all the radiation values are backwards sign. 

%% compute IR radiation - no PIR2 for MISOBOB
% This is how you would adjust Epply PIR's for cal constants in June 2017 ESRL calibration
% CR2 logger file values were 3.25 & 4.16 for k1 and k2, which are the c3
% coefficient for PIR #1 and #2 units.
% pir1 = therm1 + sig_sb*(Tc1+adj(4)+C2K).^4 - k1*sig_sb*((Td1+adj(3)+C2K).^4-(Tc1+adj(4)+C2K).^4);  %% how it was done in the field
% pir2 = therm2 + sig_sb*(Tc2+adj(6)+C2K).^4 - k2*sig_sb*((Td2+adj(5)+C2K).^4-(Tc2+adj(6)+C2K).^4);  %% how it was done in the field

% downwelling IR = (-) measured positive net IR + upwelling IR absolute value - dome heating effect.
pir1 = therm1 + sig_sb*(Tc1+C2K).^4 - k1*sig_sb*((Td1+C2K).^4-(Tc1+C2K).^4);
pir2 = therm2 + sig_sb*(Tc2+C2K).^4 - k1*sig_sb*((Td2+C2K).^4-(Tc2+C2K).^4);

%% filter out bad values in T, RH and P
Rhvais(Rhvais>110) = NaN;
Rhvais = replace_NaN_nearest_neighbor(Rhvais);
[Rhvais, ~] = despike2(Rhvais);

Tvais(Tvais>45) = NaN;
Tvais = replace_NaN_nearest_neighbor(Tvais);
[Tvais, ~] = despike2(Tvais);

slp = replace_NaN_nearest_neighbor(slp);
[slp,~] = despike2(slp);

pa = replace_NaN_nearest_neighbor(pa);
[pa,~] = despike2(pa);

%% interp to exact 1-min timestamp and format output
met = NaN(60,20);
met(:,1) = jd_ref;
met(:,2) = interp1(jd_pc_1',Tvais,jd_ref','nearest','extrap');
met(:,3) = interp1(jd_pc_1',Rhvais,jd_ref','nearest','extrap');
met(:,4) = interp1(jd_pc_2',Tsea,jd_ref','nearest','extrap');
met(:,5) = interp1(jd_pc_2',psp1,jd_ref','nearest','extrap');
met(:,6) = interp1(jd_pc_2',pir1,jd_ref','nearest','extrap');
met(:,7) = interp1(jd_pc_2',psp2,jd_ref','nearest','extrap');
met(:,8) = interp1(jd_pc_2',pir2,jd_ref','nearest','extrap');
met(:,9) = interp1(jd_pc_2',Tc1,jd_ref','nearest','extrap');
met(:,10) = interp1(jd_pc_2',Td1,jd_ref','nearest','extrap');
met(:,11) = interp1(jd_pc_2',Tc2,jd_ref','nearest','extrap');
met(:,12) = interp1(jd_pc_2',Td2,jd_ref','nearest','extrap');
met(:,13) = interp1(jd_pc_1',org,jd_ref','nearest','extrap');
met(:,14) = interp1(jd_pc_2',slp,jd_ref','nearest','extrap');
met(:,15) = interp1(jd_pc_1',aspir_on,jd_ref','nearest','extrap');
met(:,16) = interp1(jd_pc_1',org_carrier,jd_ref','nearest','extrap');
met(:,17) = interp1(jd_pc_1',org_V,jd_ref','nearest','extrap');
met(:,18) = interp1(jd_pc_2',therm1,jd_ref','nearest','extrap');
met(:,19) = interp1(jd_pc_2',therm2,jd_ref','nearest','extrap');
met(:,20) = interp1(jd_pc_2',pa,jd_ref','nearest','extrap');

end
