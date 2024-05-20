function licm = read_licor_day(path_working_ddd,ddd,yyyy,...
                flicor,prtit,path_raw_images,graphformat,graphdevice)
%{
Reads hourly licor files for one day
Loops through 24 hours, calling read_licor for each hour
Raw 10Hz data is averaged to 1Hz output

Inputs: path_working_ddd: path to daily folder containing hourly files
        ddd: day-of-year variable
        yyyy: year string

Output: licm: 86400 x num_vars array of 1-Hz avg gps data
        Returns NaN for data gaps

        1   jd_ref          10 Hz timestamp
        2   Licor_H2O       H2O water vapor density, mmol/m3
        3   Licor_T         box temperature
        4   Licor_P         box pressure kPa
        5   Licor_AGC       diagnostic
        6   Licor_CO2       CO2 vapor density, mmol/m3
        
        % output below is not used or included... because it's calculated
        using noisy and wrong values of Tbox and Pbox. 
        output below is all for dry air only instead of the
        true wet+dry air mixture. And it's not very accurate because it's
        using raw Tbox and Pbox... which are subject to deck and box
        heating and not corrected for dynamic pressure effects. Best
        practice is to recompute it in the driver program using the Vaisala
        RH/T and dynamic pressure measurements. 
        7   Licor_H2O_mr        H2O water vapor mixing ratio, mol/mol
        8   Licor_rhoair_dry    dry air density, kg/m3
        9   Licor_q_dry         dry air H2O specific humidity, g/kg
        10  Licor_CO2_mr        CO2 vapor mixing ratio, ppm or micromol CO2 / mol dry air
                
%}

fclose all;

%%% number of variables to be saved
num_vars = 6;

licm = zeros(86400,10)*NaN;
delta = double(1.0/86400);
last = ddd + 86400*delta;
jd_ref = ddd:delta:last;  jd_ref = jd_ref(:); % ref 1 Hz timestamp
licm(:,1) = jd_ref(1:end-1)';
jd = sprintf('%03i',ddd);

% %%% these are repeats from setup_cruise.m but are not accessible in this
% %%% form of a function, so are repeated here if calculations below are
% used. They are not currently because they use the erroneous values of
% Pbox and Tbox, which have crazy values.
% Rgas = 287.1;
% Rgas_universal = 8.314472; % Pa m3 /K /mol
Mw = 18.01528;  % molar mass of H20 g/mol
% Md = 28.964;    % molar mass dry air
% epsilon = Mw/Md; % mass of water to mass of dry air = 0.622
% C2K = 273.15;  % deg C to Kelvin

for hhh = 0:23               % cycle thru 24 hourly gprm files
    hr = sprintf('%02i',hhh);
    dfl = fullfile(path_working_ddd,['lic0' yyyy(3:4),jd,hr,'_raw.txt']);
    if exist(dfl,'file')==2
        lic = read_licor(dfl,ddd,hhh);
    else
        continue % no data this hour
    end
    jdlic = lic(:,1);

%     % this is neither used in this program, nor is it up to date. It is 
%     % computed in main program instead. 
%     % convert to mol/mol and specific humidity
%     H2O_mr = (lic(:,2)*1e-3).*Rgas_universal.*(lic(:,3)+C2K)./(lic(:,4)*1e3); % mmol/m3 to mol water / mol dry air
%     rhoa_dry = (lic(:,4)*1e3).*(1-H2O_mr*(1-epsilon))./(Rgas.*(lic(:,3)+C2K));% dry air density [kg/m3]
%     q_dry = (lic(:,2)*1e-3).*Mw./rhoa_dry;  % mmol/m3 to g/kg for dry air only
%     CO2_mr = (lic(:,2)*1e3).*Rgas_universal.*(lic(:,3)+C2K)./(lic(:,4)*1e3); % C02 MR converted from mmol/m3 to ppm... micromol C02 / mol dry air
%     lic = [lic, H2O_mr, rhoa_dry, q_dry, CO2_mr];      % add columns to lic array


    % average this hour's data into 1Hz licm array
    start = ddd + hhh/24.0; % jd start time this hour
    diff = jd_ref - start;  % look for closest time stamp to start
    [~,ii] = min(abs(diff));
    temp = interval_avg(jdlic, lic(:,2:num_vars), jd_ref(ii:ii+3600));
    licm(ii:ii+3599,2:num_vars) = temp(:,2:num_vars);

    % for plots below only.... recalculated and saved in driver program
    % manualflux_eval
%     rhoai = interp1(jd_ref, rhoa, jdlic);
%     q = (licm(:,2)*1e-3).*Mw./rhoai;  % mmol/m3 to g/kg for dry+moist air mix

    
    % not sure if this is necessary... it gets redone in runmotcorr I
    % think. try seeing if this program gets a lot faster if this is
    % skipped. it's not necessary for simply producing mean 1-min and 10-min
    % data files and plots. I am not doing it because qa is technically
    % wrong here... it includes Tbox and Pbox which are bad.
    if prtit && mod(hhh,6)==0 % spectra plot
        Npts = length(lic(:,2));
        [Pqq,~] = psd2(detrend(lic(:,2)),Npts,flicor,hamming(Npts));
        [Pqq,Fxx] = specsmoo(Pqq,flicor);

        figure;loglog(Fxx,Fxx.*Pqq,'b-');
        hold on;
        loglog(Fxx(17:48),((Fxx(17:48)).^(-2/3))./(Fxx(35).^(-2/3)).*Pqq(35).*Fxx(35),'r');
        title([jd,' ',hr,' Spectrum of H2O vapor density (mmol m^{-3}, for calculating q) ']); grid;
        ppath = fullfile(path_raw_images,['qSpectra'],['qSpectrum_',jd,'_',hr,graphformat]);
        print(graphdevice,ppath);

        close all;
    end
end
