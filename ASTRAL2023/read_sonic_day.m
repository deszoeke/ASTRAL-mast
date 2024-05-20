function [sonm,QC] = read_sonic_day(path_working_ddd,ddd,yyyy,sonicmodel,...
         rotationsonic,fsonic,prtit,path_raw_images,graphformat,graphdevice)
%{
Reads hourly sonic files for one day
Loops through 24 hours, calling read_sonic for each hour
Raw 10Hz data is averaged to 1Hz output

Inputs: path_working_ddd: path to daily folder containing hourly files
        ddd: day-of-year variable
        yyyy: year string
        sonicmodel: string (WindMasterPro, R3, R3A, R2, R2A)
        rotationsonic: true for 30 deg rotation on some models
        fsonic: sampling frequency, Hz
        prtit: set 'true' to save spectra plots
        path_raw_images: path for saving plots
        graphformat: for plots, e.g. .png
        graphdevice: print driver, e.g. -dpng

Output: sonm: 86400 x 5 array of 1Hz avg data
        Returns NaN for data gaps

        1   jd_ref      10 Hz timestamp
        2   U           x axis acceleration, m/s
        3   V           y axis acceleration, m/s
        4   W           z axis acceleration, m/s
        5   Tsonic      x axis rotation rate, rad/s

        QC array, 144 x 2, # bad, # missing per 10 min

%}

fclose all;

sonm = zeros(86400,5)*NaN;
delta = double(1.0/86400);
last = ddd + 86400*delta;
jd_ref = ddd:delta:last; jd_ref=jd_ref';	% ref 1 Hz timestamp
sonm(:,1) = jd_ref(1:end-1)';
jd = sprintf('%03i',ddd);
QC = zeros(144,2)*NaN;      % for saving count of bad/missing points

for hhh = 0:23               % cycle thru 24 hourly sonm files
    hr = sprintf('%02i',hhh);
    start = ddd + hhh/24.0; % jd start time this hour
    dfl = fullfile(path_working_ddd,['son0' yyyy(3:4),jd,hr,'_raw.txt']);
    [son,badPnts,missing] = read_sonic(dfl,ddd,hhh,...
                         sonicmodel,rotationsonic);
    jdson = son(:,1);
    % average this hour's data into 1Hz array

    diff = jd_ref - start;  % look for closest time stamp to start
    [~,ii] = min(abs(diff));
    temp = interval_avg(jdson, son(:,2:5), jd_ref(ii:ii+3600));
    sonm(ii:ii+3599,2:5) = temp(:,2:5);
    QC(1+hhh*6:1+hhh*6+5,:) = [badPnts',missing'];

    if prtit && mod(hhh,6)==0 && all(isfinite(son(:,2))) % a few spectra plots
        Npts = length(son(:,2));
        [Puu,~] = psd2(detrend(son(:,2)),Npts,fsonic,hamming(Npts));
        [Puu,Fxx] = specsmoo(Puu,fsonic);
        [Pww,~] = psd2(detrend(son(:,4)),Npts,10,hamming(Npts));
        [Pww,~] = specsmoo(Pww,fsonic);
        [Ptt,~] = psd2(detrend(son(:,5)),Npts,10,hamming(Npts));
        [Ptt,~] = specsmoo(Ptt,fsonic);

        figure;loglog(Fxx,Fxx.*Puu,'b-');
        hold on;
        loglog(Fxx(17:48),((Fxx(17:48)).^(-2/3))./(Fxx(35).^(-2/3)).*Puu(35).*Fxx(35),'r');
        title([jd,' ',hr,' Spectrum of Velocity U']); grid;
        ppath = fullfile(path_raw_images,['uSpectra'],['uSpectrum_',jd,'_',hr,graphformat]);
        print(graphdevice,ppath);

        figure;loglog(Fxx,Fxx.*Pww,'b-');
        hold on;
        loglog(Fxx(17:48),((Fxx(17:48)).^(-2/3))./(Fxx(35).^(-2/3)).*Pww(35).*Fxx(35),'r');
        title([jd,' ',hr,' Spectrum of Velocity W']); grid;
        ppath = fullfile(path_raw_images,['wSpectra'],['wSpectrum_',jd,'_',hr,graphformat]);
        print(graphdevice,ppath);

        figure;loglog(Fxx,Fxx.*Ptt,'b-');
        hold on;
        loglog(Fxx(17:48),((Fxx(17:48)).^(-2/3))./(Fxx(35).^(-2/3)).*Ptt(35).*Fxx(35),'r');
        title([jd,' ',hr,' Spectrum of Sonic Temperature']); grid;
        ppath = fullfile(path_raw_images,['tSpectra'],['tSpectrum_',jd,'_',hr,graphformat]);
        print(graphdevice,ppath);

        close all;
    end
end

if prtit
    jd10 = ddd:10/1440:ddd+1-10/1440;
    figure; plot(jd10,QC(:,1),'bo-',jd10,QC(:,2),'r.-');
    title([jd,' Sonic Bad and Missing Counts / 10min']); grid;
    legend('Bad Pnts','Missing','location','best');
    ppath = fullfile(path_raw_images,['sonicQC'],['SonicQC_',jd,graphformat]);
    print(graphdevice,ppath);
    close all;

    jd10bin = ddd:10/1440:ddd+1; jd10bin = jd10bin';
    rwspd = sqrt(sonm(:,2).^2 + sonm(:,3).^2);
    temp = interval_avg(sonm(:,1), rwspd, jd10bin);
    figure; plot(temp(:,1),temp(:,2),'bo-'); grid;
    title([jd,' Sonic Relative Wind Speed, 10min']);
    xlabel('DOY'); ylabel('Rel Wind Spd, m/s');
    ppath = fullfile(path_raw_images,['Rwspd'],['Rwspd_',jd,graphformat]);
    print(graphdevice,ppath);
    close all;
end

end

