% Correct Riegl for doppler shift following Collins et al. 2017

%get the ship heading and speed, wave direction
load('/Volumes/NOAA_Ldisk/ATOMIC_2020/RHB/flux/Processed_Images_motcorr3_ok/riegl_wave_plots/ATOMIC2020_wavedirection_nshipnav.mat')
thetar= abs(head-Wave_dir_composite); %relative angle between ship sheading and wave direction
thetar(thetar>180)=360-thetar(thetar>180);
da=load('/Volumes/NOAA_Ldisk/ATOMIC_2020/RHB/flux/Processed_motcorr3_ok/ATOMIC_2020_da_decorr.txt');
yday = da(:,1);  
sig_h   =   da(:,113); % std dev of heading, >5 indicates ship maneuver
sig_uim =   da(:,98);  % std dev of ship speed, m/s (>.2 indicates maneuver)
lat=da(:,151);
u10n=ncread('/Volumes/NOAA_Ldisk/ATOMIC_2020/RHB/flux/Processed/EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc','wspd_10N');

%%
sig_uim_lim=0.2;sig_h_lim=5;
figure('Position', [10 10 1200 800]);
subplot(3,3,1);%plot(yday,sog_ship,'b.:');
    errorbar(yday,sog_ship,sig_uim,'b+:');
    hold on;errorbar(yday(sig_uim>sig_uim_lim),sog_ship(sig_uim>sig_uim_lim),sig_uim(sig_uim>sig_uim_lim),'ro');
    axis([ddd+hhh/24 ddd+(hhh+50/60)/24 0 10])
    title('ship speed (m/s)');grid;
subplot(3,3,4);plot(yday,wdir,'r.:',yday,Wave_dir_composite,'b.:')
    hold on;errorbar(yday,head,sig_h,'k+:');
    hold on;errorbar(yday(sig_h>sig_h_lim),head(sig_h>sig_h_lim),sig_h(sig_h>sig_h_lim),'ro');
    legend('wind direction (from)','wave direction (from)','ship heading (to)','Location','best');
    axis([ddd+hhh/24 ddd+(hhh+50/60)/24 0 360]);grid;
    title('ship heading & wave/wind direction (deg)')
subplot(3,3,7);plot(yday,thetar,'b.:',[ddd+hhh/24 ddd+(hhh+1)/24],[90 90],'k--')
    axis([ddd+hhh/24 ddd+(hhh+50/60)/24 0 180]);grid;
    title('relative angle heading - wave (deg)')
    text(ddd+(hhh+0.9)/24,90,'with waves','Rotation',90,'FontWeight','Bold') 
    text(ddd+(hhh+0.9)/24,0,'into waves','Rotation',90,'FontWeight','Bold') 

% add spectra
%raw unfiltered
[Sw,F] = psd2(detrend(h),1024,f_mot,win); %Swraw(kk,:) = Sw;
[Sws,Fsw,Df] = specsmoo4wave(Sw,f_mot); %Sws1(kk,:) = Sws;    
ff = find(Fsw>0.01,1,'first');
subplot(3,3,[2,3,5,6,8,9]);
loglog(Fsw(ff:end)',Sws(ff:end),'bo-');
hold on;loglog(Fsw(ff:end)',1.5e-4.*Fsw(ff:end)'.^(-4),'k--');

[Sw,F] = psd2(detrend(h2),1024,f_mot,win); %Swraw(kk,:) = Sw;
[Sws,Fsw,Df] = specsmoo4wave(Sw,f_mot); %Sws1(kk,:) = Sws;            
hold on; loglog(Fsw(ff:end)',Sws(ff:end),'ro-');
% ylim([max(min(Sws(ff:end)),10^-6),10^(ceil(max(log10(Sws(ff:end)))))]);
xlabel('Frequency (Hz)'); ylabel('S_{hh} (m^2/Hz)'); grid;
title(['Wave Spectrum jd ',sprintf('%03i',ddd),'; hr ',sprintf('%02i',hhh), '; min:', sprintf('%02i',floor(jj/120))]);


g=9.81;
segind=find((abs(yday-jdseg(1)))==min(abs(yday-jdseg(1))));
if isnan(thetar(segind))
    thetar(segind)=thetar(segind-1);
end

try
    [fincombined, fobscombined, fcr, fin1, fobs1] = Collins_solutions(sog_ship(segind),thetar(segind),g,Fsw);

dfin=[diff(fin1)];dfin=[dfin(1) dfin];
dfobs=[diff(fobs1)];dfobs=[dfobs(1) dfobs];
try
    Sin=Sws(Fsw<fcr(:,2)).*dfobs./dfin;
catch
    Sin=Sws(Fsw<=fcr(:,2)).*dfobs./dfin;
end

% hold on;loglog(fin1,Sin,'m')
% [maxSws,idx] = max(Sin);
% [fincombined(idx) 1./fincombined(idx)]

[hf]=freqz(bhigh,ahigh,fin1,f_mot);  %Frequency response of digital filter
hold on;loglog(fin1,Sin.*abs(hf),'g')
% ylim([max(min(Sws(ff:end)),10^-6),10^(ceil(max(log10(Sin.*abs(hf)))))]);
ylim([1e-3 500])
% legend('raw','f^{-4}','raw high-pass filtered','DC corrected')

orient landscape

catch
    fcr = NaN(1,2);
    fobs1 = NaN(1,length(Fsw));
    fin1 = NaN(1,length(Fsw));
    fincombined = NaN(1,length(Fsw));
    fobscombined= NaN(1,length(Fsw));
    Sin= NaN(1,length(Fsw));
end

% ppath = fullfile(path_wave_plots,'DC_correction',['ship_motion_diagnostics_and_spectra_',jd,'_',hr,'_',sprintf('%02i',floor(jj/120)), '.png']);
% print('-dpng',ppath);   
% 
% close all
