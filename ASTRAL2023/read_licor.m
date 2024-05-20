function lic = read_licor(dfl,ddd,hhh)
%{
NO CO2 data for PISTON_MISOBOB 2019
reads licor file specified by dfl and returns array with the
following columns:

1   jd_ref      10 Hz timestamp
2   Licor_H2O   H2O water vapor density, mmol/m3
3   Licor_T     Licor box temperature, C
4   Licor_P     Licor box pressure, kPa
5   AGC         Licor diagnostic
6   Licor_CO2   CO2 vapor density, mmol/m3

input parameters: dfl = string path to scs file
                  ddd = julian date
                  hhh = hour
%}

%% reference timestamp
disp(['Reading lic file for hour ',int2str(hhh)]);
start = (ddd+hhh/24);
delta = double(1.0/864000);
last = start + 35999*delta;
tref = start:delta:last;

%% read data if file exists, if not return NaNs
if exist(dfl,'file')==2
    %% read data
    flist = fopen(dfl,'r');
    temp = {};
    z = [];

    while feof(flist)==0
        try    % read file into cell array. temp is 1x8 cell array
            temp = textscan(flist,'%2f%2f%3f %f %f %f %f %f %*[^\n]','delimiter', ', ','headerlines', 1,'emptyvalue',NaN,'treatAsEmpty','Overflow');
            z = [z cell2mat(temp)'];
        catch
            for ii=1:10  % length of last cell reflects missing values from any/all fields
                if ii <= size(temp,2) %%% EJT fix... something wrong with numbers 8 vs. 10 here
                    if length(temp{1,ii})~=length(temp{1,8})
                        temp{1,ii}(length(temp{1,8})+1)=[]; %truncate length of all cells to length(temp{1,10})
                    end
                end
            end
            if ~isempty(temp{1,1})
                z=[z cell2mat(temp)'];
            end
        end
    end
    fclose(flist);
    z(:,any(isnan(z),1))=[]; % kill lines beginning with NaN

    tlic =  ddd + (hhh + (z(1,:)+(z(2,:)+z(3,:)/1000)/60)/60)/24;

    Licor_diag = z(4,:);
    Licor_CO2 = z(5,:); % mmol/m3 differs by factor of 10^3. ppm = micromoles Co2 per mole of moist air
    %%% used for quality controlling licor data - if it's reading 400 ppm,
    %%% that means C02 is working. if co2 is working, then humidity values
    %%% are good. If it's reading 100 ppm, then CO2 is bad and it's
    %%% unlikely that H20 is good. We can recreate it. The mean C02 needs
    %%% to be corrected for temperature and pressure in driver program.
    Licor_H2O = z(6,:); % mmol/m^3
    Licor_T = z(7,:); % C
    Licor_P = z(8,:); % kPa

    %% interpret status diagnostic
    xx = Licor_diag;
    b7 = floor(xx/128);
    y = b7*128;
    b6 = floor((xx-y)/64);
    y = y+b6*64;
    b5 = floor((xx-y)/32);
    y = y+b5*32;
    b4 = floor((xx-y)/16);
    y = y+b4*16;
    AGC = (xx-y)*6.25;

    %% Format output to exactly 10 Hz
    lic = NaN(36000,6);
    lic(:,1) = tref;

    % interpolate to 10 Hz
%     [~,ii,~] = unique(tlic);
%     lic(:,2) = interp1(tlic(ii),Licor_H2O(ii),tref,'linear','extrap');
%     lic(:,3) = interp1(tlic(ii),Licor_T(ii),tref,'linear','extrap');
%     lic(:,4) = interp1(tlic(ii),Licor_P(ii),tref,'linear','extrap');
%     lic(:,5) = interp1(tlic(ii),AGC(ii),tref,'linear','extrap');
%     lic(:,6) = interp1(tlic(ii),Licor_CO2(ii),tref,'linear','extrap');

    % or, copy over original points to the 10 Hz grid - this is slower
    temp = tref;
    for ii = 1:length(tlic)                          % loop over all data points
        if any(isfinite(temp))                       % if there are unused grid points
            [~,jj] = nanmin1(abs(temp - tlic(ii)));  % find closest grid point for this datum
            dT = (tlic(ii)-temp(jj))*86400*1000;     % compute deltaT in mSec
            if dT>=-90 && dT<90                      % only if dT is within range of +/- 90 ms...
                lic(jj,2) = Licor_H2O(ii);           % copy data over to gridded array point
                lic(jj,3) = Licor_T(ii);
                lic(jj,4) = Licor_P(ii);
                lic(jj,5) = AGC(ii);
                lic(jj,6) = Licor_CO2(ii);
                temp(jj) = NaN;                      % then remove this grid point from further consideration
            end
        end
    end

    % fill in NaNs with nearest neighbor
    lic(:,2) = replace_NaN_nearest_neighbor(lic(:,2));
    lic(:,3) = replace_NaN_nearest_neighbor(lic(:,3));
    lic(:,4) = replace_NaN_nearest_neighbor(lic(:,4));
    lic(:,5) = replace_NaN_nearest_neighbor(lic(:,5));
    lic(:,6) = replace_NaN_nearest_neighbor(lic(:,6));

    %% remove spikes in each 10-min segment of H2O - set rejection criteria in despike2.m
    for ii = 1:6000:30001
        [lic(ii:ii+5999,2),~] = despike2(lic(ii:ii+5999,2));
    end
    
%     disp(lic(1,:));
else
    lic = NaN(36000,6);
    lic(:,1) = tref';
end

end
