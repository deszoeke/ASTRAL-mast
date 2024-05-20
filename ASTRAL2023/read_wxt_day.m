function out = read_wxt_day(path_working_ddd,ddd,yyyy,zwxt)
%{
Reads hourly wxt files for one day
Loops through 24 hours, calling read_wxt for each hour
Raw data is averaged to 1-min output

Inputs: path_working_ddd: path to daily folder containing hourly files
        ddd: day-of-year variable
        yyyy: year string

Output: gprm: 1440 x 9 array of 1-min avg wxt data
        Returns NaN for data gaps

        1  jd_ref         1 min timestamp
        2  rwdir_wxt       rel wind direction ,0-360 deg from bow
        3  wspd_wxt        wind speed, m/s
        4  Tair_wxt        air temperature, C
        5  RH_wxt          RH, %
        6  Pmb_wxt         air pressure, mb or hPa
        7  rain_wxt        rain rate, mm/hr
        8  U_wxt           U wind component, WXT
        9  V_wxt           V wind component, WXT
%}

fclose all;

out = NaN(1440,9);
delta = double(1.0/1440);
last = ddd + 1440*delta;
jd_ref = ddd:delta:last;	% ref 1 min time bins
out(:,1) = jd_ref(1:end-1)';
jd = sprintf('%03i',ddd);

for hhh = 0:23              % cycle thru 24 hourly gprm files
    hr = sprintf('%02i',hhh);
    dfl = fullfile(path_working_ddd,['wxt0' yyyy(3:4),jd,hr,'_raw.txt']);
    if exist(dfl,'file')
        wxt = read_wxt(dfl,ddd,hhh,zwxt);
        jdgps = wxt(:,1);
        
        % copy this hours data into output array
        for jj = 1:length(jdgps)
            diff = jd_ref - jdgps(jj);
            [val,idx] = min(abs(diff));
            % allow 2 sec tolerance for timestamp match
            if val<double(2.0/86400) && length(idx)==1
                out(idx,2:9) = wxt(jj,2:9);
            end
        end
    end
end
