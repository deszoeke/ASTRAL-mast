function metm = read_met_day(path_working_ddd,ddd,yyyy,zp)
%{
Reads hourly met files for one day
Loops through 24 hours, calling read_met for each hour.

Inputs: path_working_ddd: path to daily folder containing hourly files
        ddd: day-of-year variable
        yyyy: year string
        adj: temperature adjustmemts array

Output: metm: 1440x20 array of bulk met data columns
        Returns NaN for data gaps

    1   jd_met          1 min timestamp, decimal DOY @ start of interval
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
    14  press           psd air pressure, mb
    15  aspir_on        T/RH diagnostic
    16  org_carrier     org diagnostic
    17  org_V           ORG signal voltage
    18  therm1          pir1 thermopile output, W/m2
    19  therm1          pir2 thermopile output, W/m2
    20  P               measured pressure, mb
%}

fclose all;

metm = zeros(1440,20)*NaN;
delta = double(1.0/1440);
last = ddd + 1439*delta;
jd_ref = ddd:delta:last;	% ref 1-min timestamp
metm(:,1) = jd_ref';
jd = sprintf('%03i',ddd);

for hhh = 0:23              % cycle thru 24 hourly met files
    hr = sprintf('%02i',hhh);
    dfl1 = fullfile(path_working_ddd,['me1' yyyy(3:4),jd,hr,'_raw.txt']);
    dfl2 = fullfile(path_working_ddd,['me2' yyyy(3:4),jd,hr,'_raw.txt']);
    % read this hours 1-min avg data from both systems if files exist
    % met is 60x11 array
    met = read_met(dfl1,dfl2,ddd,hhh,zp);
    jdmet = met(:,1);

    % copy this hours data into metm array
    for jj = 1:length(jdmet)
        diff = jd_ref - jdmet(jj);
        [val,idx] = min(abs(diff));
        % allow 2 sec tolerance for timestamp match
        if val<double(2.0/86400) && length(idx)==1
            metm(idx,2:20) = met(jj,2:20);
        end
    end
end

end
