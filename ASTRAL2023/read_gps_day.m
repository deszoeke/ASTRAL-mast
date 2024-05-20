function gprm = read_gps_day(path_working_ddd,ddd,yyyy,PosLims)
%{
Reads hourly gprm files for one day
Loops through 24 hours, calling read_gprm for each hour
Raw 10Hz data is averaged to 1Hz output

Inputs: path_working_ddd: path to daily folder containing hourly files
        ddd: day-of-year variable
        yyyy: year string

Output: gprm: 86400 x 7 array of 1-Hz avg gps data
        Returns NaN for data gaps

        1   jd_ref      1 Hz timestamp
        2   cog         cog, deg
        3   sog         sog, kts
        4   Ngps        filtered N component of ship speed, m/s
        5   Egps        filtered E component of ship speed, m/s
        6   lat         latitude (N positive)
        7   lon         longitude (0-360 deg, east)
%}

fclose all;

gprm = zeros(86400,7)*NaN;
delta = double(1.0/86400);
last = ddd + 86400*delta;
jd_ref = ddd:delta:last;	% ref 1 Hz time bins
gprm(:,1) = jd_ref(1:end-1)';
jd = sprintf('%03i',ddd);

for hhh = 0:23              % cycle thru 24 hourly gprm files
% for hhh = 0:23              % cycle thru 24 hourly gprm files
    hr = sprintf('%02i',hhh);
    dfl = fullfile(path_working_ddd,['gprm' yyyy(3:4),jd,hr,'_raw.txt']);
    gps = read_gprm(dfl,ddd,hhh,PosLims);
    jdgps = gps(:,1);
    % average this hour's data into 1Hz gprm array
    start = ddd + hhh/24.0; % jd start time this hour
    diff = jd_ref - start;  % look for closest time stamp to start
    [~,ii] = min(abs(diff));
    temp = interval_avg(jdgps, gps(:,2:7), jd_ref(ii:ii+3600)');
    gprm(ii:ii+3599,2:7) = temp(:,2:7);
end
