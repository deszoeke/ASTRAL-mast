function hedm = read_hed_day(path_working_ddd,ddd,yyyy)
%{
Reads hourly heading files for one day
Loops through 24 hours, calling read_hed for each hour.
Raw 10Hz data is averaged to 1Hz output

Inputs: path_working_ddd: path to daily folder containing hourly files
        ddd: day-of-year variable
        yyyy: year string

Output: hedm: 86400x3 array of heading and angle
        Returns NaN for data gaps

    1   jd_ref      10 Hz timestamp
    2   hed         heading
    3   roll/pitch  angle

Choose appropriate line below for roll or pitch setup
%}

fclose all;

hedm = zeros(86400,3)*NaN;
delta = double(1.0/86400);
last = ddd + 86400*delta;
jd_ref = ddd:delta:last;	% ref 1 Hz timestamp
hedm(:,1) = jd_ref(1:end-1)';
jd = sprintf('%03i',ddd);

for hhh = 0:23             % cycle thru 24 hourly gprm files
    hr = sprintf('%02i',hhh);
    dfl = fullfile(path_working_ddd,['hed0' yyyy(3:4),jd,hr,'_raw.txt']);
%    hed = read_hed_pitch(dfl,ddd,hhh);
    hed = read_hed_roll(dfl,ddd,hhh); 

    jdhed = hed(:,1);
    % average this hour's data into 1Hz hedm array
    start = ddd + hhh/24.0; % jd start time this hour
    diff = jd_ref - start;  % look for closest time stamp to start
    [~,ii] = min(abs(diff));
    temp = interval_avg(jdhed, hed(:,2:3), jd_ref(ii:ii+3600)');
    hedm(ii:ii+3599,2:3) = temp(:,2:3);
end

