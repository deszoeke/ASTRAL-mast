function [son,badPnts,missing] = read_sonic(dfl,ddd,hhh,sonicmodel,rotationsonic)
%{
reads son file specified by dfl and returns array with the
following columns:

1   jd_ref      10 Hz timestamp
2   U           bow-stern axis wind velocity, m/s
3   V           port-starboard axis wind velocity, m/s
4   W           z axis wind velocity, m/s
5   Tsonic      sonic temperature, C

also returns the number of bad data points (NaNs) and the number
of missing data lines for each 10-min data segment

input parameters: dfl = string path to scs file
                  ddd = julian date
                  hhh = hour
                  sonicmodel = string parameter
                  rotationsonic = boolean
%}

%% reference timestamp
disp(['Reading son file for hour ',int2str(hhh)]);
start = (ddd+hhh/24);
delta = double(1.0/864000);
last = start + 35999*delta;
tref = start:delta:last;  % exact 10Hz timestamp

%% read data if file exists, if not return NaNs
if exist(dfl,'file')==2
    %% read file
    flist = fopen(dfl,'r');
    st1 = [];               % data array
    temp = {};              % cell array for textscan output
    badPnts = zeros(1,6);   % initial value for # bad points in each 10-min seg
    missing = zeros(1,6);   % initial value for # missing lines in each 10-min seg
    
    while feof(flist)==0   % read entire file into cell array
        try
            % temp will be 1x8 cell array, cells will be ~36000x1 arrays
%           temp = textscan(flist,'%2f%2f%3f %*2c %f %f %f %*c %f %*[^\n]','delimiter',', ','headerlines',1,'emptyvalue',NaN,'treatAsEmpty','Sonic');
            temp = textscan(flist,'%2f%2f%3f %*2c %f %f %f %*c %f %f %*[^\n]','delimiter',', ','headerlines',1,'emptyvalue',NaN,'treatAsEmpty','Sonic');

            st1 = [st1 cell2mat(temp)'];  % 8x~36000 array
        catch
            for ii=1:8  % length of last cell reflects missing values from any/all fields
                if length(temp{1,ii})~=length(temp{1,8})
                    temp{1,ii}(length(temp{1,8})+1) = [];  %truncate length of all cells to length(temp{1,8})
                end
            end
            if ~isempty(temp{1,1})
                st1 = [st1 cell2mat(temp)']; % append truncated cells to st1 and return for more
            end
        end
    end
    fclose (flist);
    tsons =  ddd + (hhh + (st1(1,:)+(st1(2,:)+st1(3,:)/1000)/60)/60)/24;

    %% count missing data lines - gaps in timeseries
    %  negative values mean a few extra points...which is not a problem
    segBins = (ddd+hhh/24):10/1440:(ddd+(hhh+1)/24); % 10-min bin edges
    for ii=1:6
        zz = tsons>=segBins(ii) & tsons < segBins(ii+1) & isfinite(tsons);
        missing(ii) = 6000 - sum(zz); % points missing in each 10min segment
    end

    %% clean data
    u = st1(4,:);
    v = st1(5,:);
    W = st1(6,:);
    Tsonic = st1(7,:);

    % Count NaNs in each 10-min segment
    % i.e. when ice or water interferes with sonic transducers...
    first = (ddd+hhh/24);
    last = first + 50/1440;
    delta10 = 1.0/144;
    st = (first:delta10:last);
    cc = zeros(1,6);
    dd = zeros(1,6);
    if (all(isfinite(tsons)))
        for bb=1:6     % get indices in tsons for start/end of each 10-min seg
            cc(bb) = find(tsons > st(bb),1,'first');
            dd(bb) = find(tsons < st(bb)+10/1440,1,'last');
            badPnts(bb) = sum(isnan(u(cc(bb):dd(bb)))); % sum Trues (==1) from isnan
        end
    end

    % Remove all NaNs so interp will fill in the gaps
    ee = isfinite(u);
    tsons = tsons(ee);
    u = u(ee); v = v(ee); W = W(ee); Tsonic = Tsonic(ee);

    %% Apply rotation as required by particular model and mounting config
    if rotationsonic % true when sonic coordinates are rotated by 30 deg
        switch sonicmodel
            case {'R3','WindMasterPro'}     % omnidirectional  (N symbol pointing backward!)
                U = -u*cos(30/180*pi)+v*sin(30/180*pi);     %to North
                V = -u*sin(30/180*pi)-v*cos(30/180*pi);     %to West
            case 'R3A'     % asymmetric      (N symbol pointing forward!)
                U = u*cos(30/180*pi)-v*sin(30/180*pi);      %to North
                V = u*sin(30/180*pi)+v*cos(30/180*pi);      %to West
            case 'R2'      % omnidirectional  (N symbol pointing backward!)
                U = u*cos(30/180*pi)-v*sin(30/180*pi);      %to North
                V = u*sin(30/180*pi)+v*cos(30/180*pi);      %to West
            case 'R2A'     % asymmetric      (N symbol pointing forward!)
                U = -u*cos(30/180*pi)+v*sin(30/180*pi);     %to North
                V = -u*sin(30/180*pi)-v*cos(30/180*pi);     %to West
        end
    else
        switch sonicmodel
            case 'WindMasterPro'	% (N symbol pointing backward!)
                U = -u; % to North
                V = -v; % to West
            case 'R3'               % omnidirectional (N symbol pointing backward!)
                U = -u; % to North
                V = -v; % to West
            case 'R3A'              % asymmetric (N symbol pointing forward!)
                U = u; % to North
                V = v; % to West
        end
    end

    % correct error in Winmaster Pro W velocity, Gill formula
    if strcmp(sonicmodel,'WindMasterPro')
        ppp = W>0; W(ppp) = W(ppp)*1.166;
        nnn = W<0; W(nnn) = W(nnn)*1.289;
    end

    %% Format output to exactly 10 Hz

    % Interpolate to 10Hz grid, linear interp over missing and NaN data
%     [~,zz,~] = unique(tsons);
%     su = interp1(tsons(zz),U(zz),tref,'linear','extrap');
%     sv = interp1(tsons(zz),V(zz),tref,'linear','extrap');
%     sw = interp1(tsons(zz),W(zz),tref,'linear','extrap');
%     st = interp1(tsons(zz),Tsonic(zz),tref,'linear','extrap');

    % or, copy over original points to the 10 Hz grid - this is slower
    temp = tref;
    su = NaN(1,length(tref)); sv = su; sw = su; st = su;
    for ii = 1:length(tsons)                         % loop over all data points
        if any(isfinite(temp))                       % if there are unused grid points
            [~,jj] = nanmin1(abs(temp - tsons(ii))); % find closest grid point for this datum
            dT = (tsons(ii)-temp(jj))*86400*1000;    % compute deltaT in mSec
            if dT>=-90 && dT<90                      % only if dT is within range of +/- 90 ms...
                su(jj) = U(ii);                      % copy data over to gridded array point
                sv(jj) = V(ii);
                sw(jj) = W(ii);
                st(jj) = Tsonic(ii);
                temp(jj) = NaN;                      % then remove this grid point from further consideration
            end
        end
    end

    %% despike - set rejection criteria in despike2.m
    % This will also fill all data gaps with nearest neighbor
    % Use the missing/bad data count to exclude 10-min periods when data
    % loss is unacceptable and gaps are to frequent/large
    for ii = 1:6000:30001  % remove spikes in each 10-min segment
       [su(ii:ii+5999),~] = despike2(su(ii:ii+5999));
       [sv(ii:ii+5999),~] = despike2(sv(ii:ii+5999));
       [sw(ii:ii+5999),~] = despike2(sw(ii:ii+5999));
       [st(ii:ii+5999),~] = despike2(st(ii:ii+5999));
    end
    son = [tref' su' sv' sw' st'];
else
    son = NaN(36000,5);
    son(:,1) = tref';
    badPnts = zeros(1,6);
    missing = ones(1,6)*6000;   % all data missing
end

end

