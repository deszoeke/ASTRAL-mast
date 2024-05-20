function gps = read_gprm(dfl,ddd,hhh,PosLims)
%{
reads gprm file specified by dfl and returns array with the
following columns:

1   tref      10 Hz timestamp
2   cog         cog, deg
3   sog         sog, m/s
4   sogN        filtered N component of ship speed, m/s
5   sogE        filtered E component of ship speed, m/s
6   Latitude
7   Longitude

input parameters: dfl = string path to scs file
                  ddd = julian date
                  hhh = hour

gps files are reformatted to gprm using a python script.
GPRMC lines are stripped from the raw file and saved to a new gprm file.

Format for $GPRMC lines in the raw file:

0000078 $GPRMC,000000.00,A,5822.95924260,N,04515.92776422,W,1.58,223.96,151013,22.0,W,A*04
1 2 3                      4             5 6              7 8    9      10     11
PSD timecodes              lat             lon              spd  hdg

PYTHON code:

"""HiWinGS PSD gps file parse: parseGpsFiles.py
Parses data from Hemisphere gps files saved on the PSD data system.
saves GPRMC lines with the following format to an output file.
0000078 $GPRMC,000000.00,A,5822.95924260,N,04515.92776422,W,1.58,223.96,151013,22.0,W,A*04

source file names follow the convention: gps0yydddhh_raw.txt
output file names follow the convention: gprmyydddhh_raw.txt

Using Anaconda iPython installation.
Versions:anaconda 1.6.0, python 2.7.5

Byron Blomquist, NOAA/ESRL/PSD3, Aug 2014
"""
import os

# assume current directory contains raw gps files,
# dump list of data file names into a list variable
files = [f for f in os.listdir('.') if f.startswith('gps')]

# parse data files and save into txt files
for file in files:
    # open output file, new name begins with gprmc
    outp = 'gprm'+file[4:]
    fgprmc = open(outp, 'w')
    print file
    fin = open(file)
    for line in fin:
        cleanLine = line.strip()
        lineOut = cleanLine + '\n'
        fields = cleanLine.split(',')
        if '$GPRMC' in fields[0]:
            fgprmc.write(lineOut)

    print outp+' finished'

fgprmc.close()


Test with: 
path_working_ddd = '/Users/eliz/DATA/ATOMIC/Brown/flux/Raw/20011';
ddd = 11;
yyyy = '2020';
PosLims = [-63   -45     8    20];

char = '0000078 $GPRMC,000000.00,A,5822.95924260,N,04515.92776422,W,1.58,223.96,151013,22.0,W,A*04';


%}

% reference timestamp
disp(['Reading gprm file for hour ',int2str(hhh)]);
start = (ddd+hhh/24);
delta = double(1.0/864000);
last = start + 35999*delta;
tref = start:delta:last;

%% read data if file exists, if not return NaNs
if exist(dfl,'file')==2
    %% read file
    flist = fopen(dfl);
    temp = {};    % empty arrays for textscan
    stx = [];
%     disp(dfl);
    
    while feof(flist)==0    % read entire file into cell array
        try
            
%             test = '5545733 $GPRMC,125545.00,A,1449.29376461,N,05103.11997654,W,0.77,195.96,150120,,,A*4E'
            % temp will be 1x11 cell array, cells will be ~36000x1 arrays
            temp = textscan(flist,'%2f%2f%3f %*6c %*9c %*c %f %c %f %c %f %f %f %*[^\n]','delimiter',', ','headerlines',0,'emptyvalue',NaN);
            temp{5} = double(temp{5}); % will be ascii code for N (78) or S (83)
            temp{7} = double(temp{7}); % will be ascii code for E (69) or W (87)
            stx = [stx cell2mat(temp)'];  % 11x~36000 array
%             disp('main')
        catch
            for ii=1:10      % length of last cell reflects missing values from any/all fields
                 if length(temp{1,ii})~=length(temp{1,10})
                     temp{1,ii}(length(temp{1,10})+1)=[];  %truncate length of all cells to length(temp{1,10})
%                      disp('catch 1')
                    fgetl(flist);
                 end
            end
            if ~isempty(temp{1,1})  % save good values from last read, if any
                temp{5} = double(temp{5});
                temp{7} = double(temp{7});
                stx = [stx cell2mat(temp)'];
%                 disp('catch 2')
            end
        end
    end
    
    disp('done with list');
    fclose(flist);
    jd_gps = ddd + (hhh + (stx(1,:)+(stx(2,:)+stx(3,:)/1000)/60)/60)/24;
    
    %% format variables
    glt = stx(4,:);
    gpslat_sign = stx(5,:);             % ascii code for "N" or "S"
    gln = stx(6,:);
    gpslon_sign = stx(7,:);             % ascii code for "W" or "E"
    gpsspeed = stx(8,:);                % in knots
    gpsspeed(gpsspeed>100) = NaN;       % kill wierd values
    gpshead = stx(9,:);                 % in deg
    gpshead(gpshead>361) = NaN;         % kill wierd values
    gpshead = unwrap(gpshead*pi/180);   % unwrapped in radians for now...

    bb = floor(glt/100);                % convert lat format to decimal degrees
    glt = (glt-bb*100)/60;
    glt = glt+bb;
    ii = find(gpslat_sign==83); glt(ii) = -glt(ii);  % s hemisphere = negative values

    bbn = floor(gln/100);               % convert lon format to decimal degrees
    gln = (gln-bbn*100)/60;
    gln = gln+bbn;
    ii = find(gpslon_sign==87);         % W longitudes
%     gln(ii) = -gln(ii) + 360;           % reformat to 0-360 deg *** choose based on experiment needs
    gln(ii) = -gln(ii);                 % reformat to +/- 180 deg *** choose based on experiment needs

    % discard unreasonable values
    gln(gln<PosLims(1)) = NaN;
    gln(gln>PosLims(2)) = NaN;
    glt(glt<PosLims(3)) = NaN;
    glt(glt>PosLims(4)) = NaN;
   
    
    %% despike variables
   disp('starting despike')
    [glt,~] = despike2(glt);    % despike & replace NaNs
    [gln,~] = despike2(gln);

    % interp
    disp('starting interp')
    [~,zz,~] = unique(jd_gps);
    gpslon1 = interp1(jd_gps(zz),gln(zz),tref,'linear','extrap');
    gpslat1 = interp1(jd_gps(zz),glt(zz),tref,'linear','extrap');
    gpshead1 = interp1(jd_gps(zz),gpshead(zz),tref,'linear','extrap');
    gpsspeed1 = interp1(jd_gps(zz),gpsspeed(zz),tref,'linear','extrap');   % 10 Hz sog

    %% look for lat/lon spikes again
    ii = abs(diff(gpslon1))>0.01; gpslon1(ii) = NaN;
    ii = abs(diff(gpslat1))>0.01; gpslat1(ii) = NaN;

    % Replace NaNs with median of +/- 20 point range around NaN
    ii=find(isnan(gpslon1));    % indices of NaNs
    if ~isempty(ii)
        jj = (ii<21);       % NaNs in the first 20 points
        kk = (ii>36000-21); % NaNs in the last 20 points
        ll = ~xor(jj,kk);   % all others
        gpslon1(ii(jj)) = nanmedian1(gpslon1(1:ii(jj)+20)); % median of 1st 20 pnts only
        gpslon1(ii(ll)) = nanmedian1(gpslon1(ii(ll)-20:ii(ll)+20));
        gpslon1(ii(kk)) = nanmedian1(gpslon1(ii(kk)-20:end)); % median of last 20 pnts only
    end

    ii=find(isnan(gpslat1));    % same for lat
    if ~isempty(ii)
        jj = (ii<21);
        kk = (ii>36000-21);
        ll = ~xor(jj,kk);
        gpslat1(ii(jj)) = nanmedian1(gpslat1(1:ii(jj)+20));
        gpslat1(ii(ll)) = nanmedian1(gpslat1(ii(ll)-20:ii(ll)+20));
        gpslat1(ii(kk)) = nanmedian1(gpslat1(ii(kk)-20:end));
    end

    ii=find(isnan(gpsspeed1));    % same for speed
    if ~isempty(ii)
        jj = (ii<21);
        kk = (ii>36000-21);
        ll = ~xor(jj,kk);
        gpsspeed1(ii(jj)) = nanmedian1(gpsspeed1(1:ii(jj)+20));
        gpsspeed1(ii(ll)) = nanmedian1(gpsspeed1(ii(ll)-20:ii(ll)+20));
        gpsspeed1(ii(kk)) = nanmedian1(gpsspeed1(ii(kk)-20:end));
    end

    ii=find(isnan(gpshead1));    % same for heading
    if ~isempty(ii)
        jj = (ii<21);
        kk = (ii>36000-21);
        ll = ~xor(jj,kk);
        gpshead1(ii(jj)) = nanmedian1(gpshead1(1:ii(jj)+20));
        gpshead1(ii(ll)) = nanmedian1(gpshead1(ii(ll)-20:ii(ll)+20));
        gpshead1(ii(kk)) = nanmedian1(gpshead1(ii(kk)-20:end));
    end

    cog = gpshead1.*180/pi;    % back to degrees
    cog = mod(cog+360,360);    % be sure range is 0-360
    sog = gpsspeed1*0.5144;    % in m/s

    %% compute ship N & E velocity components in m/s
    sogN = sog.*cos(cog./180.*pi); % 1x36000 arrays
    sogE = sog.*sin(cog./180.*pi);

    %------filter ship velocities
    fs = 10;         % GPS output @ 10 Hz
    fc = 1/200;      % Cutoff frequency
    wp = 1/25/(fs/2);
    ws = 0.8*wp;
    % [n,wn]=buttord(wp,ws,10,25);
    [n,wn] = buttord(wp,ws,3,7);
    [bhi,ahi] = butter(n,wn,'low');
    sogN = filtfilt(bhi,ahi,sogN);
    sogE = filtfilt(bhi,ahi,sogE);

    %% output
    disp('outputting')
    gps = [tref' cog' sog' sogN' sogE' gpslat1' gpslon1'];
else
    gps = NaN(36000,7);
    gps(:,1) = tref';
end

end
