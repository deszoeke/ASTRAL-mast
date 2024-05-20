function hed = read_hed_roll(dfl,ddd,hhh)
%{
reads heading/roll file from PSD Hemisphere system

outputs:
1   jd_ref      10 Hz timestamp
2   hed         heading, deg
3   roll        roll angle, deg

input parameters: dfl = string, full path to hed0 file
                  ddd = variable, day of year (julian date)
                  hh = variable, hour, 0-23
%}

% reference timestamp
disp(['Reading hed roll file for hour ',int2str(hhh)]);
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

    while feof(flist)==0	% read entire file into cell array
        try
            % temp will be 1x6 cell array, cells will be ~36000x1 arrays
            % 0000106 $PSAT,HPR,000002.00,184.19,,-1.73,N*3D
            temp = textscan(flist,'%2f%2f%3f %*5c %*3c %f %f %*f %f %*[^\n]','delimiter',', ','headerlines',1,'emptyvalue', NaN);
            stx = [stx cell2mat(temp)'];  % 6x~36000 array
        catch
            for ii=1:6      % length of last cell reflects missing values from any/all fields
                 if length(temp{1,ii})~=length(temp{1,6})
                     temp{1,ii}(length(temp{1,6})+1)=[];  %truncate length of all cells to length(temp{1,11})
                 end
            end
            if ~isempty(temp{1,1}) % save good values from last read, if any
                stx = [stx cell2mat(temp)'];
            end
        end
    end
    fclose(flist);
    jd_hed = ddd + (hhh + (stx(1,:)+(stx(2,:)+stx(3,:)/1000)/60)/60)/24;
    heading = replace_NaN_nearest_neighbor(stx(5,:));
    roll = replace_NaN_nearest_neighbor(stx(6,:));

    %% format to 3x36000 output array
    hed = NaN(36000,3);
    hed(:,1) = tref;

    % interp to 10 Hz grid
    hed(:,2) = interp1(jd_hed,heading,tref,'linear','extrap');
    hed(:,3) = interp1(jd_hed,roll,tref,'linear','extrap');

else
    hed = NaN(36000,3);
    hed(:,1) = tref';
end

end

