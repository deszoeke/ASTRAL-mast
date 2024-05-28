"general PSL DAS readers"
module NoaaDas

using Dates
using Printf

export get_das_filenames, get_das_pathfiles
export yday, psldatetime, decimal_hour
export read_das_dict

# utility functions
# m2n(x) = ismissing(x) ? NaN : x
# pd = permutedims

# general PSL DAS readers

# parameters for ASTRAL - feed defaults to functions
baseyear = 2024
mastdir = "./data/PSL/"
# ncolumn = 22 # change for different file types

"yearday from Date or DateTime"
yday(dt) = Dates.value(Date(dt) - Date(year(dt)-1,12,31))
yday(ddd::AbstractString) = parse(Int, ddd)

"return all filenames with prefix and date dt"
function get_das_filenames(prefix, yd::Integer, mastdir=mastdir)
    filter(startswith(prefix), 
           readdir(joinpath(mastdir, @sprintf("%03d",yd))))
end
function get_das_filenames(prefix, dt::Date, mastdir=mastdir)
    get_das_filenames(prefix, yday(dt), mastdir)
end
"return filenames with prefix matching hour of dt"
function get_das_filenames(prefix, dt::DateTime, mastdir=mastdir)
    hh = @sprintf("%02d", Dates.value(Hour(dt)))
    filter( endswith(hh*"_raw.txt"),
            get_das_filenames(prefix, Date(dt), mastdir) )
end

"return vector of full paths of all files with prefix and date"
function get_das_pathfiles(prefix, dt, mastdir=mastdir)
    joinpath.(mastdir,
              @sprintf("%03d",yday(dt)),
              get_das_filenames(prefix, dt, mastdir))
end

"Convert a PSL date to a DateTime."
function psldatetime(d::Date, hr::Integer, psltime::AbstractString)
    mm  = Minute(     parse(Int32, psltime[1:2]))
    SS  = Second(     parse(Int32, psltime[3:4]))
    sss = Millisecond(parse(Int32, psltime[5:7]))
    DateTime(Date(d)) + Hour(hr) + mm + SS + sss
end
function psldatetime(d::Date, hr::Integer=0, minute::Real=0)
    DateTime(Date(d)) + Hour(hr) + Minute(minute) + Second(secnd) + Millisecond(secnd)
end
function psldatetime(d::Date, hr::Integer=0, minute::Integer=0, secnd::Integer=0, millisecnd::Integer=0)
    DateTime(Date(d)) + Hour(hr) + Minute(minute) + Second(secnd) + Millisecond(millisecnd)
end
function psldatetime(yearday::Integer, hr::Integer=0; minute::Real=0, baseyear=baseyear)
    DateTime(baseyear-1,12,31) + Day(yearday) + Hour(hr) + Minute(minute)
end
# psldatetime(Date(2024,5,8), hr, "0606111")
# To get filename from the date, hour
# psldatefile(prefix, psldatetime(d, hr, psltime))

# file format default parameters
nheader = 1
# dt = DateTime(2024,5,8,11,0,0)
# pathfilename = get_das_pathfiles(prefix, dt::DateTime, mastdir=mastdir)

"read das data and put in a dictionary"
function read_das_dict(pathfilename, keys;
    nheader=1,
    nsample=sum(countlines.(pathfilename) .- nheader),
    ncolumn=length(keys)-1 ) # default columns in data excluding datestamp

    # read data from file
    datatime, X = read_das_data(pathfilename; 
        nheader=nheader,
        nsample=nsample,
        ncolumn=ncolumn ) # specifying manages memory best
    
    return das_dict(keys, datatime, X) # returns dict
end

parseblank2missing(T, s) = isempty(s) ? missing : parse(T, s)

"read and parse one das data file"
function read_das_data(pathfilename::AbstractString;
    nheader=1,
    nsample=countlines(pathfilenames) - nheader, 
    ncolumn=30 )

    # trivially iterate over 1-vector
    dt, X = read_das_data(pathfilename[1:1];  # changed to read_das_data
        nheader=nheader, nsample=nsample, ncolumn=ncolumn )
    return dt, X
end

"read and concatenate data from multiple files"
function read_das_data(pathfilename::Vector{<:AbstractString};
    nheader=1,
    nsample=sum( countlines.(pathfilenames) .- nheader ),
    ncolumn=30 ) # data, not including timestamp

    # preallocate the data
    # psltime = Vector{String}(undef, nsample) # will point to data as it is read
    psldt = Vector{DateTime}(undef, nsample) # will point to data as it is read
    X = Array{Union{Float32,Missing}, 2}(undef, nsample, ncolumn)
    fill!(X, missing)

    nl = 0
    maxcol=0
    for pfile in pathfilename
        # find hour from the filename
        shortfilename = last(splitpath(pfile))
        ddd = shortfilename[end-12:end-10]
        basedt = Date(baseyear-1,12,31) + Day(yday(ddd))
        hr = parse(Int32, shortfilename[end-9:end-8])

        open(pfile) do file
            for _ in 1:nheader
                readline(file) # skip header
            end
            for line in readlines(file)
                nl += 1
                splt = split(line, r"[\s,]+")

                nx = min(ncolumn, length(splt[2:end]))
                if nx > 0 # skip empty lines
                    psltime = splt[1]
                    psldt[nl] = psldatetime(basedt, hr, psltime)

                    dataline = try
                        parseblank2missing.(Float32, splt[2:end])
                    catch
                        error("failed to parse: $(splt[2:end])")
                    end
                    maxcol = max(maxcol, nx) # data in longest line
                    X[nl, 1:nx] .= dataline[1:nx]
                end
            end
        end
    end

    return psldt[1:nl], X[1:nl, 1:maxcol]
end

"assign data to a dict of keys"
function das_dict(keys, datatime, X)
    D = Dict{eltype(keys), Any}()

    # special DateTime
    D[keys[1]] = datatime 
    # fill rest of dictionary
    for (ik, ky) in enumerate(keys[2:end])
        D[ky] = X[:, ik]
    end
    return D
end

"compute the decimal hour of day from a datetime"
function decimal_hour(dt::DateTime)::Float64
    hour = Dates.hour(dt)
    minute = Dates.minute(dt)
    second = Dates.second(dt)
    millisecond = Dates.millisecond(dt)
    
    # Calculate the decimal hour
    decimal_hour = hour + minute / 60 + second / 3600 + millisecond / 3600000
    return decimal_hour
end

end

# GPS and heading-roll readers
module DasGps

using Dates
using ..NoaaDas

export read_gps_dict, read_hed_dict

# parameters for ASTRAL - feed defaults to functions
baseyear = 2024
mastdir = "./data/PSL/"

"yearday from Date or DateTime"
yday(dt) = Dates.value(Date(dt) - Date(year(dt)-1,12,31))
yday(ddd::AbstractString) = parse(Int, ddd)

"read gps file GPRMC data into a dict"
function read_gps_dict(pathfilenames;
    nheader=1,
    nsample=sum(countlines.(pathfilenames) .- nheader),
    ncolumn=7 )

    # method: psldatetime(d::Date, hr::Integer, psltime::String)

    # doesn't use das_dict() but it could.

    # preallocate the data
    psldt  = Vector{DateTime}(undef, nsample) # will point to data as it is read
    gpsdt  = Vector{DateTime}(undef, nsample) # will point to data as it is read
    lat    = Vector{Float64}(undef, nsample)
    lon    = Vector{Float64}(undef, nsample)
    speed  = Vector{Float64}(undef, nsample)
    course = Vector{Float64}(undef, nsample)

    nl = 0
    # handle a vector of pathfilenames or a single pathfilename
    p = pathfilenames isa Vector ? pathfilenames : Ref(pathfilenames)
    for pathfilename in p
        # find datetime from the filename
        shortfilename = last(splitpath(pathfilename))
        ddd = shortfilename[end-12:end-10]
        basedt = Date(baseyear-1,12,31) + Day(yday(ddd))
        hr = parse(Int32, shortfilename[end-9:end-8])

        open(pathfilename) do file
            for _ in 1:nheader
                readline(file) # skip header
            end
            for line in readlines(file)
                fields = split(line, [' ','\t'])
                if startswith(fields[2], "\$GPRMC")
                    nl += 1
                    psltime = fields[1] # string
                    psldt[nl] = psldatetime(Date(basedt), hr, psltime)

                    gpsdt[nl], lat[nl], lon[nl], speed[nl], course[nl] = parse_gprmc(fields[2])
                end
            end
        end
    end

    # Return a dictionary with the parsed fields
    D = Dict{Symbol, Any}()
    # DateTime
    D[:time]     = psldt[1:nl]
    D[:gpstime]  = gpsdt[1:nl]
    # Float
    D[:lat]      = lat[1:nl]
    D[:lon]      = lon[1:nl]
    D[:sog_kts]  = speed[1:nl]
    D[:cog_deg]  = course[1:nl]

    return D
end

"""
Parse one GPRMC NMEA message.

Arguments:
- `message::String`: GPRMC NMEA message

Returns:
- A dictionary containing the parsed fields
"""
function parse_gprmc(message::AbstractString)
    # Split the message into fields
    fields = split(message, ",")
    
    # Check that the message is a GPRMC message
    # and position is valid.
    if fields[1] != "\$GPRMC"
        error("Not a GPRMC message")
    elseif fields[3] != "A"
        error("Invalid position status")
    end
    
    # Extract the fields
    time_str = fields[2]
    # status = fields[3]
    lat_str = fields[4]
    lat_dir = fields[5]
    lon_str = fields[6]
    lon_dir = fields[7]
    speed_str = fields[8]
    course_str = fields[9]
    date_str = fields[10]
    
    # Parse the time (hhmmss.sss) and date (ddmmyy)
    time = tryparse(Time, time_str[1:6], dateformat"HHMMSS")
    date = tryparse(Date, date_str, dateformat"ddmmyy")

    # Combine date and time into a DateTime object
    gpstime = try 
        DateTime(date, time) 
    catch 
        @warn "Invalid date or time format"
    end

    # Parse latitude and longitude in degrees and minutes (ddmm.mmmm)
    lat_deg = tryparse(Float64, lat_str[1:2])
    lat_min = tryparse(Float64, lat_str[3:end])
    lat = lat_deg + lat_min / 60
    if lat_dir == "S"
        lat = -lat
    end

    lon_deg = tryparse(Float64, lon_str[1:3])
    lon_min = tryparse(Float64, lon_str[4:end])
    lon = lon_deg + lon_min / 60
    if lon_dir == "W"
        lon = -lon
    end

    # Parse speed (knots) and course (degrees)
    speed = tryparse(Float64, speed_str)
    course = tryparse(Float64, course_str)

    # Return a dictionary with the parsed fields
    return gpstime, lat, lon, speed, course
end
# Example usage
# message = "\$GPRMC,123519,A,4807.038,N,01131.000,E,022.4,084.4,230394,003.1,W*6A"
# parsed_data = parse_gprmc(message)
# println(parsed_data)

# "read and parse one hed data file"
# function read_hed_data(pathfilename::AbstractString;
#     nheader=1,
#     nsample=countlines(pathfilenames) - nheader, 
#     ncolumn=30 )

#     # trivially iterate over 1-vector
#     dt, X = read_hed_data(pathfilename[1:1]; 
#         nheader=nheader, nsample=nsample, ncolumn=ncolumn )
#     return dt, X
# end

"read differential gps heading data into a dict"
# reads into array _and_ puts in dictionary
function read_hed_dict(pathfilenames;
    nheader=1,
    nsample=sum(countlines.(pathfilenames) .- nheader),
    ncolumn=2 )

    # preallocate the data
    psldt   = Vector{DateTime}(undef, nsample) # will point to data as it is read
    heading = Vector{Float32}(undef, nsample)
    roll    = Vector{Float32}(undef, nsample)

    nl = 0
    # handle scalars or vectors with 1 line!
    p = pathfilenames isa Vector ? pathfilenames : Ref(pathfilenames)
    for pathfilename in p
        # find datetime from the filename
        shortfilename = last(splitpath(pathfilename))
        ddd = shortfilename[end-12:end-10]
        basedt = Date(baseyear-1,12,31) + Day(yday(ddd))
        hr = parse(Int32, shortfilename[end-9:end-8])
        # method: psldatetime(d::Date, hr::Integer, psltime::String)

        open(pathfilename) do file
            for _ in 1:nheader
                readline(file) # skip header
            end
            for line in readlines(file)
                fields = split(line, [' ','\t'])
                if startswith(fields[2], "\$PSAT,HPR")
                    nl += 1
                    psltime = fields[1] # string
                    psldt[nl] = psldatetime(Date(basedt), hr, psltime)
                    heading[nl], roll[nl] = parse_hpr(fields[2])
                end
            end
        end
    end

    # Return a dictionary with the parsed fields
    D = Dict{Symbol, Any}()
    D[:time]     = psldt[1:nl]    # DateTime
    D[:heading]  = heading[1:nl]    # Float
    D[:roll]     = roll[1:nl]
    return D
end

"""
Parse a PSAT,HPR NMEA message of the form
\$PSAT,HPR,000002.00,184.19,,-1.73,N*3D

Arguments:  message::String`: PSAT,HPR NMEA message
Returns:    heading::Float32, roll::Float:32
"""
function parse_hpr(message::AbstractString)
    # example
    # $PSAT,HPR,000002.00,184.19,,-1.73,N*3D

    # Split the message into fields
    fields = split(message, ",")
    
    # Check that the message is a PSAT message
    if !startswith(fields[1], "\$PSAT") || !startswith(fields[2], "HPR")
        error("Not a $PSAT,HPR message")
    end

    # Parse heading and roll
    heading_deg = tryparse(Float32, fields[4])
    roll_deg    = tryparse(Float32, fields[6])
    
    # Return a dictionary with the parsed fields
    return heading_deg, roll_deg
end

end # DasGps

module DasScs

using Dates
using ..NoaaDas
using ..NoaaDas: das_dict

export read_scs_dict

baseyear = 2024

# SCS readers process lat,lon into decimal degrees

"read das data and put in a dictionary"
function read_scs_dict(pathfilename, keys;
    nheader=1,
    nsample=sum(countlines.(pathfilename) .- nheader),
    ncolumn=length(keys)-1 ) # default columns in data excluding datestamp

    # read data from file
    datatime, X = read_scs_data(pathfilename; 
        nheader=nheader, nsample=nsample, ncolumn=ncolumn )
    return das_dict(keys, datatime, X) # returns dict
end

# parse gps lat,lon to decimal degrees
NSEWsgn(s) = endswith(uppercase(s),r"S|W") ? -1 : 1
declat(s) = NSEWsgn(s)*( parse(Float32, s[1:2]) + parse(Float32, s[3:end-1])/60 )
declon(s) = NSEWsgn(s)*( parse(Float32, s[1:3]) + parse(Float32, s[4:end-1])/60 )

parseblank2missing(T, s) = isempty(s) ? missing : parse(T, s)

"read and parse one file"
function read_scs_data(pathfilename::AbstractString;
    nheader=1, 
    nsample=countlines(pathfilenames) - nheader, 
    ncolumn=26 )

    # trivially iterate over 1-vector
    dt, X = read_scs_data(pathfilename[1:1]; 
        nheader=nheader, nsample=nsample, ncolumn=ncolumn )
    return dt, X
end

"read and concatenate data from multiple files"
function read_scs_data(pathfilename::Vector{<:AbstractString};
    nheader=1,
    nsample=sum( countlines.(pathfilenames) .- nheader ),
    ncolumn=26 ) # data, not including timestamp

    # preallocate the data
    # psltime = Vector{String}(undef, nsample) # will point to data as it is read
    psldt = Vector{DateTime}(undef, nsample) # will point to data as it is read
    X = Array{Union{Float32,Missing}, 2}(undef, nsample, ncolumn)
    fill!(X, missing)

    nl = 0
    maxcol=0
    for pfile in pathfilename
        # find hour from the filename
        shortfilename = last(splitpath(pfile))
        ddd = shortfilename[end-12:end-10]
        basedt = Date(baseyear-1,12,31) + Day(yday(ddd))
        hr = parse(Int32, shortfilename[end-9:end-8])

        open(pfile) do file
            for _ in 1:nheader
                readline(file) # skip header
            end
            for line in readlines(file)
                nl += 1
                splt = split(line, r"[\s,]+")

                nx = min(ncolumn, length(splt[2:end]))
                if nx > 0 # skip empty lines
                    psltime = splt[1]
                    psldt[nl] = psldatetime(basedt, hr, psltime)

                    lat = declat( splt[2] )
                    lon = declon( splt[3] )

                    dataline = parseblank2missing.(Float32, splt[4:end])
                    maxcol = max(maxcol, nx) # data in longest line

                    X[nl, 1:nx] .= cat(lat,lon,dataline[1:nx-2], dims=1)
                end
            end
        end
    end

    return psldt[1:nl], X[1:nl, 1:maxcol]
end

end


"read POSMV data, esp. PASHR pitch, roll, heave messages"
module ShipPosmv

using Dates

export get_posmv_file
export read_pashr_dict

function get_posmv_file(msg, dt; path="./data/scs/NAV") 
    # POSMV-V5-PASHR-RAW_20240510-000001.Raw
    ds = Dates.format(dt, dateformat"yyyymmdd")
    joinpath.( path, 
        filter(s -> startswith(s,"POSMV-V5-$(uppercase(msg))-RAW_$ds") 
                    && endswith(s,".Raw"), 
            readdir(path)) )
end

"expand scalars to make iterable, keep vectors as they are"
itr_expand(x) = x isa Vector ? x : Ref(x)

"read PSDMV-PASHR into a vector of columns X"
function read_pashr_data(fullfiles; nheader=0, nsample=sum(countlines.(fullfiles) .- nheader))

    # column designations for POSMV PASHR file
    colnames = split("date time nmeastring gpstime heading trueheading roll pitch heave roll_accuracy pitch_accuracy heading_accuracy gps_update_qualiy_flag ins_Status_flag checksum")
    coltypes = [Date; Time; String; Time; Float32; Bool; fill(Float32,6); UInt8; UInt8; String]
    colparsers = [  s -> Date(s, dateformat"mm/dd/yyyy");
                    s -> Time(s, dateformat"HH:MM:SS.sss");
                    s -> s;
                    s -> Time(s, dateformat"HHMMSS.sss");
                    s -> parse(Float32, s)
                    s -> s == "T"
                    fill(s -> parse(Float32, s), 6);
                    s -> parse(UInt8, s)
                    s -> parse(UInt8, s)
                    s -> s ]
    ncol = length(colnames)

    # initialize vectors in X to receive data
    X = Vector{Any}(undef, ncol)
    for (ci, ct) in enumerate(coltypes)
        X[ci] = Vector{ct}(undef, nsample)
    end
    # index as X[col][line]

    # read the file
    nl = 0
    for file in itr_expand(fullfiles)
        open(file) do f
            for line in readlines(f)
                s = split(line, r",|\*")
                if s[3] == "\$PASHR"
                    nl += 1
                    # parse and write the data to X
                    for (ci, cp) in enumerate(colparsers)
                        X[ci][nl] = cp(s[ci])
                    end
                end
            end
        end
    end
    return X, nl
end

function read_pashr_dict(fullfiles; nheader=0, nsample=sum(countlines.(fullfiles) .- nheader))
    # read the data into the omnivector
    X, nl = read_pashr_data(fullfiles; nheader=nheader, nsample=nsample)
    #       1    2    3          4       5       6           7    8     9     10            11             12               13                     14              15
    #       date time nmeastring gpstime heading trueheading roll pitch heave roll_accuracy pitch_accuracy heading_accuracy gps_update_qualiy_flag ins_Status_flag checksum
    #            1               2       3       4           5    6     7     8             9              10               11                     12              13            
    
    pashrkeys = Symbol.(
    split("     time            gpstime heading trueheading roll pitch heave roll_accuracy pitch_accuracy heading_accuracy gps_update_qualiy_flag ins_Status_flag checksum"))
   
    D = Dict{Symbol, Any}()

    D[:time] = X[1][1:nl] .+ X[2][1:nl] # make a datetime
    # assign the data columns to X, truncating to number of lines read
    for (i, k) in enumerate(pashrkeys[2:end])
        D[k] = X[i+3][1:nl]
    end
    return D
end

end


"read fast pressure sensors"
module DasFps

using Dates
using ..NoaaDas
export read_fps_dict

baseyear = 2024

"read fps das data and put in a dictionary"
function read_fps_dict(pathfilename, keys;
    nheader=1,
    nsample=sum(countlines.(pathfilename) .- nheader),
    ncolumn=length(keys)-1 ) # default columns in data excluding datestamp

    # read data from file
    datatime, X = read_fps_data(pathfilename; 
        nheader=nheader, nsample=nsample, ncolumn=ncolumn )
    return NoaaDas.das_dict(keys, datatime, X) # returns dict
end

"read and parse one file"
function read_fps_data(pathfilename::AbstractString;
    nheader=1, 
    nsample=countlines(pathfilenames) - nheader, 
    ncolumn=26 )

    # trivially iterate over 1-vector
    dt, X = read_fps_data(pathfilename[1:1]; 
        nheader=nheader, nsample=nsample, ncolumn=ncolumn )
    return dt, X
end

parseblank2missing(T, s) = isempty(s) ? missing : parse(T, s)

"read and concatenate data from multiple files"
function read_fps_data(pathfilename::Vector{<:AbstractString};
    nheader=1,
    nsample=sum( countlines.(pathfilenames) .- nheader ),
    ncolumn=1 ) # data, not including timestamp

    # preallocate the data
    # psltime = Vector{String}(undef, nsample) # will point to data as it is read
    psldt = Vector{DateTime}(undef, nsample) # will point to data as it is read
    X = Array{Union{Float32,Missing}, 2}(undef, nsample, ncolumn)
    fill!(X, missing)

    nl = 0
    maxcol=0
    for pfile in pathfilename
        # find hour from the filename
        shortfilename = last(splitpath(pfile))
        ddd = shortfilename[end-12:end-10]
        basedt = Date(baseyear-1,12,31) + Day(yday(ddd))
        hr = parse(Int32, shortfilename[end-9:end-8])

        open(pfile) do file
            for _ in 1:nheader
                readline(file) # skip header
            end
            for line in readlines(file)
                nl += 1
                splt = split(line, r"[\s\*]+") # as many delims as possible, counting literal *

                nx = min(ncolumn, length(splt[2:end]))
                if nx > 0 # skip empty lines
                    psltime = splt[1]
                    psldt[nl] = psldatetime(basedt, hr, psltime)

                    dataline = parseblank2missing.(Float32, splt[2:end])
                    maxcol = max(maxcol, nx) # data in longest line

                    X[nl, 1:nx] .= dataline[1:nx]
                end
            end
        end
    end

    return psldt[1:nl], X[1:nl, 1:maxcol]
end

end