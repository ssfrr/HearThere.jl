module HearThere

# currently using my hacked version of OSC.jl to enable 1-based indexing and
# the path method. PR has been filed
using OSC
using DataFrames
using PyPlot
using Quaternions
using LightXML
using DataFramesMeta
using Compat

export getrawdata, getcookeddata


# full-scale range in g's (based on 16-bit signed ADC word)
const accelFSR = 16
const accelScale = accelFSR / (2^15)
# full-scale range in degrees / s (based on 16-bit signed ADC word)
const gyroFSR = 2000
const gyroScale = gyroFSR / (2^15)
# full-scale range in μT (based on 14-bit signed ADC word)
const magFSR = 4800
const magScale = magFSR / (2^13)

# magnetic field strength in Cambridge in μT
const magField = 52

"""
Listen on the given port for OSC messages of the form:

    /rawimu iiiiiiiii <accel_x> <accel_y> <accel_z> <gyro_x> <gyro_y> <gyro_z> <mag_x> <mag_y> <mag_z>

and return a DataFrame of the timestamped values.
"""
function getrawdata(duration::Real, oscport::Integer=10001)
    timestamps = Float64[]
    data = [Int32[] for _ in 1:9]
    oscpath = "/rawimu"
    sock = UdpSocket()
    if !bind(sock, ip"0.0.0.0", oscport)
        println("Couldn't bind to port $oscport")
        return makerawdf(timestamps, data)
    end
    try
        # wait for the first message
        println("Waiting for first message...")
        while true
            msg = OscMsg(recv(sock))
            if path(msg) == oscpath
                break
            end
        end
        println("Capturing Data...")
        starttime = time()
        now = time()
        lastreported = 0
        while(now < starttime + duration)
            if(now - lastreported >= 5)
                println("$(round(starttime + duration - now)) seconds left")
                lastreported = now
            end
            msg = OscMsg(recv(sock))
            now = time()
            if path(msg) == oscpath
                push!(timestamps, now)
                for i in 1:9
                    push!(data[i], msg[i])
                end
            end
        end
    finally
        close(sock)
    end
    makerawdf(timestamps, data)
end

"""
Listen on the given port for OSC messages of the form:

    /body_orientation ihffff id timestamp q1 q2 q3 q0
    /body_location  ihfff id timestamp x y z
and
    /ranges ffff anchor0 anchor1 anchor2 anchor3

and return a tuple of DataFrames of the timestamped values.
"""
function getoptidata(duration::Real, oscport::Integer=10001)
    orientationpath = "/body_orientation"
    locationpath = "/body_location"
    orientationdf = DataFrame(
        id=Int32[],
        timestamp=Float64[],
        q0=Float32[],
        q1=Float32[],
        q2=Float32[],
        q3=Float32[])
    locationdf = DataFrame(
        id=Int32[],
        timestamp=Float64[],
        x=Float32[],
        y=Float32[],
        z=Float32[])
    sock = UdpSocket()
    if !bind(sock, ip"0.0.0.0", oscport)
        println("Couldn't bind to port $oscport")
        return (locationdf, orientationdf)
    end
    try
        # wait for the first message
        println("Waiting for first message...")
        while true
            msg = OscMsg(recv(sock))
            if path(msg) in [orientationpath, locationpath]
                break
            else
                println("got unrecognized path: $(path(msg))")
            end
        end
        println("Capturing Data...")
        starttime = time()
        now = time()
        lastreported = 0
        while(now < starttime + duration)
            if(now - lastreported >= 5)
                println("$(round(starttime + duration - now)) seconds left")
                lastreported = now
            end
            msg = OscMsg(recv(sock))
            now = time()
            # opti-track data has an inverted z, so we need to convert
            if path(msg) == orientationpath
                push!(orientationdf, [msg[1], now, msg[6], -msg[3], -msg[4], msg[5]])
            elseif path(msg) == locationpath
                push!(locationdf, [msg[1], now, msg[3], msg[4], -msg[5]])
            else
                println("Got unrecognized path $(path(msg))")
            end
        end
    catch InterruptException
        println("Interrupted.")
    end
    close(sock)
    (locationdf, orientationdf)
end

"""
Listen on the given port for OSC messages of the form:

    /orientation    ffff    q0 q1 q2 q3
    /ranges         ffff     anchor0 anchor1 anchor2 anchor3
    /gpsgeo         ddddd   latitude longitude acceleration horiz_accuracy, vert_accuracy
    /gpsxyz         fff     x y z
    /uwblocalxyz    ffffff  x y z std_x std_y std_z
    /uwbglobalxyz   ffff    x y z confidence
    /fusedxyz       fff     x y z

and return a dictionary of DataFrames of the timestamped values.
"""
function getcookeddata(duration::Real, fileprefix::AbstractString, oscport::Integer=10001)
    dfs = @compat Dict{UTF8String, Any}(
        "/orientation" => DataFrame(
            timestamp=Float64[],
            q0=Float32[], q1=Float32[], q2=Float32[], q3=Float32[]),
        "/ranges" => DataFrame(
            timestamp=Float64[],
            anchor0=Float32[], anchor1=Float32[], anchor2=Float32[], anchor3=Float32[]),
        "/gpsgeo" => DataFrame(
            timestamp=Float64[],
            latitude=Float64[], longitude=Float64[], elevation=Float64[],
            horizAccuracy=Float32[],
            vertAccuracy=Float32[]),
        "/gpsxyz" => DataFrame(
            timestamp=Float64[],
            x=Float32[], y=Float32[], z=Float32[]),
        "/uwblocalxyz" => DataFrame(
            timestamp=Float64[],
            x=Float32[], y=Float32[], z=Float32[],
            stdX=Float32[], stdY=Float32[], stdZ=Float32[]),
        "/uwbglobalxyz" => DataFrame(
            timestamp=Float64[],
            x=Float32[], y=Float32[], z=Float32[],
            confidence=Float32[]),
        "/fusedxyz" => DataFrame(
            timestamp=Float64[],
            x=Float32[], y=Float32[], z=Float32[])
        )

    # returns the filled-in dataframes
    collect_data(dfs, duration, oscport)

    for metric in keys(dfs)
        writetable("$(fileprefix)_$(metric[2:end]).tsv", dfs[metric])
    end

    dfs
end

function collect_data{T<:AbstractString}(dfs::Dict{T, DataFrame}, duration, oscport)
    sock = UdpSocket()
    if !bind(sock, ip"0.0.0.0", oscport)
        println("Couldn't bind to port $oscport")
        return dfs
    end
    try
        # wait for the first message
        println("Waiting for first message...")
        while true
            msg = OscMsg(recv(sock))
            if path(msg) in keys(dfs)
                break
            else
                println("Got unrecognized message path \"$(path(msg))\"")
            end
        end
        println("Capturing Data...")
        starttime = time()
        now = time()
        lastreported = 0
        while(now < starttime + duration)
            if(now - lastreported >= 5)
                println("$(round(starttime + duration - now)) seconds left")
                lastreported = now
            end
            msg = OscMsg(recv(sock))
            now = time()
            msgPath = path(msg)
            df = DataFrame()
            try
                df = dfs[msgPath]
            catch e
                isa(e, KeyError) || rethrow(e)
                println("Got unrecognized message path \"$(path(msg))\"")
                continue
            end

            args = [msg[i] for i in 1:(size(df, 2)-1)]
            row = [now, args...]
            push!(df, [now, args...])
        end
    catch e
        isa(e, InterruptException) || rethrow(e)
        println("Interrupted.")
    finally
        close(sock)
    end

    dfs
end

function rotate(q::Quaternion, v::Vector)
    q = normalize(q)

    imag(q*Quaternion(v)*conj(q))
end

function plotpath(data)
    PyPlot.plot3D(convert(Array, data[:x], NaN),
        convert(Array, data[:y], NaN),
        convert(Array, data[:z], NaN),
        alpha=0.5)
    xlabel("x")
    ylabel("y")
    zlabel("z")
end

function plotorientation(data)
    orientations = Quaternion{Float32}[Quaternion(d...) for d in zip(data[:q0], data[:q1], data[:q2], data[:q3])]
    # this is probably really slow
    path = hcat([rotate(o, Float32[0, 0, 1]) for o in orientations]...)
    PyPlot.plot3D(vec(path[1, :]), vec(path[2, :]), vec(path[3, :]))
    xlim([-1, 1])
    ylim([-1, 1])
    zlim([-1, 1])
    xlabel("x")
    ylabel("y")
    zlabel("z")
end

"""
Takes a set of turntable data collections and returns the x, y, and z
scaling coefficients. Each collection should be with the relevant axis
pointing down, i.e. rotation should be positive. bias_data should be
a measurement when the device is not moving.
"""
function calibrate_gyro_scale(x33, x45, y33, y45, z33, z45, bias_data)
    # convert rpm into gyro ADC units
    units33 = 33 * 360 / 60 / gyroScale
    units45 = 45 * 360 / 60 / gyroScale
    xs = [0, units33, units45]

    # as in y = mx + b
    as = Float64[]
    bs= Float64[]
    x0mean = mean(bias_data[:gyro_x])
    x33mean = mean(x33[:gyro_x])
    x45mean = mean(x45[:gyro_x])
    a, b = linreg(xs, [x0mean, x33mean, x45mean])

    # linreg gives us a, b s.t. x ~= a + b*x', (minimizing the squared error)
    # where x is the measured value and x' is the actual value. Reversing that
    # gives us x' = (x - a) / b. We define our conversion values in terms of addition
    # and multiplication respectively, so we invert them here.
    push!(as, -a)
    push!(bs, 1/b)

    y0mean = mean(bias_data[:gyro_y])
    y33mean = mean(y33[:gyro_y])
    y45mean = mean(y45[:gyro_y])
    a, b = linreg(xs, [y0mean, y33mean, y45mean])
    push!(as, -a)
    push!(bs, 1/b)

    z0mean = mean(bias_data[:gyro_z])
    z33mean = mean(z33[:gyro_z])
    z45mean = mean(z45[:gyro_z])
    a, b = linreg(xs, [z0mean, z33mean, z45mean])
    push!(as, -a)
    push!(bs, 1/b)

    (as, bs)
end

"""
A very simplistic ellipsoid fit of the given data. It just uses the min/max values
of each axis to generate a containing ellipse. It also ignores any rotation.

It returns a tuple of 3-element vectors, the first is the offset α and the 2nd is the
scaling β.

They are meant to be applied as x′ = β(x+α), where x′ is the corrected value
and x is the raw value. We apply the offset first so it can be done in the
integer domain more efficiently.

"""
function calibrate_mag(data)
    xmax = maximum(data[:mag_x])
    xmin = minimum(data[:mag_x])
    ymax = maximum(data[:mag_y])
    ymin = minimum(data[:mag_y])
    zmax = maximum(data[:mag_z])
    zmin = minimum(data[:mag_z])

    α = [
        (xmin + xmax) / -2,
        (ymin + ymax) / -2,
        (zmin + zmax) / -2]
    β = [
         2magField / ((xmax - xmin) * magScale),
         2magField / ((ymax - ymin) * magScale),
         2magField / ((zmax - zmin) * magScale)]

    (α, β)
end

function correct_mag!(data, α, β)
    data[:mag_x] = (data[:mag_x] + α[1]) * β[1]
    data[:mag_y] = (data[:mag_y] + α[2]) * β[2]
    data[:mag_z] = (data[:mag_z] + α[3]) * β[3]
end

function plot_mag(data)
    scatter3D(data[:mag_x], data[:mag_y], data[:mag_z])
    xlabel("x")
    ylabel("y")
    zlabel("z")

    range = maximum(abs([
        minimum(data[:mag_x]),
        minimum(data[:mag_y]),
        minimum(data[:mag_z]),
        maximum(data[:mag_x]),
        maximum(data[:mag_y]),
        maximum(data[:mag_z])]))
    xlim([-range, range])
    ylim([-range, range])
    zlim([-range, range])
end

function plot_gyro(data, bias=[0.0, 0.0, 0.0])
    plot(data[:gyro_x] + bias[1], label="x")
    plot(data[:gyro_y] + bias[2], label="y")
    plot(data[:gyro_z] + bias[3], label="z")
end

"""
Make a dataframe from the timestamps and collected data
"""
function makerawdf(timestamps, data)
    DataFrame(
        timestamp=timestamps,
        accel_x=data[1],
        accel_y=data[2],
        accel_z=data[3],
        gyro_x=data[4],
        gyro_y=data[5],
        gyro_z=data[6],
        mag_x=data[7],
        mag_y=data[8],
        mag_z=data[9])
end

function parseoptitrack(filename)
    bodies = Dict{Int, AbstractString}()
    rootElement = root(parse_file(filename))
    for bodyDesc in get_elements_by_tagname(rootElement, "RigidBodyDesc")
        name = content(find_element(bodyDesc, "Name"))
        id = int(content(find_element(bodyDesc, "ID")))
        bodies[id] = name
    end

    df = DataFrame(
        timestamp=Float64[],
        id=Int[],
        x=Float64[],
        y=Float64[],
        z=Float64[],
        q0=Float64[],
        q1=Float64[],
        q2=Float64[],
        q3=Float64[]
    )
    for dfElement in get_elements_by_tagname(rootElement, "DataFrame")
        # the Latency field apparently means the capture computer's timestamp. See
        # https://forums.naturalpoint.com/viewtopic.php?p=57338
        timestamp = float(content(find_element(dfElement, "Latency")))
        bodydataElement = find_element(dfElement, "RigidBodyData")
        # the rigid body properties aren't contained inside an XML tag, they're all
        # just repeated, so we iterate through the zipped list
        # we're also switching qz and qy as described here:
        # https://forums.naturalpoint.com/viewtopic.php?p=42839
        for (id, x, y, z, q0, q1, q2, q3) in zip(
                get_elements_by_tagname(bodydataElement, "RigidBodyID"),
                get_elements_by_tagname(bodydataElement, "x"),
                get_elements_by_tagname(bodydataElement, "y"),
                get_elements_by_tagname(bodydataElement, "z"),
                get_elements_by_tagname(bodydataElement, "qw"),
                get_elements_by_tagname(bodydataElement, "qx"),
                get_elements_by_tagname(bodydataElement, "qz"),
                get_elements_by_tagname(bodydataElement, "qy"))
            # the optitrack data is 1/2 what it should be because of a calibration
            # error. Also the Z-axis is reversed
            push!(df, [
                timestamp,
                int(content(id)),
                float(content(x)) * 2,
                float(content(y)) * 2,
                -float(content(z)) * 2,
                float(content(q0)),
                float(content(q1)),
                float(content(q2)),
                -float(content(q3))])
        end
    end
    # The OptiTrack reports a location of 0, 0, 0 if the trackable isn't present
    @byrow! df begin
        if :x == 0.0 && :y == 0.0 && :z == 0.0
            :x = NA
            :y = NA
            :z = NA
            :q0 = NA
            :q1 = NA
            :q2 = NA
            :q3 = NA
        end
    end
    df
end

"""
Calculate the ranges for the given path (with x, y, z columns) and add them
to the dataframe
"""
function addranges!(optipath, anchors)
    for (anc, i) in [(:anchor0, 1), (:anchor1, 2), (:anchor2, 3), (:anchor3, 4)]
        optipath[anc] = sqrt(
            (anchors[:x][i] - optipath[:x]).^2 +
            (anchors[:y][i] - optipath[:y]).^2 +
            (anchors[:z][i] - optipath[:z]).^2)
    end
    optipath
end

"""
Takes 2 points (a1 and a2) and two ranges to those points (r1 and r2) and finds
the 0, 1 or 2 points that satisfy the ranges. The are returned as a 2xN matrix,
where N is the number of solutions.
"""
function trilaterate(a1, a2, r1, r2)
    # we first solve in a coordinate system where a1 is at the origin
    x2 = a2[1] - a1[1]
    y2 = a2[2] - a1[2]

    α = r1^2 - r2^2 + x2^2
    a = 1 + x2^2/y2^2
    b = -(α*x2+x2*y2^2)/y2^2
    c = (α^2 + 2α*y2^2+y2^4)/4y2^2 - r1^2

    det = b^2-4a*c
    if det < 0
        sol = Array(Float64, 2, 0)
    elseif det == 0.0
        x0 = -b/2a
        y0 = (r1^2-r2^2+x2^2-2x2*x0+y2^2)/2y2
        sol = [x0, y0]
    else
        # possibly 2 solutions to quadratic
        x0 = [(-b + √det)/2a, (-b - √det)/2a]
        y0 = (r1^2-r2^2+x2^2-2x2*x0+y2^2)/2y2
        sol = vcat(x0', y0')
    end
    # now re-apply the offset. Note this will broadcast in the 2x2 case
    sol .+ a1
end

function cleanranges(ranges)
    for anc in [:anchor0, :anchor1, :anchor2, :anchor3]
        ranges[ranges[anc] .== -1.0, anc] = NA
    end
    ranges
end

function getanchorlocations(rangefile, locationfile)
    anchorDF = DataFrame(
        anchor=0:3,
        x=zeros(4),
        y=zeros(4),
        z=zeros(4),
    )

    locdata = readtable(locationfile)
    counts = zeros(Int, 4)
    for (id, x, y, z) in zip(locdata[:anchor], locdata[:x], locdata[:y], locdata[:z])
        # the values are doubled because the optitrack was calibrated incorrectly
        # so all the reported values are half what they should be. We also
        # divide by 1000 to convert mm to m. Also the Z axis is reversed
        anchorDF[id+1, :x] += 2x / 1000
        anchorDF[id+1, :y] += 2y / 1000
        anchorDF[id+1, :z] += -2z / 1000
        counts[id+1] += 1
    end

    anchorDF[:x][3:4] = anchorDF[:x][3:4] ./ counts[3:4]
    anchorDF[:y][3:4] = anchorDF[:y][3:4] ./ counts[3:4]
    anchorDF[:z][3:4] = anchorDF[:z][3:4] ./ counts[3:4]

    ranges = readtable(rangefile)
    r20 = mean(dropna(ranges[:r20]))
    r21 = mean(dropna(ranges[:r21]))
    r30 = mean(dropna(ranges[:r30]))
    r31 = mean(dropna(ranges[:r31]))

    # we start out considering anchor 3 (index 4) to be the origin
    a0 = trilaterate(
        [anchorDF[4, :x], anchorDF[4, :z]],
        [anchorDF[3, :x], anchorDF[3, :z]],
        r30, r20)
    @assert size(a0) == (2, 2)
    a1 = trilaterate(
        [anchorDF[4, :x], anchorDF[4, :z]],
        [anchorDF[3, :x], anchorDF[3, :z]],
        r31, r21)
    @assert size(a1) == (2, 2)

    # we should get 2 solutions and we happend to know that the correct one is
    # in front of (more positive Z) a3 and a2
    for i in 1:2
        if a0[2, i] > anchorDF[4, :z]
            anchorDF[1, :x] = a0[1, i]
            anchorDF[1, :y] = anchorDF[4, :y]
            anchorDF[1, :z] = a0[2, i]
        end
        if a1[2, i] > anchorDF[3, :z]
            anchorDF[2, :x] = a1[1, i]
            anchorDF[2, :y] = anchorDF[3, :y]
            anchorDF[2, :z] = a1[2, i]
        end
    end

    anchorDF
end


"""
Takes a dataframe of (possibly irregularly) timestamped data and resamples to the given resolution as
if the data had been sampled at a regular interval. Returns an array
of the resampled data. This assumes the data is sorted by timestamp and timestamp.
is in seconds.
"""
function samplereg(data::DataFrame, resolution::AbstractFloat, timefield=:timestamp)
    colnames = filter(n -> n != timefield, names(data))
    endidx = length(data[timefield])
    out = DataFrame()
    out[timefield] = data[timefield][1]:resolution:data[timefield][endidx]
    n = length(out[timefield])
    for name in colnames
        out[name] = Array(Float64, n)
    end

    quat = false
    quatcols = [:q0, :q1, :q2, :q3]
    # we treat quaternions columns specially if present
    if all([col in colnames for col in quatcols])
        quat = true
        # remove the quaternions columns from the regular list
        colnames = filter(n -> !(n in quatcols), colnames)
        for col in quatcols
            out[col] = Array(Float64, n)
        end
    end

    # set up the current interval start/end times
    intstart = data[timefield][1]
    intend = data[timefield][2]

    idx = 1 # input index
    odx = 1 # output index
    while odx <= n
        sampletime = out[odx, timefield]
        # find the data timestamps that bracket our current sample time
        while sampletime > intend
            idx += 1
            intstart = intend
            intend = data[timefield][idx+1]
        end
        if intend == intstart
            α = 0.0
        else
            α = (sampletime - intstart) / (intend - intstart)
        end
        for name in colnames
            if isna(data[idx, name]) || isna(data[idx+1, name])
                out[odx, name] = NA
            else
                out[odx, name] = lerp(data[idx, name], data[idx+1, name], α)
            end
        end
        if quat
            if any([isna(data[idx+off, c]) for off in [0, 1], c in quatcols])
                for col in quatcols
                    out[odx, col] = NA
                end
            else
                p = Quaternion([data[idx, col] for col in quatcols]...)
                q = Quaternion([data[idx+1, col] for col in quatcols]...)
                lerped = nlerp(p, q, α)
                out[odx, :q0] = lerped.s
                out[odx, :q1] = lerped.v1
                out[odx, :q2] = lerped.v2
                out[odx, :q3] = lerped.v3
            end
        end
        odx += 1
    end
    out
end

Quaternion([1, 0, 0, 0]...)
# HearThere.lerp(10, 15, 1)
# p = Quaternion(1, 0, 0, 0.0)
# q = Quaternion(0, 1, 0, 0.0)
# HearThere.nlerp(p, q, 0.99)

"""
returns the lag (in samples) that vector a needs to be delayed relative to
vector b to maximize their cross-correlation, or alternatively the number of
sample of vector b that need to be removed to align the vectors. Assumes the
vectors are at the same sample rate.
"""
function getlag(a, b)
    xc = xcorr(a - mean(a), b - mean(b))
    _, i = findmax(xc)

    max(length(a), length(b)) - i
end

"""
Takes two DataFrames at the same sample rate and returns new ones that are
time-aligned and the same length. The intersection of the column names
excluding the time field are used for the cross-correlation if no offset
is given.
"""
function align(a, b; offset=nothing, timefield=:timestamp)
    if offset == nothing
        acolnames = filter(n -> n != timefield, names(a))
        bcolnames = filter(n -> n != timefield, names(b))
        colnames = intersect(acolnames, bcolnames)

        lagsum = 0
        for colname in colnames
            lagsum += getlag(convert(Array, a[colname], 0.0), convert(Array, b[colname], 0.0))
        end
        avglag = int(round(lagsum / length(colnames)))
    else
        avglag = offset
    end
    if avglag >= 0
        # b needs to be trimmed at the beginning
        bstart = avglag + 1
        astart = 1
    else
        println("WARNING: negative lag ($avglag) untested")
        # a needs to be trimmed at the beginning
        astart = -avglag + 1
        bstart = 1
    end

    n = min(size(a[astart:end, :], 1), size(b[bstart:end, :], 1))

    a = a[astart:astart+n-1, :]
    b = b[bstart:bstart+n-1, :]

    # we don't know what the correct times are anyways, so align them starting
    # at zero
    a[timefield] -= a[1, timefield]
    b[timefield] = a[timefield]
    (a, b)
end


function combine_range_calib{T <: Real, S <: AbstractString}(
    distances::Array{T},
    measurementfiles::Array{S},
    undobias=false)
    # these are the horizontal offsets of the anchors during calibration in meters
    offsets = @compat Dict(
        :anchor0 => -0.165,
        :anchor1 => -0.055,
        :anchor2 => 0.055,
        :anchor3 => 0.165
    )
    df = DataFrame(
        timestamp=Float64[],
        anchor=UTF8String[],
        actual=Float64[],
        measured=Float64[]
    )
    for (r, f) in zip(distances, measurementfiles)
        measuredf = readtable(f)
        for anc in [:anchor0, :anchor1, :anchor2, :anchor3]
            actual = sqrt(offsets[anc]^2 + r^2)
            anchordf = DataFrame(
                timestamp=measuredf[:timestamp],
                anchor=string(anc),
                actual=ones(size(measuredf, 1)) * actual,
                measured=measuredf[anc]
            )
            @byrow! anchordf begin
                if :measured == -1.0 || :measured > 1000
                    :measured = NA
                elseif undobias
                    # reverse the bias that the headtracker is doing internally
                    :measured = reversebias(:measured)
                end
            end
            df = vcat(df, anchordf)
        end
    end

    df
end


"""
Takes 2 dataframes and returns a new dataframe with the errors between
matching columns
"""
function geterr(a, b, timefield=:timestamp)
    acolnames = filter(n -> n != timefield, names(a))
    bcolnames = filter(n -> n != timefield, names(b))
    colnames = intersect(acolnames, bcolnames)

    errdf = DataFrame()
    for colname in colnames
        errdf[colname] = a[colname] - b[colname]
    end
    errdf
end

# this is ported over from the decawave code and is used to correct the range
# measurements
function getrangebias(range)
    # bias table from the decawave driver, used in the tag code to correct the distance
    biastable = [ 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 15, 18, 20, 22, 24, 27,
                  29, 32, 35, 38, 41, 44, 47, 51, 55, 58, 62, 66, 71, 78, 85,
                  96, 111, 135, 194, 240, 255 ]
    biasoff = -23

    range25 = range * 4 # in units of 25cm
    if range25 > 255
        range25 = 255
    end

    i = 1
    while range25 > biastable[i]
        i += 1
    end
    bias = i - 1 + biasoff # bias is in cm
    bias * 0.01 # convert to m
end

function adjustrange(range)
    range - getrangebias(range)
end

"""
Try to reverse the bias applied in the tag. This can have an error of 1cm
"""
function reversebias(range)
    range + getrangebias(range)
end

"""
Convert a quaternion to the Aerospace Euler Angle representation (Z-Y'-X'').
This representation considers X to be forward (towards the nose), Y is to the
right, and Z is down. This sequence is also known as Yaw-Pitch-Roll. Note that
the primes indicate that the rotations are about the object's local coordinate
system, which is transformed at each step. In the global coordinate system
the same rotation is achieved by rotating in the reverse order, X-Y-Z. Also
note that rotations follow the right-hand rule.

See Kuipers "Quaternions and Rotation Sequences" (1999) for more info.
"""
function aerospace(q::Quaternion)
    # TODO: handle gimbal lock situations
    q = normalize(q)
    # test = 2(q.v1*q.v3+q.s*q.v2)
    test = 2(-q.v1*q.v3+q.s*q.v2)
    if test > 0.998
        ϕ = 0;
        θ = pi/2
        ψ = 2 * atan2(q.v1, q.s);
    elseif test < -0.998
        ϕ = 0;
        θ = -pi/2
        ψ = -2 * atan2(q.v1, q.s);
    else
        ϕ = -atan2(2(q.s*q.v1 + q.v2*q.v3), 1 - 2(q.v1^2+q.v2^2))
        θ = asin(2(q.s*q.v2 - q.v3*q.v1))
        ψ = -atan2(2(q.s*q.v3 + q.v1*q.v2), 1 - 2(q.v2^2+q.v3^2))
    end

    # yaw, pitch, roll
    (ψ, θ, ϕ)
end

function addaerospace!(df::DataFrame)
    yaw = DataArray(Float64, 0)
    pitch = DataArray(Float64, 0)
    roll = DataArray(Float64, 0)
    for i in 1:size(df, 1)
        # convert to aerospace axis convention
        q0 = df[i, :q0]
        q1 = -df[i, :q3]
        q2 = -df[i, :q1]
        q3 = df[i, :q2]

        if isna(q0) || isna(q1) || isna(q2) || isna(q3)
            y = NA
            p = NA
            r = NA
        else
            y, p, r = aerospace(Quaternion(q0, q1, q2, q3, true))
        end
        push!(yaw, y)
        push!(pitch, p)
        push!(roll, r)
    end
    df[:yaw] = yaw
    df[:pitch] = pitch
    df[:roll] = roll

    nothing
end

"""
Linearly interpolate between two values.
"""
lerp(lhs, rhs, t::Real) = (1-t)*lhs + t*rhs

"""
Dot product of two quaternions, as if they're vectors in 4D
"""
dot(p::Quaternion, q::Quaternion) = p.s*q.s + p.v1*q.v1 + p.v2*q.v2 + p.v3*q.v3

"""
Interpolate between two quaternions by linearly interpolating their
components and ensuring it's still a valid rotation by normalizing. This is
still guaranteed to be on the great-circle path between the two quaternions,
(like slerp), but will not travel at a constant rotational velocity with respect
to t.
"""
function nlerp(p::Quaternion, q::Quaternion, t::Real)
    if dot(p, q) < 0
        # negating the quaternion doesn't change the rotation it represents
        # but this ensures we take the short way around instead of the long way
        q = -q
    end
    res = Quaternion(
        lerp(p.s, q.s, t),
        lerp(p.v1, q.v1, t),
        lerp(p.v2, q.v2, t),
        lerp(p.v3, q.v3, t))

    normalize(res)
end


"""
Calculates the error between two unit quaternions representing rotations.
Uses metric Φ6 from Huynh's "Metrics for 3D Rotations: Comparison and Analysis",
which has units of radians.
"""
err(p::Quaternion, q::Quaternion) = 2acos(abs(dot(p, q)))

"""
Calculates the quaternion error for the given dataframes. Assumes they're the
same length
"""
function err(df1::DataFrame, df2::DataFrame)
    @assert(size(df1, 1) == size(df2, 1))
    len = size(df1, 1)
    er = DataFrame(timestamp=Array(Float64, len), err=Array(Float64, len))
    quatcols = [:q0, :q1, :q2, :q3]
    for i in 1:len
        er[i, :timestamp] = df1[i, :timestamp]
        if any(Bool[isna(df[i, col]) for df in DataFrame[df1, df2], col in quatcols])
            er[i, :err] = NA
        else
            p = Quaternion([df1[i, col] for col in quatcols]...)
            q = Quaternion([df2[i, col] for col in quatcols]...)
            er[i, :err] = err(p, q)
        end
    end
    er
end

end # module
