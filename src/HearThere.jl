module HearThere

# currently using my hacked version of OSC.jl to enable 1-based indexing and
# the path method. PR has been filed
using OSC
using DataFrames
using PyPlot
using Quaternions
using LightXML
using DataFramesMeta

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

    /orientation ffff q0 q1 q2 q3
and
    /ranges ffff anchor0 anchor1 anchor2 anchor3

and return a tuple of DataFrames of the timestamped values.
"""
function getcookeddata(duration::Real, oscport::Integer=10001)
    orientationpath = "/orientation"
    rangepath = "/ranges"
    rangedf = DataFrame(
        timestamp=Float64[],
        anchor0=Float32[],
        anchor1=Float32[],
        anchor2=Float32[],
        anchor3=Float32[])
    orientationdf = DataFrame(
        timestamp=Float64[],
        q0=Float32[],
        q1=Float32[],
        q2=Float32[],
        q3=Float32[])
    sock = UdpSocket()
    if !bind(sock, ip"0.0.0.0", oscport)
        println("Couldn't bind to port $oscport")
        return (rangedf, orientationdf)
    end
    try
        # wait for the first message
        println("Waiting for first message...")
        while true
            msg = OscMsg(recv(sock))
            if path(msg) in [orientationpath, rangepath]
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
            if path(msg) == orientationpath
                push!(orientationdf, [now, msg[1], msg[2], msg[3], msg[4]])
            elseif path(msg) == rangepath
                push!(rangedf, [now, msg[1], msg[2], msg[3], msg[4]])
            end
        end
    finally
        close(sock)
    end
    (rangedf, orientationdf)
end

function rotate(q::Quaternion, v::Vector)
    q = normalize(q)

    imag(q*Quaternion(v)*conj(q))
end

function plotpath(data)
    PyPlot.plot3D(data[:x], data[:y], data[:z])
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
    bodies = Dict{Int, String}()
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
        for (id, x, y, z, q0, q1, q2, q3) in zip(
                get_elements_by_tagname(bodydataElement, "RigidBodyID"),
                get_elements_by_tagname(bodydataElement, "x"),
                get_elements_by_tagname(bodydataElement, "y"),
                get_elements_by_tagname(bodydataElement, "z"),
                get_elements_by_tagname(bodydataElement, "qw"),
                get_elements_by_tagname(bodydataElement, "qx"),
                get_elements_by_tagname(bodydataElement, "qy"),
                get_elements_by_tagname(bodydataElement, "qz"))
            push!(df, [
                timestamp,
                int(content(id)),
                float(content(x)),
                float(content(y)),
                float(content(z)),
                float(content(q0)),
                float(content(q1)),
                float(content(q2)),
                float(content(q3))])
        end
    end
    # The OptiTrack reports a location of 0, 0, 0 if the trackable isn't present
    @where df ((:x .!= 0) | (:y .!= 0) | (:z .!= 0))
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
        # divide by 1000 to convert mm to m
        anchorDF[id+1, :x] += 2x / 1000
        anchorDF[id+1, :y] += 2y / 1000
        anchorDF[id+1, :z] += 2z / 1000
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

end # module
