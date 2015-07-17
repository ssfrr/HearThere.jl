module IMUCalib

# currently using my hacked version of OSC.jl to enable 1-based indexing and
# the path method. PR has been filed
using OSC
using DataFrames
using PyPlot

export getrawdata

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
        return makedf(timestamps, data)
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
    makedf(timestamps, data)
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
function makedf(timestamps, data)
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

end # module
