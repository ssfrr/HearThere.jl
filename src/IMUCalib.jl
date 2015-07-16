module IMUCalib

# currently using my hacked version of OSC.jl to enable 1-based indexing and
# the path method. PR has been filed
using OSC
using DataFrames

export collect_stationary

function collect_stationary(duration::Real, oscport::Integer=9999)
    timestamps = Float64[]
    data = [Float32[] for _ in 1:9]
    oscpath = "/rawimu"
    sock = UdpSocket()
    if !bind(sock, ip"0.0.0.0", oscport)
        println("Couldn't bind to port $oscport")
        return makedf(timestamps, data)
    end
    # wait for the first message
    println("Waiting for first message...")
    while true
        msg = OscMsg(recv(sock))
        println("Got path: $(path(msg))")
        if path(msg) == oscpath
            break
        end
    end
    println("Capturing Data...")
    starttime = time()
    now = time()
    while(now < starttime + duration)
        msg = OscMsg(recv(sock))
        now = time()
        if path(msg) == oscpath
            push!(timestamps, now)
            for i in 1:9
                push!(data[i], msg[i])
            end
        end
    end
    close(sock)
    makedf(timestamps, data)
end

function makedf(timestamps, data)
    DataFrame(
        timestamp=timestamps,
        accel_x=data[1],
        accel_y=data[2],
        accel_z=data[3],
        gyro_x=data[5],
        gyro_y=data[6],
        gyro_z=data[7],
        mag_x=data[7],
        mag_y=data[8],
        mag_z=data[9])
end

end # module
