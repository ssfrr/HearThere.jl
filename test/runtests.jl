module HearThereTests

using HearThere
using DataFrame
using Base.Test

function generate_rawdata(n::Integer)
    generate_random_rawdata(n, 100, 0, 100)
end

function generate_rawdata(n::Integer, μ::Real, σ::Real)
    generate_random_rawdata(n, fill(μ, 9), fill(σ, 9))
end

function generate_rawdata(n::Integer, μ::Array{Real}, σ::Array{Real})
    DataFrame(
        timestamp=linspace(0, 10, n),
        accel_x=convert(Array{Int32}, map(round, randn(n) * σ[1] + μ[1])),
        accel_y=convert(Array{Int32}, map(round, randn(n) * σ[2] + μ[2])),
        accel_z=convert(Array{Int32}, map(round, randn(n) * σ[3] + μ[3])),
        gyro_x=convert(Array{Int32}, map(round, randn(n) * σ[4] + μ[4])),
        gyro_y=convert(Array{Int32}, map(round, randn(n) * σ[5] + μ[5])),
        gyro_z=convert(Array{Int32}, map(round, randn(n) * σ[6] + μ[6])),
        mag_x=convert(Array{Int32}, map(round, randn(n) * σ[7] + μ[7])),
        mag_y=convert(Array{Int32}, map(round, randn(n) * σ[8] + μ[8])),
        mag_z=convert(Array{Int32}, map(round, randn(n) * σ[9] + μ[9])))
end

function test_gyro_bias()
    n = 1000
    df = generate_rawdata(n, 2, 100)
end

test_gyro_bias()

end
