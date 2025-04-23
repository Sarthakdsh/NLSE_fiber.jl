using FFTW


"""
    create_time_grid(t_min, t_max, N)

Create a time grid with N points from t_min to t_max.

# Arguments
- `t_min::Float64`: Minimum time [s]
- `t_max::Float64`: Maximum time [s]
- `N::Int`: Number of grid points

# Returns
- `t::Vector{Float64}`: Time grid
"""
function create_time_grid(t_min::Float64, t_max::Float64, N::Int)
    return range(t_min, t_max, length=N) |> collect
end



"""
    create_frequency_grid(t)

Create a frequency grid from a time grid.

# Arguments
- `t::Vector{Float64}`: Time grid

# Returns
- `ω::Vector{Float64}`: Angular frequency grid [rad/s]
"""
function create_frequency_grid(t::Vector{Float64})
    N = length(t)
    dt = t[2] - t[1]
    return fftshift(fftfreq(N, 2π / dt))
end

