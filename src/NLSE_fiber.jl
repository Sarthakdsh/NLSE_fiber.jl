module NLSE_fiber

# function demo(x=1, y=2)
#     println("Hello, World!")
#     return x + y
# end

export create_time_grid, create_frequency_grid
export Pulse, gaussian_pulse, sech_pulse, energy, peak_power, get_field_ω, fwhm
export plot_pulse

include("grid.jl")
include("pulse.jl")
include("visualization.jl")

end
