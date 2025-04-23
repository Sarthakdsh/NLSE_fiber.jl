# Example usage
using NLSE_fiber

function example()
    # Create time grid
    t_min = -10e-12  # -10 ps
    t_max = 10e-12   # 10 ps
    N = 2^12         # Number of points (power of 2 for FFT efficiency)
    t = create_time_grid(t_min, t_max, N)

    # Create a Gaussian pulse
    width = 1e-12  # 1 ps
    power = 1e3    # 1 kW peak power
    chirp = 0.5 * ones(Float64, N)    # Chirp parameter
    pulse = gaussian_pulse(t, width, power, chirp)

    # Calculate some properties
    e = energy(pulse)
    pp = peak_power(pulse)
    # width = fwhm(pulse)

    println("Pulse energy: ", e, " J")
    println("Peak power: ", pp, " W")
    # println("FWHM: ", width * 1e12, " ps")

    # Plot the pulse
    plot_pulse(pulse)

end