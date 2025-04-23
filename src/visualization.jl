using Plots


# export plot_pulse

"""
    plot_pulse(pulse::Pulse)

Create a plot of the pulse in both time and frequency domains.
"""
function plot_pulse(pulse::Pulse)
    # Check if pulse is invalid (empty containts)
    if isempty(pulse.t) || isempty(pulse.field_t)
        error("Pulse does not contain any data. Please provide a valid pulse.")
    end

    # Time domain plot
    p1 = plot(pulse.t, abs2.(pulse.field_t),
        xlabel="Time (s)", ylabel="Power (W)",
        title="Temporal Profile", legend=false)

    # Frequency domain plot
    ω = create_frequency_grid(pulse.t)
    field_ω = get_field_ω(pulse)
    p2 = plot(ω, abs2.(field_ω),
        xlabel="Angular frequency(rad/s)", ylabel="Spectral Power (a.u.)",
        title="Spectral Profile", legend=false)

    plot(p1, p2, layout=(2, 1), size=(800, 600))
end