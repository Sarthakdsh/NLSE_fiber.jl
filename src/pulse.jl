using FFTW
using LinearAlgebra


# export Pulse, gaussian_pulse, sech_pulse, energy, peak_power, fwhm

"""
A class representing an optical pulse in time domain.

# Fields
- `t`: Time grid [s]
- `field_t`: Complex field in time domain
"""
mutable struct Pulse
    t::Vector{Float64}       # Time grid [s]
    field_t::Vector{ComplexF64}  # Field in time domain
    """
        Pulse(t, field_t)
    
    Create a pulse from time domain data.
    
    # Arguments
    - `t::Vector{Float64}`: Time grid
    - `field_t::Vector{ComplexF64}`: Complex field in time domain
    """
    function Pulse(t::Vector{Float64}, field_t::Vector{ComplexF64})
        return new(t, field_t)
    end
end # Pulse


# ----------- Pulse Functions -----------

"""
    get_field_ω(pulse::Pulse)

Get the field in frequency domain.
"""
function get_field_ω(pulse::Pulse)
    N = length(pulse.t)
    dt = pulse.t[2] - pulse.t[1]
    return fftshift(ifft(pulse.field_t)) .* (N * dt / sqrt(2π))
end

"""
    energy(pulse::Pulse)

Calculate the pulse energy.
"""
function energy(pulse::Pulse)
    dt = pulse.t[2] - pulse.t[1]
    return dt * sum(abs2.(pulse.field_t))
end

"""
    peak_power(pulse::Pulse)

Calculate the pulse peak power.
"""
function peak_power(pulse::Pulse)
    return maximum(abs2.(pulse.field_t))
end

"""
    fwhm(pulse::Pulse)

Calculate the full width at half maximum of the pulse intensity.
"""
function fwhm(pulse::Pulse)
    # To do
end

# ---------- Derived Pulse ----------


"""
    gaussian_pulse(t, t0, power, chirp=zeros(Float64, length(t)))

Create a Gaussian pulse with optional chirp.

# Arguments
- `t::Vector{Float64}`: Time grid
- `t0::Float64`: Pulse duration (1/e intensity half-width)
- `power::Float64`: Peak power
- `chirp::Vector{Float64}`: Chirp parameter
"""
function gaussian_pulse(t::Vector{Float64}, t0::Float64, power::Float64, chirp::Vector{Float64}=zeros(Float64, length(t)))
    @assert length(t) == length(chirp) "t and chirp must have the same length"
    field_t = @. sqrt(power) * exp(-0.5 * (1 + 1im * chirp) * (t / t0)^2)
    return Pulse(t, field_t)
end

"""
    sech_pulse(t, t0, power, chirp=zeros(Float64, length(t)))

Create a Hyperbolic-Secant (Sech) pulse with optional chirp.

# Arguments
- `t::Vector{Float64}`: Time grid
- `t0::Float64`: Pulse duration
- `power::Float64`: Peak power
- `chirp::Vector{Float64}`: Chirp parameter
"""
function sech_pulse(t::Vector{Float64}, t0::Float64, power::Float64, chirp::Vector{Float64}=zeros(Float64, length(t)))
    @assert length(t) == length(chirp) "t and chirp must have the same length"
    field_t = @. sqrt(power) * sech(t / t0) * exp(-0.5 * (1im * chirp * t^2) / t0^2)
    return Pulse(t, field_t)
end

