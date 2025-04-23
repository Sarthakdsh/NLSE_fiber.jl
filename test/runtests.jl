using NLSE_fiber
using Test
using FFTW
using Plots

@testset "NLSE_fiber.jl" begin
    @testset "Grid Functions" begin
        # Test time grid creation
        t_min = -1.0
        t_max = 1.0
        N = 1024
        t = create_time_grid(t_min, t_max, N)
        dt = t[2] - t[1]

        @test length(t) == N
        @test t[1] ≈ t_min
        @test t[end] ≈ t_max
        @test dt > 0
        @test issorted(t)

        # Test frequency grid creation
        ω = create_frequency_grid(t)
        dω = ω[2] - ω[1]

        @test length(ω) == N
        @test issorted(ω)
        @test sum(ω) ≈ ω[1]  # Should be symmetric around zero for odd N and for even N 1 term will be extra
        @test dω > 0

        @test dt * dω ≈ 2π / N   # Should satisfy the uncertainty relation

    end

    @testset "Pulse Creation" begin
        t = create_time_grid(-10e-12, 10e-12, 1024)
        t0 = 1e-12  # 1 ps
        power = 1e3  # 1 kW

        # Test Gaussian pulse
        pulse_gauss = gaussian_pulse(t, t0, power)
        @test length(pulse_gauss.field_t) == length(t)
        @test pulse_gauss.t == t

        # Test chirped Gaussian pulse
        chirp = ones(length(t))
        pulse_chirped = gaussian_pulse(t, t0, power, chirp)
        @test length(pulse_chirped.field_t) == length(t)

        # Test Sech pulse
        pulse_sech = sech_pulse(t, t0, power)
        @test length(pulse_sech.field_t) == length(t)
        @test pulse_sech.t == t
    end

    @testset "Pulse Analysis" begin
        t = create_time_grid(-10e-12, 10e-12, 1024)
        t0 = 1e-12
        power = 1e3
        pulse = gaussian_pulse(t, t0, power)

        # Test energy calculation
        E = energy(pulse)
        @test E > 0
        @test typeof(E) == Float64

        # Test peak power calculation
        P = peak_power(pulse)
        @test abs(P - power) / power < 1e-3  # Allow some tolerance
        @test typeof(P) == Float64

        # Test frequency domain conversion
        field_ω = get_field_ω(pulse)
        @test length(field_ω) == length(pulse.field_t)
        @test all(isfinite.(field_ω))
    end

    @testset "Error Handling" begin
        t = create_time_grid(-10e-12, 10e-12, 1024)
        t0 = 1e-12
        power = 1e3

        # Test chirp length mismatch
        wrong_chirp = ones(length(t) + 1)
        @test_throws AssertionError gaussian_pulse(t, t0, power, wrong_chirp)
        @test_throws AssertionError sech_pulse(t, t0, power, wrong_chirp)
    end

    @testset "Visualization" begin
        @testset "plot_pulse" begin
            # Create a test pulse
            t = create_time_grid(-10e-12, 10e-12, 1024)
            t0 = 1e-12
            power = 1e3
            pulse = gaussian_pulse(t, t0, power)

            # Test normal plotting
            p = plot_pulse(pulse)
            @test p isa Plots.Plot
            @test length(p.subplots) == 2  # Should have 2 subplots


            # Test invalid pulse
            invalid_pulse = Pulse(Float64[], Complex{Float64}[])
            @test_throws ErrorException plot_pulse(invalid_pulse)
        end
    end
end
