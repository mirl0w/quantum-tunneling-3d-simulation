#=
tunneling_solver.jl

3D Time-Dependent Schrödinger Equation solver using the split-operator
(split-step Fourier) method.

Physical setup (atomic units: ħ = m = 1):
    A 3D Gaussian wave packet is launched with momentum k0 along +x
    toward a rectangular potential barrier of height V0 and width w.
    Some probability tunnels through even though the packet's mean
    kinetic energy is LESS than V0 — this is the headline quantum effect.

Split-operator propagator (accurate to O(dt^3)):
    psi(t+dt) = exp(-i*V*dt/2) * IFFT( exp(-i*k^2*dt/2) * FFT( exp(-i*V*dt/2) * psi(t) ) )

Output:
    tunneling_frames.mat containing:
        psi2   -> (Nx,Ny,Nz,Nframes) array of |psi|^2
        x,y,z  -> coordinate vectors
        t      -> time vector (one entry per saved frame)
        V      -> (Nx,Ny,Nz) potential array (for drawing the barrier)
        barrier_x -> [x1 x2] barrier extent, for MATLAB slab plotting
=#

using FFTW
using MAT
using Printf
using LinearAlgebra

function build_grid(N::Int, L::Float64)
    x = range(-L/2, L/2, length=N) |> collect
    dx = x[2] - x[1]
    # FFT-ordered wavenumbers
    k = 2π * fftfreq(N, 1/dx)
    return x, k
end

function gaussian_wavepacket(x, y, z, x0, y0, z0, kx0, sigma)
    Nx, Ny, Nz = length(x), length(y), length(z)
    psi = Array{ComplexF64}(undef, Nx, Ny, Nz)
    norm_const = (1/(2π*sigma^2))^(3/4)
    for i in 1:Nx, j in 1:Ny, k in 1:Nz
        r2 = (x[i]-x0)^2 + (y[j]-y0)^2 + (z[k]-z0)^2
        envelope = norm_const * exp(-r2/(4*sigma^2))
        phase = exp(im*kx0*x[i])
        psi[i,j,k] = envelope * phase
    end
    return psi
end

function rectangular_barrier(x, y, z, V0, x1, x2)
    Nx, Ny, Nz = length(x), length(y), length(z)
    V = zeros(Float64, Nx, Ny, Nz)
    for i in 1:Nx
        if x1 <= x[i] <= x2
            V[i, :, :] .= V0
        end
    end
    return V
end

function run_tunneling_simulation(; N=64, L=40.0, dt=0.01, nsteps=800, save_every=20,
                                     x0=-10.0, kx0=2.0, sigma=1.5,
                                     V0=2.5, barrier_half_width=0.5)

    x, kx = build_grid(N, L)
    y, ky = build_grid(N, L)
    z, kz = build_grid(N, L)

    # 3D kinetic term  K = kx^2 + ky^2 + kz^2
    K2 = Array{Float64}(undef, N, N, N)
    for i in 1:N, j in 1:N, k in 1:N
        K2[i,j,k] = kx[i]^2 + ky[j]^2 + kz[k]^2
    end

    x1, x2 = -barrier_half_width, barrier_half_width
    V = rectangular_barrier(x, y, z, V0, x1, x2)

    # Mean kinetic energy of the packet, for a sanity printout
    E_mean = 0.5 * kx0^2
    @printf("Wave packet mean energy  E = %.3f\n", E_mean)
    @printf("Barrier height           V0 = %.3f  (E < V0 => classically forbidden)\n", V0)

    psi = gaussian_wavepacket(x, y, z, x0, 0.0, 0.0, kx0, sigma)
    psi ./= sqrt(sum(abs2.(psi)) * (x[2]-x[1]) * (y[2]-y[1]) * (z[2]-z[1]))  # normalize

    half_V_phase  = exp.(-im .* V .* (dt/2))
    kinetic_phase = exp.(-im .* K2 .* (dt/2))

    nframes = div(nsteps, save_every) + 1
    psi2_frames = Array{Float64}(undef, N, N, N, nframes)
    t_frames = Array{Float64}(undef, nframes)

    dV = (x[2]-x[1])*(y[2]-y[1])*(z[2]-z[1])

    psi2_frames[:,:,:,1] = abs2.(psi)
    t_frames[1] = 0.0
    norm_frames = Array{Float64}(undef, nframes)
    norm_frames[1] = sum(psi2_frames[:,:,:,1]) * dV
    frame_idx = 2

    for step in 1:nsteps
        psi .= half_V_phase .* psi
        psi_k = fft(psi)
        psi_k .= kinetic_phase .* psi_k
        psi .= ifft(psi_k)
        psi .= half_V_phase .* psi

        if step % save_every == 0
            psi2_frames[:,:,:,frame_idx] = abs2.(psi)
            t_frames[frame_idx] = step * dt
            norm_frames[frame_idx] = sum(psi2_frames[:,:,:,frame_idx]) * dV
            frame_idx += 1
            @printf("step %4d / %4d saved (frame %d/%d)\n", step, nsteps, frame_idx-1, nframes)
        end
    end

    # Report transmission probability (fraction of |psi|^2 past the barrier)
    dV = (x[2]-x[1])*(y[2]-y[1])*(z[2]-z[1])
    transmitted = sum(psi2_frames[x .> x2, :, :, end]) * dV
    @printf("\nApprox. tunneling transmission probability: %.4f\n", transmitted)

    matwrite(joinpath(@__DIR__, "tunneling_frames.mat"),
        Dict(
            "psi2" => psi2_frames,
            "x" => x, "y" => y, "z" => z,
            "t" => t_frames,
            "V" => V,
            "barrier_x" => [x1, x2],
            "transmission_prob" => transmitted,
            "norm_frames" => norm_frames,
            "E_mean" => E_mean,
            "V0" => V0,
            "barrier_width" => (x2 - x1)
        )
    )
    println("Saved tunneling_frames.mat")
end

run_tunneling_simulation()
