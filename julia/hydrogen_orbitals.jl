#=
hydrogen_orbitals.jl

Computes the EXACT analytical hydrogen-atom wavefunctions
    psi_nlm(r, theta, phi) = R_nl(r) * Y_lm(theta, phi)
on a 3D Cartesian grid, for a gallery of famous orbitals
(1s, 2s, 2p_z, 3d_z2, 4f_z3), and saves the real part and |psi|^2
to orbitals.mat for MATLAB isosurface rendering.

R_nl uses the generalized (associated) Laguerre polynomials L_{n-l-1}^{2l+1}.
Y_lm uses the associated Legendre polynomials P_l^m.
We use REAL orbitals (the standard chemistry/physics textbook shapes with
lobes), built from linear combinations of complex Y_lm, since those are
what give the recognizable "dumbbell"/"cloverleaf" shapes.

Units: Bohr radius a0 = 1 (atomic units).
=#

using MAT
using Printf

# ---- Generalized Laguerre polynomial L_n^alpha(x) via recurrence ----
function laguerre(n::Int, alpha::Float64, x::Float64)
    if n == 0
        return 1.0
    elseif n == 1
        return 1.0 + alpha - x
    end
    Lkm2 = 1.0
    Lkm1 = 1.0 + alpha - x
    Lk = 0.0
    for k in 2:n
        Lk = ((2k - 1 + alpha - x) * Lkm1 - (k - 1 + alpha) * Lkm2) / k
        Lkm2, Lkm1 = Lkm1, Lk
    end
    return Lk
end

# ---- Associated Legendre polynomial P_l^m(x), m >= 0, via recurrence ----
function assoc_legendre(l::Int, m::Int, x::Float64)
    # Standard stable recurrence (Numerical Recipes style)
    pmm = 1.0
    if m > 0
        somx2 = sqrt((1.0 - x) * (1.0 + x))
        fact = 1.0
        for i in 1:m
            pmm *= -fact * somx2
            fact += 2.0
        end
    end
    if l == m
        return pmm
    end
    pmmp1 = x * (2m + 1) * pmm
    if l == m + 1
        return pmmp1
    end
    pll = 0.0
    for ll in (m+2):l
        pll = ((2ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m)
        pmm, pmmp1 = pmmp1, pll
    end
    return pll
end

factorial_f(n::Int) = n <= 1 ? 1.0 : prod(Float64.(2:n))

function radial_wavefunction(n::Int, l::Int, r)
    a0 = 1.0
    rho = 2r / (n * a0)
    norm_const = sqrt((2/(n*a0))^3 * factorial_f(n - l - 1) / (2n * factorial_f(n + l)))
    L = laguerre(n - l - 1, Float64(2l + 1), rho)
    return norm_const * exp(-rho/2) * rho^l * L
end

function real_spherical_harmonic(l::Int, m::Int, theta, phi)
    # Real (tesseral) spherical harmonics, standard normalization
    absm = abs(m)
    Plm = assoc_legendre(l, absm, cos(theta))
    norm_const = sqrt((2l+1)/(4π) * factorial_f(l-absm) / factorial_f(l+absm))
    if m == 0
        return norm_const * Plm
    elseif m > 0
        return sqrt(2) * norm_const * Plm * cos(m*phi)
    else
        return sqrt(2) * norm_const * Plm * sin(absm*phi)
    end
end

function orbital_grid(n::Int, l::Int, m::Int; N=80, L=30.0)
    coords = range(-L/2, L/2, length=N) |> collect
    psi = Array{Float64}(undef, N, N, N)
    for i in 1:N, j in 1:N, k in 1:N
        xv, yv, zv = coords[i], coords[j], coords[k]
        r = sqrt(xv^2 + yv^2 + zv^2)
        theta = r < 1e-9 ? 0.0 : acos(clamp(zv/r, -1.0, 1.0))
        phi = atan(yv, xv)
        R = radial_wavefunction(n, l, r)
        Y = real_spherical_harmonic(l, m, theta, phi)
        psi[i,j,k] = R * Y
    end
    return coords, psi
end

function main()
    # (name, n, l, m) — a gallery spanning s, p, d, f shapes
    orbitals = [
        ("1s",   1, 0,  0),
        ("2s",   2, 0,  0),
        ("2pz",  2, 1,  0),
        ("3dz2", 3, 2,  0),
        ("4fz3", 4, 3,  0),
    ]

    result = Dict{String,Any}()
    coords_ref = nothing

    for (name, n, l, m) in orbitals
        # larger orbitals need a bigger box; scale roughly with n^2
        box = 14.0 * n^2 / 1  # heuristic extent in Bohr radii
        box = clamp(box, 20.0, 90.0)
        @printf("Computing orbital %s (n=%d, l=%d, m=%d) on box size %.1f a0\n", name, n, l, m, box)
        coords, psi = orbital_grid(n, l, m; N=70, L=box)
        coords_ref = coords
        result["psi_$name"]  = psi
        result["coords_$name"] = coords
    end

    result["orbital_names"] = [o[1] for o in orbitals]
    matwrite(joinpath(@__DIR__, "orbitals.mat"), result)
    println("Saved orbitals.mat")
end

main()
