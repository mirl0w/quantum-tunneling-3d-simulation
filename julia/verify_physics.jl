#=
verify_physics.jl

Sanity-checks the outputs of tunneling_solver.jl and hydrogen_orbitals.jl
against known, closed-form analytical results. This is the "did I get
real physics or just a pretty picture" check.

Checks performed:
  1. Norm conservation: total probability ∫|psi|^2 dV should stay ≈ 1.0
     at every saved timestep. The split-operator method is unitary, so
     any drift here indicates a numerical bug (grid too coarse, dt too
     large causing aliasing, etc.) rather than real physics.

  2. Tunneling transmission vs. the exact 1D analytical formula for a
     rectangular barrier (E < V0):
         T_analytical = 1 / (1 + V0^2 * sinh^2(k2*L) / (4*E*(V0-E)))
         k2 = sqrt(2*(V0-E))
     This is the textbook result (Griffiths, "Introduction to Quantum
     Mechanics"). Our simulation is 3D and the packet has finite width
     in momentum, so exact agreement isn't expected -- but the two
     numbers should be the same order of magnitude and move the same
     direction as you vary V0 or barrier width. Large disagreement
     (>2x) usually means the grid box is too small (packet interacting
     with its own periodic image) or dt is too large.

  3. Orbital normalization: ∫|psi_nlm|^2 dV over all space must equal 1
     for a properly normalized wavefunction. This directly checks that
     the analytical normalization constant in radial_wavefunction() is
     correct.

  4. Radial node count: the hydrogen radial wavefunction R_nl(r) has
     EXACTLY (n - l - 1) zero crossings (nodes) for r > 0 -- this is a
     hard theorem, not an approximation. Counting sign changes of R_nl
     and comparing to n-l-1 is a strong correctness check independent
     of the 3D grid/normalization code entirely.

Writes verification_report.txt summarizing all results with PASS/FAIL
flags (using generous tolerances appropriate for a discretized grid).
=#

using MAT
using Printf

include("hydrogen_orbitals.jl")  # reuse radial_wavefunction(), laguerre(), etc.
# NOTE: hydrogen_orbitals.jl's main() already ran when included, printing
# its own progress -- that's expected, it's just recomputing the grids.

function analytical_transmission(E, V0, L)
    if E >= V0
        return NaN  # formula below assumes classically forbidden regime
    end
    k2 = sqrt(2*(V0 - E))
    return 1 / (1 + (V0^2 * sinh(k2*L)^2) / (4*E*(V0-E)))
end

function count_radial_nodes(n, l; rmax=60.0, npts=20000)
    rs = range(1e-6, rmax, length=npts)
    vals = [radial_wavefunction(n, l, r) for r in rs]
    nodes = 0
    for i in 2:length(vals)
        if sign(vals[i]) != sign(vals[i-1]) && abs(vals[i]) > 1e-12
            nodes += 1
        end
    end
    return nodes
end

function main()
    lines = String[]
    push!(lines, "="^60)
    push!(lines, "PHYSICS VERIFICATION REPORT")
    push!(lines, "="^60)

    # ---- 1 & 2: Tunneling checks ----
    tpath = joinpath(@__DIR__, "tunneling_frames.mat")
    if isfile(tpath)
        T = matread(tpath)
        norm_frames = T["norm_frames"]
        norm_drift = maximum(abs.(norm_frames .- norm_frames[1]))
        push!(lines, "")
        push!(lines, "-- Norm conservation (split-operator unitarity) --")
        push!(lines, @sprintf("  Initial norm: %.6f", norm_frames[1]))
        push!(lines, @sprintf("  Max deviation over all frames: %.2e", norm_drift))
        normOK = norm_drift < 1e-3
        push!(lines, "  Result: " * (normOK ? "PASS (norm conserved to <0.1%)" : "FAIL (check dt/grid resolution)"))

        E = T["E_mean"]; V0 = T["V0"]; L = T["barrier_width"]
        Tsim = T["transmission_prob"]
        Tana = analytical_transmission(E, V0, L)
        push!(lines, "")
        push!(lines, "-- Tunneling transmission vs. 1D analytical formula --")
        push!(lines, @sprintf("  E = %.3f, V0 = %.3f, barrier width L = %.3f", E, V0, L))
        push!(lines, @sprintf("  Simulated (3D) transmission probability: %.4f", Tsim))
        push!(lines, @sprintf("  Analytical (1D, Griffiths formula) T:    %.4f", Tana))
        ratio = Tsim / Tana
        push!(lines, @sprintf("  Ratio sim/analytical: %.2f", ratio))
        tunnelOK = 0.2 <= ratio <= 5.0
        push!(lines, "  Result: " * (tunnelOK ? "PASS (same order of magnitude, as expected for 3D vs 1D)" :
                                              "CHECK (large disagreement -- try a bigger box or smaller dt)"))
    else
        push!(lines, "\n[!] tunneling_frames.mat not found -- run tunneling_solver.jl first.")
    end

    # ---- 3 & 4: Orbital checks ----
    opath = joinpath(@__DIR__, "orbitals.mat")
    if isfile(opath)
        O = matread(opath)
        orbital_defs = [("1s",1,0), ("2s",2,0), ("2pz",2,1), ("3dz2",3,2), ("4fz3",4,3)]
        push!(lines, "")
        push!(lines, "-- Orbital normalization (∫|psi|^2 dV should = 1) --")
        for (name, n, l) in orbital_defs
            psi = O["psi_$name"]
            coords = O["coords_$name"]
            dx = coords[2] - coords[1]
            dV = dx^3
            normval = sum(psi.^2) * dV
            ok = abs(normval - 1.0) < 0.05
            push!(lines, @sprintf("  %-6s n=%d l=%d : integral = %.4f  [%s]",
                  name, n, l, normval, ok ? "PASS" : "CHECK (grid may be too small/coarse)"))
        end

        push!(lines, "")
        push!(lines, "-- Radial node count (exact theorem: nodes = n - l - 1) --")
        for (name, n, l) in orbital_defs
            expected = n - l - 1
            actual = count_radial_nodes(n, l)
            ok = actual == expected
            push!(lines, @sprintf("  %-6s n=%d l=%d : expected %d nodes, found %d  [%s]",
                  name, n, l, expected, actual, ok ? "PASS" : "FAIL"))
        end
    else
        push!(lines, "\n[!] orbitals.mat not found -- run hydrogen_orbitals.jl first.")
    end

    push!(lines, "")
    push!(lines, "="^60)

    report = join(lines, "\n")
    println(report)
    open(joinpath(@__DIR__, "verification_report.txt"), "w") do f
        write(f, report)
    end
    println("\nSaved verification_report.txt")
end

main()
