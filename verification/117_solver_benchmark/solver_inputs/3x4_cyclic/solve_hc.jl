using HomotopyContinuation
@var x[1:3]
F = System([(-0.80740582853223486-0.32165471124566825im) + (1+0im)*x[1]^4 + (0.16508174810698892+0.071749679037196118im)*x[2],
    (-0.94624681390357368+0.99280755863269787im) + (1+0im)*x[2]^4 + (0.028818776079869626+0.17767801818249307im)*x[3],
    (0.67688917150439254-0.49690207506311412im) + (1+0im)*x[3]^4 + (-0.13394004577048665+0.12025000681496836im)*x[1]])
t0 = time()
result = solve(F; show_progress=false)
sols = solutions(result)
println("HC_SOLUTIONS=", length(sols))
println("HC_TIME=", time() - t0)
