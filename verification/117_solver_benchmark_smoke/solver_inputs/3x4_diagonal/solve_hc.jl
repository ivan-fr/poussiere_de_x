using HomotopyContinuation
@var x[1:3]
F = System([(-0.091814798735985914-1.0908336690405624im) + (1+0im)*x[1]^4,
    (-0.82925683810487183+0.067086120166431631im) + (1+0im)*x[2]^4,
    (-0.13203282344697168-1.2573505938393033im) + (1+0im)*x[3]^4])
t0 = time()
result = solve(F; show_progress=false)
sols = solutions(result)
println("HC_SOLUTIONS=", length(sols))
println("HC_TIME=", time() - t0)
