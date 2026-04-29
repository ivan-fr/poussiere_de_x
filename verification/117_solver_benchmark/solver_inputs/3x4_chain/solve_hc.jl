using HomotopyContinuation
@var x[1:3]
F = System([(0.36572469128584523-1.26760118292579im) + (1+0im)*x[1]^4 + (0.16508174810698892+0.071749679037196118im)*x[2],
    (-0.90244624466806955+0.44029874260899998im) + (1+0im)*x[2]^4 + (0.028818776079869626+0.17767801818249307im)*x[3],
    (-0.68077727010818401+1.2573132719701732im) + (1+0im)*x[3]^4])
t0 = time()
result = solve(F; show_progress=false)
sols = solutions(result)
println("HC_SOLUTIONS=", length(sols))
println("HC_TIME=", time() - t0)
