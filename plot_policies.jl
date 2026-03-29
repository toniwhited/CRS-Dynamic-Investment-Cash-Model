using DelimitedFiles, Plots

outdir = joinpath(@__DIR__, "Outputcomp")

cg = vec(readdlm(joinpath(outdir, "cg.txt")))
zg = vec(readdlm(joinpath(outdir, "zg.txt")))

V = readdlm(joinpath(outdir, "v.txt"))[:, 2:end]
C = readdlm(joinpath(outdir, "c.txt"))[:, 2:end]
I = readdlm(joinpath(outdir, "i.txt"))[:, 2:end]

# V, C, I are (nz × nc). Surface plots want f(x, y) where x and y are vectors.
surface(cg, zg, V, xlabel="Cash/Capital", ylabel="z", zlabel="V(c,z)",
        title="Value Function", camera=(30, 30), colormap=:viridis)
savefig(joinpath(outdir, "value_function.png"))

surface(cg, zg, C, xlabel="Cash/Capital", ylabel="z", zlabel="c'(c,z)",
        title="Cash Policy", camera=(30, 60), colormap=:viridis)
savefig(joinpath(outdir, "cash_policy.png"))

surface(cg, zg, I, xlabel="Cash/Capital", ylabel="z", zlabel="i(c,z)",
        title="Investment Policy", camera=(30, 30), colormap=:viridis)
savefig(joinpath(outdir, "investment_policy.png"))

println("Plots saved to $outdir")
