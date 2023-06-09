#######################################################################################################

# Exemplary usage of the model

#######################################################################################################


include(pwd() * "\\methods\\methods.jl")

@time sim_single = TO_GO(500, 2, 200, 2000, [0.6, 0.6], [1.1, 1.1], "random", .8, .8, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[1., 2.], [1., 2.]], 0.10, true, 1, [0.85, 1.], "softmax", ["market research", "market research"], [0.1, 0.1], 3, "dist", TriangularDist(0,1,0.5))

quality_expectation_buyers = [[] for _ = 1:500] # to T

for b in sim_single.buyers
    ubsh = b.unit_buying_selling_history
    for item in ubsh
        if item.decision == "buy, primary market"
            qe = getindex.(b.quality_expectation_history, 1)[item.t]
            push!(quality_expectation_buyers[item.t], qe)
        end
    end
end

p = Plots.plot(xlabel = "T", ylabel = "Jakość / przekonanie o jakości", legend = :bottomleft)

for i in 1:200
    p = Plots.plot!(getindex.(sim_single.buyers[i].quality_expectation_history, 1), color = "grey", linealpha = 0.10, label = nothing)
end

p

Plots.plot!(sim_single.sellers[1].quality_history, label = "Średnia jakość, producent", linewidth = 2, color = "blue")
Plots.plot!(mean([getindex.(x,1) for x in getfield.(sim_single.buyers, :quality_expectation_history)]), color = "red", linewidth = 2, label = "Przekonanie o jakości, cała populacja")
Plots.plot!(mean_nothing.(quality_expectation_buyers), label = "Przekonanie o jakości, kupujący w t", linewidth = 2, color = "orange")