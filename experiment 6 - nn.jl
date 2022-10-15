# Experiment 3

include(pwd() * "\\methods\\methods.jl")

# Impact on social communication of simulation dynamics

ex6_v1_total_surplus = []
ex6_v1_producer_surplus = []
ex6_v1_consumer_surplus = []
ex6_v1_price = []
ex6_v1_quantity_produced = []
ex6_v1_quantity_sold = []
ex6_v1_quantity_leased = []
ex6_v1_reselling = []
ex6_v1_producer_surplus_singleton = []
ex6_v1_buying_history = []
ex6_v1_quality = []
ex6_v1_durability = []
ex6_v1_margin = []

ex6_v2_total_surplus = []
ex6_v2_producer_surplus = []
ex6_v2_consumer_surplus = []
ex6_v2_price = []
ex6_v2_quantity_produced = []
ex6_v2_quantity_sold = []
ex6_v2_quantity_leased = []
ex6_v2_reselling = []
ex6_v2_producer_surplus_singleton = []
ex6_v2_buying_history = []
ex6_v2_quality = []
ex6_v2_durability = []
ex6_v2_margin = []

ex6_ss = []

for i in 1:500

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    ss = sample(LinRange(0.01:0.01:0.10))
    push!(ex6_ss, ss)

    ex6_v1_sim = TO_GO(200, 2, 200, 300, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.5, 0.95], [0.5, 0.95]], [[.8, 2.], [.8, 2.]], 0.10, true, true, 1, [0.7, 1.], "softmax", [true, false], ss)

    push!(ex6_v1_total_surplus, calculate_surplus(ex6_v1_sim, "total", true))
    push!(ex6_v1_producer_surplus, calculate_surplus(ex6_v1_sim, "producer", true))
    push!(ex6_v1_consumer_surplus, calculate_surplus(ex6_v1_sim, "consumer,total",true))
    push!(ex6_v1_quality, getfield.(ex6_v1_sim.sellers, :quality_history))
    push!(ex6_v1_durability, getfield.(ex6_v1_sim.sellers, :durability_history))
    push!(ex6_v1_margin, getfield.(ex6_v1_sim.sellers, :margin_history))
    push!(ex6_v1_price, calculate_price_history.(ex6_v1_sim.sellers))
    push!(ex6_v1_quantity_produced, getfield.(ex6_v1_sim.sellers, :quantity_produced_history))
    push!(ex6_v1_quantity_sold, getfield.(ex6_v1_sim.sellers, :quantity_sold_history))
    push!(ex6_v1_quantity_leased, getfield.(ex6_v1_sim.sellers, :quantity_leased_history))
    push!(ex6_v1_producer_surplus_singleton, calculate_profit_history.(ex6_v1_sim.sellers))
    push!(ex6_v1_reselling, getfield.(ex6_v1_sim.sellers, :reselling_history))
    push!(ex6_v1_buying_history, getfield.(ex6_v1_sim.buyers, :unit_buying_selling_history))

end

countmap(ex6_ss)

ex6_p1 = Plots.scatter(sort(unique(ex6_ss)), [mean(sum.(getindex.(ex6_v1_producer_surplus_singleton, 1))[ex6_ss .== ss]) for ss in sort(unique(ex6_ss))], label = "Nadwyżka producenta", legend = :bottomleft, title = "Nadwyżka producena a szczegółowość badania")

ss_model = GLM.lm(@formula(y~logx), DataFrame(y = [mean(sum.(getindex.(ex6_v1_producer_surplus_singleton, 1))[ex6_ss .== ss]) for ss in sort(unique(ex6_ss))], logx = log.(sort(unique(ex6_ss)))))

Plots.plot!(sort(unique(ex6_ss)), fitted(ss_model), xlabel = "% populacji w badaniu", ylabel = "Zysk producenta wykorzystującego badanie", label = "Dopasowana funkcja logarytmiczna")

Plots.plot!(sort(unique(ex6_ss)), fitted(ss_model) .- sort(unique(ex6_ss)) * 300 * 0.4, label = "Zysk producenta - koszt badania")

Plots.savefig(ex6_p1, pwd() * "\\Plots\\research perc vs income.svg")

Plots.scatter(sort(unique(ex6_ss)), [mean(mean.(getindex.(ex6_v1_producer_surplus_singleton, 2))[ex6_ss .== ss]) for ss in sort(unique(ex6_ss))], smooth = true)

Plots.scatter(ex6_nn, mean.(getindex.(ex6_v1_producer_surplus_singleton, 1)), smooth = true)
Plots.scatter(ex6_nn, mean.(getindex.(ex6_v1_producer_surplus_singleton, 2)), smooth = true)

e_k_all = getindex.(getfield.(sim_single.buyers, :quality_expectation), 1)
e_d_all = getindex.(getfield.(sim_single.buyers, :durability_expectation), 1)
ρ_all = getfield.(sim_single.buyers, :future_discount)

kdρ = [x for x in zip(e_k_all,e_d_all,ρ_all)]

ss = 0.5

kdρ_idx = sample(1:length(kdρ), Int(ss * 300), replace = false)

observed_kdρ = kdρ[kdρ_idx]

assumed_kdρ = vcat(fill(observed_kdρ, Int(1 / ss))...)

assumed_kdρ[1]
assumed_kdρ[151]

e_k = getindex.(assumed_kdρ, 1)
e_d = getindex.(assumed_kdρ, 2)
ρ_mean = getindex.(assumed_kdρ, 3)

mean(e_k)

Plots.plot(e_k)