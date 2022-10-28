# EX1 CONT'D

include(pwd() * "\\methods\\methods.jl")

# 

ex2_v1_total_surplus = []
ex2_v1_producer_surplus = []
ex2_v1_consumer_surplus = []
ex2_v1_price = []
ex2_v1_quantity_produced = []
ex2_v1_quantity_sold = []
ex2_v1_quantity_leased = []
ex2_v1_reselling = []
ex2_v1_producer_surplus_singleton = []
ex2_v1_buying_history = []
ex2_v1_quality = []
ex2_v1_durability = []
ex2_v1_margin = []
ex2_v1_quality_exp = []
ex2_v1_durability_exp = []
ex2_v1_quality_std = []
ex2_v1_durability_std = []

ex2_v1_pl = []

for i in 1:500

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    pl = sample(1:10)
    push!(ex2_v1_pl, pl)

    ex2_v1_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], 0.10, true, 1, [0.7, 1.], "softmax", [false, true], [0, 0.1], pl, false)
    push!(ex2_v1_total_surplus, calculate_surplus(ex2_v1_sim, "total", true))
    push!(ex2_v1_producer_surplus, calculate_surplus(ex2_v1_sim, "producer", true))
    push!(ex2_v1_consumer_surplus, calculate_surplus(ex2_v1_sim, "consumer,total",true))
    push!(ex2_v1_quality, trim_first.(getfield.(ex2_v1_sim.sellers, :quality_history); trimmed = 5))
    push!(ex2_v1_durability, trim_first.(getfield.(ex2_v1_sim.sellers, :durability_history); trimmed = 5))
    push!(ex2_v1_margin, trim_first.(getfield.(ex2_v1_sim.sellers, :margin_history); trimmed = 5))
    push!(ex2_v1_price, trim_first.(calculate_price_history.(ex2_v1_sim.sellers; product_life = 5); trimmed = 5))
    push!(ex2_v1_quantity_produced, trim_first.(getfield.(ex2_v1_sim.sellers, :quantity_produced_history); trimmed = 5))
    push!(ex2_v1_quantity_sold, trim_first.(getfield.(ex2_v1_sim.sellers, :quantity_sold_history); trimmed = 5))
    push!(ex2_v1_producer_surplus_singleton, trim_first.(calculate_profit_history.(ex2_v1_sim.sellers); trimmed = 5))
    push!(ex2_v1_reselling, trim_first.(getfield.(ex2_v1_sim.sellers, :reselling_history); trimmed = 5))
    push!(ex2_v1_buying_history, trim_first.(getfield.(ex2_v1_sim.buyers, :unit_buying_selling_history); trimmed = 5))
    push!(ex2_v1_quality_exp, trim_first.(mean([b.quality_expectation_history for b in ex2_v1_sim.buyers]); trimmed = 5))
    push!(ex2_v1_durability_exp, trim_first.(mean([b.durability_expectation_history for b in ex2_v1_sim.buyers]); trimmed = 5))
    push!(ex2_v1_quality_std, [std(mean.([getindex.(x,1) for x in getfield.(ex2_v1_sim.buyers, :quality_expectation_history)])), std(mean.([getindex.(x,2) for x in getfield.(ex2_v1_sim.buyers, :quality_expectation_history)]))])
    push!(ex2_v1_durability_std, [std(mean.([getindex.(x,1) for x in getfield.(ex2_v1_sim.buyers, :durability_expectation_history)])), std(mean.([getindex.(x,2) for x in getfield.(ex2_v1_sim.buyers, :durability_expectation_history)]))])

end

ex2_v1_pl

[mean(sum.(getindex.(ex2_v1_producer_surplus_singleton, 1))[ex2_v1_pl .== x]) for x in sort(unique(ex2_v1_pl))]

ex1_p10 = StatsPlots.boxplot([sum.(getindex.(ex2_v1_producer_surplus_singleton, 1))[ex2_v1_pl .== x] for x in sort(unique(ex2_v1_pl))], legend = false, xlabel = "Maksymalna przydatność dobra [liczba okresów]", ylabel = "Zysk firmy", title = "Zysk firmy, która nie bada konsumentów", ylim = (-1500, 1500))

Plots.savefig(ex1_p10, pwd() * "\\plots\\ex1\\profit per life cycle no research.svg")

ex1_p11 = StatsPlots.boxplot([sum.(getindex.(ex2_v1_producer_surplus_singleton, 2))[ex2_v1_pl .== x] for x in sort(unique(ex2_v1_pl))], legend = false, xlabel = "Maksymalna przydatność dobra [liczba okresów]", ylabel = "Zysk firmy", title = "Zysk firmy, która bada konsumentów", ylim = (-1500, 1500))

Plots.savefig(ex1_p11, pwd() * "\\plots\\ex1\\profit per life cycle research.svg")

"""Value-at-risk."""
function value_at_risk(x::Vector{Float64}, f::Vector{Float64}, α::Float64)
    i = findfirst(p -> p≥α, cumsum(f))
    if i === nothing
        return x[end]
    else
        return x[i]
    end
end

[percentile(sum.(getindex.(ex2_v1_producer_surplus_singleton, 1))[ex2_v1_pl .== x], 5) for x in sort(unique(ex2_v1_pl))]

[percentile(sum.(getindex.(ex2_v1_producer_surplus_singleton, 2))[ex2_v1_pl .== x], 5) for x in sort(unique(ex2_v1_pl))]

normalize(v) = v ./ sum(v)
scale(v, low, high) = v * (high - low) + low
n = 10
x = sort(scale.(rand(n), -1.0, 1.0))
f = normalize(rand(n))
α = 0.05

@assert issorted(x)
@assert all(f .≥ 0)
@assert sum(f) ≈ 1
@assert 0 ≤ α ≤ 1

value_at_risk(x,f,α)

Plots.scatter(ex2_v1_pl, sum.(getindex.(ex2_v1_producer_surplus_singleton, 1)), )
Plots.scatter(ex2_v1_pl, sum.(getindex.(ex2_v1_producer_surplus_singleton, 2)), )
Plots.plot(sort(unique(ex2_v1_pl)), [mean(sum.(getindex.(ex2_v1_producer_surplus_singleton, 1))[ex2_v1_pl .== x]) for x in sort(unique(ex2_v1_pl))])
Plots.plot!(sort(unique(ex2_v1_pl)), [mean(sum.(getindex.(ex2_v1_producer_surplus_singleton, 2))[ex2_v1_pl .== x]) for x in sort(unique(ex2_v1_pl))])

ex2_p1 = plot_ecdf(true, mean.(ex2_v1_total_surplus), "Nie patrzymy na konsumentów", xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false, mean.(ex2_v2_total_surplus), "Patrzymy na konsumentów")

ex2_p1 = plot_ecdf(true, mean.(ex2_v1_consumer_surplus), "Nie patrzymy na konsumentów", xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false, mean.(ex2_v2_consumer_surplus), "Patrzymy na konsumentów")

ex2_p3 = ex2_p1 = plot_ecdf(true, ex2_v1_producer_surplus, "Nie patrzymy na konsumentów", xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false, ex2_v2_producer_surplus, "Patrzymy na konsumentów")

ex2_p3 = plot_ecdf(sum.(sum.(ex2_v1_quantity_produced)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex2_v2_quantity_produced)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(sum.(sum.(ex2_v1_quantity_sold)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex2_v2_quantity_sold)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)



mean.(mean.(ex2_v1_price))

ex2_p3 = plot_ecdf(mean.(mean.(ex2_v1_price)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex2_v2_quantity_sold)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)