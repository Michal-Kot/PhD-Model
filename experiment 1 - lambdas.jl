# Hipoteza 1

include(pwd() * "\\methods\\methods.jl")

# ZAKOŃCZONY
# Wpływ komunikacji pomiędzy konsumentami na dynamikę systemu

#%% Eksperyment 1, λ_ind, λ_wom

ex1_v1_total_surplus = []
ex1_v1_producer_surplus = []
ex1_v1_consumer_surplus = []
ex1_v1_price = []
ex1_v1_quantity_produced = []
ex1_v1_quantity_sold = []
ex1_v1_quantity_leased = []
ex1_v1_producer_surplus_singleton = []

ex1_v2_total_surplus = []
ex1_v2_producer_surplus = []
ex1_v2_consumer_surplus = []
ex1_v2_price = []
ex1_v2_quantity_produced = []
ex1_v2_quantity_sold = []
ex1_v2_quantity_leased = []
ex1_v2_producer_surplus_singleton = []

ex1_v12_λ = []

for i in 1:750

    #if (mod(i,10) == 0) | (i == 1)
        println(i)
    #end

    λ_ind = rand()
    λ_wom = rand()

    push!(ex1_v12_λ, [λ_ind, λ_wom])

    ex1_v1_sim = TO_GO(200, 2, 200, 300, [0.4, 0.4], [1.1, 1.1], "random", λ_ind, λ_wom, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.5, 0.95], [0.5, 0.95]], [[0.8, 2.], [0.8, 2.]], 0.1, true, true, 1, [0.7, 1.0], "softmax", [true, false])

    push!(ex1_v1_total_surplus, calculate_surplus(ex1_v1_sim, "total", true))
    push!(ex1_v1_producer_surplus, calculate_surplus(ex1_v1_sim, "producer", true))
    push!(ex1_v1_consumer_surplus, calculate_surplus(ex1_v1_sim, "consumer,total",true))
    push!(ex1_v1_price, calculate_price_history.(ex1_v1_sim.sellers))
    push!(ex1_v1_quantity_produced, getfield.(ex1_v1_sim.sellers, :quantity_produced_history))
    push!(ex1_v1_quantity_sold, getfield.(ex1_v1_sim.sellers, :quantity_sold_history))
    push!(ex1_v1_quantity_leased, getfield.(ex1_v1_sim.sellers, :quantity_leased_history))
    push!(ex1_v1_producer_surplus_singleton, calculate_profit_history.(ex1_v1_sim.sellers))

end


ex3_p1 = plot_ecdf(true, mean.(getindex.(ex1_v1_producer_surplus_singleton, 1)), "Producent bada oczekiwania konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, mean.(getindex.(ex1_v1_producer_surplus_singleton, 2)), "Producent nie bada oczekiwań konsumentów")


L1 = getindex.(ex1_v12_λ, 1) 
L2 = getindex.(ex1_v12_λ, 2)

L1, L1u = cut_integer(L1, 5)
L2, L2u = cut_integer(L2, 5)

#### Total market surplus

hm_eq = [mean(ex1_v1_total_surplus[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'

ex4_p1 = StatsPlots.heatmap(L1u, L2u, hm_eq, xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość", titlefontsize = 8)

Plots.savefig(ex4_p1, pwd() * "\\Plots\\ex1_total surplus equal.svg")

#%% Producer 1 surplus, worse producer

hm_worse_eq = [mean(sum.(getindex.(ex1_v1_producer_surplus_singleton, 1))[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'

hm_better_eq = [mean(sum.(getindex.(ex1_v1_producer_surplus_singleton, 2))[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'

ex4_p3 = StatsPlots.heatmap(L1u, L2u, hm_better_eq, xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka producenta, identyczna jakość", titlefontsize = 8)

ex4 = StatsPlots.heatmap(L1u, L2u, hm_worse_eq, xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka producenta, identyczna jakość", titlefontsize = 8)

Plots.savefig(ex4_p3, pwd() * "\\plots\\ex1_prod surplus equal.svg")

function split_matrix(A)
    A_size = maximum.(axes(A))
    rows = Int(floor(A_size[1] / 2))
    cols = Int(floor(A_size[2] / 2))

    upper_square = A[1:rows, 1:cols]
    lower_square = A[(end - rows + 1):end, (end - cols + 1):end]

    return (mean(upper_square) - mean(lower_square))
end

split_matrix(hm_better_eq)
split_matrix(hm_worse_eq)



#%% Eksperyment 2, liczba powiązań

ex1_v3_total_surplus = []
ex1_v3_producer_surplus = []
ex1_v3_consumer_surplus = []
ex1_v3_price = []
ex1_v3_quantity_produced = []
ex1_v3_quantity_sold = []
ex1_v3_quantity_leased = []
ex1_v3_producer_surplus_singleton = []
ex1_v3_quality = []
ex1_v3_durability = []
ex1_v3_quality_exp = []
ex1_v3_durability_exp = []

ex1_v4_total_surplus = []
ex1_v4_producer_surplus = []
ex1_v4_consumer_surplus = []
ex1_v4_price = []
ex1_v4_quantity_produced = []
ex1_v4_quantity_sold = []
ex1_v4_quantity_leased = []
ex1_v4_producer_surplus_singleton = []
ex1_v4_quality = []
ex1_v4_durability = []
ex1_v4_quality_exp = []
ex1_v4_durability_exp = []

ex1_v34_nn = []

for i in 1:2000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    nn = sample(LinRange(0:25:400))

    push!(ex1_v34_nn, nn)

    ex1_v3_sim = TO_GO(250, 2, 200, nn, [0.5, 0.5], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [0.1, 0.1], [1.1, 1.1], [[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]],[[0.8,2.], [0.8,2.]], 0.1, true, true, 0)

    push!(ex1_v3_total_surplus, calculate_surplus(ex1_v3_sim, "total", true))
    push!(ex1_v3_producer_surplus, calculate_surplus(ex1_v3_sim, "producer", true))
    push!(ex1_v3_consumer_surplus, calculate_surplus(ex1_v3_sim, "consumer,total",true))
    push!(ex1_v3_price, calculate_price_history.(ex1_v3_sim.sellers))
    push!(ex1_v3_quantity_produced, getfield.(ex1_v3_sim.sellers, :quantity_produced_history))
    push!(ex1_v3_quantity_sold, getfield.(ex1_v3_sim.sellers, :quantity_sold_history))
    push!(ex1_v3_quantity_leased, getfield.(ex1_v3_sim.sellers, :quantity_leased_history))
    push!(ex1_v3_producer_surplus_singleton, calculate_profit_history.(ex1_v3_sim.sellers))
    push!(ex1_v3_quality, getfield.(ex1_v3_sim.sellers, :quality_history))
    push!(ex1_v3_durability, getfield.(ex1_v3_sim.sellers, :durability_history))
    push!(ex1_v3_quality_exp, mean([b.quality_expectation_history for b in ex1_v3_sim.buyers]))
    push!(ex1_v3_durability_exp, mean([b.durability_expectation_history for b in ex1_v3_sim.buyers]))

    ex1_v4_sim = TO_GO(250, 2, 200, nn, [0.5, 0.5], [1.0, 1.0], "random", 0.25, 0.25, "stochastic",[0.1, 0.1], [1.1, 1.1], [[0.9, 1.5], [0.6, 1.1]], [[0.2, 0.5], [0.3, 0.6]],[[0.8,2.], [0.8,2.]], 0.1, true, true, 0)
    push!(ex1_v4_total_surplus, calculate_surplus(ex1_v4_sim, "total", true))
    push!(ex1_v4_producer_surplus, calculate_surplus(ex1_v4_sim, "producer", true))
    push!(ex1_v4_consumer_surplus, calculate_surplus(ex1_v4_sim, "consumer,total",true))
    push!(ex1_v4_price, calculate_price_history.(ex1_v4_sim.sellers))
    push!(ex1_v4_quantity_produced, getfield.(ex1_v4_sim.sellers, :quantity_produced_history))
    push!(ex1_v4_quantity_sold, getfield.(ex1_v4_sim.sellers, :quantity_sold_history))
    push!(ex1_v4_quantity_leased, getfield.(ex1_v4_sim.sellers, :quantity_leased_history))
    push!(ex1_v4_producer_surplus_singleton, calculate_profit_history.(ex1_v4_sim.sellers))
    push!(ex1_v4_quality, getfield.(ex1_v4_sim.sellers, :quality_history))
    push!(ex1_v4_durability, getfield.(ex1_v4_sim.sellers, :durability_history))
    push!(ex1_v4_quality_exp, mean([b.quality_expectation_history for b in ex1_v4_sim.buyers]))
    push!(ex1_v4_durability_exp, mean([b.durability_expectation_history for b in ex1_v4_sim.buyers]))

end

Plots.plot(ex1_v3_quality[4][1])
Plots.plot!(getindex.(ex1_v3_quality_exp[4],1))

ex1_v3_quality_diff_1 = [sum((x .- y).^2) for (x,y) in zip(getindex.(ex1_v3_quality, 1), [getindex.(x,1) for x in ex1_v3_quality_exp])]
ex1_v3_quality_diff_2 = [sum((x .- y).^2) for (x,y) in zip(getindex.(ex1_v3_quality, 2), [getindex.(x,2) for x in ex1_v3_quality_exp])]
ex1_v4_quality_diff_1 = [sum((x .- y).^2) for (x,y) in zip(getindex.(ex1_v4_quality, 1), [getindex.(x,1) for x in ex1_v4_quality_exp])]
ex1_v4_quality_diff_2 = [sum((x .- y).^2) for (x,y) in zip(getindex.(ex1_v4_quality, 2), [getindex.(x,2) for x in ex1_v4_quality_exp])]

Plots.plot(sort(unique(ex1_v34_nn)), [mean(ex1_v3_quality_diff_1[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])
Plots.plot!(sort(unique(ex1_v34_nn)), [mean(ex1_v3_quality_diff_2[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])
Plots.plot!(sort(unique(ex1_v34_nn)), [mean(ex1_v4_quality_diff_1[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])
Plots.plot!(sort(unique(ex1_v34_nn)), [mean(ex1_v4_quality_diff_2[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])

ex1_v3_durability_diff_1 = [sum((x .- y).^2) for (x,y) in zip(getindex.(ex1_v3_durability, 1), [getindex.(x,1) for x in ex1_v3_durability_exp])]
ex1_v3_durability_diff_2 = [sum((x .- y).^2) for (x,y) in zip(getindex.(ex1_v3_durability, 2), [getindex.(x,2) for x in ex1_v3_durability_exp])]
ex1_v4_durability_diff_1 = [sum((x .- y).^2) for (x,y) in zip(getindex.(ex1_v4_durability, 1), [getindex.(x,1) for x in ex1_v4_durability_exp])]
ex1_v4_durability_diff_2 = [sum((x .- y).^2) for (x,y) in zip(getindex.(ex1_v4_durability, 2), [getindex.(x,2) for x in ex1_v4_durability_exp])]

Plots.plot(sort(unique(ex1_v34_nn)), [mean(ex1_v3_durability_diff_1[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))], xlabel = "Liczba połączeń w grafie", ylabel = "Różnica w oczekiwanej trwałości")
Plots.plot!(sort(unique(ex1_v34_nn)), [mean(ex1_v3_durability_diff_2[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])
Plots.plot!(sort(unique(ex1_v34_nn)), [mean(ex1_v4_durability_diff_1[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])
Plots.plot!(sort(unique(ex1_v34_nn)), [mean(ex1_v4_durability_diff_2[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])

StatsPlots.histogram(ex1_v34_nn)

Plots.plot(sort(unique(ex1_v34_nn)), [mean(ex1_v3_total_surplus[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])
Plots.plot!(sort(unique(ex1_v34_nn)), [mean(ex1_v4_total_surplus[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])

Plots.plot(sort(unique(ex1_v34_nn)), [mean(ex1_v3_consumer_surplus[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])
Plots.plot!(sort(unique(ex1_v34_nn)), [mean(ex1_v4_consumer_surplus[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])

Plots.plot(sort(unique(ex1_v34_nn)), [mean(ex1_v3_producer_surplus[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])
Plots.plot!(sort(unique(ex1_v34_nn)), [mean(ex1_v4_producer_surplus[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])

Plots.plot(sort(unique(ex1_v34_nn)), [mean(mean.(mean.(ex1_v3_quantity_sold))[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])
Plots.plot!(sort(unique(ex1_v34_nn)), [mean(ex1_v4_producer_surplus[ex1_v34_nn .== x]) for x in sort(unique(ex1_v34_nn))])
