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
ex1_v1_reselling = []
ex1_v1_producer_surplus_singleton = []
ex1_v1_buying_history = []
ex1_v1_quality = []
ex1_v1_durability = []
ex1_v1_margin = []
ex1_v1_quality_exp = []
ex1_v1_durability_exp = []
ex1_v1_quality_std = []
ex1_v1_durability_std = []

ex1_v12_λ = []

JULIA_NUM_THREADS = 8

for i in 1:1000

    #if (mod(i,10) == 0) | (i == 1)
        println(i)
    #end

    λ_ind = sample(0.1:0.1:0.9)
    λ_wom = sample(0.1:0.1:0.9)

    push!(ex1_v12_λ, [λ_ind, λ_wom])

    ex1_v1_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1., 1.], "random", λ_ind, λ_wom, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], 0.50, true, 1, [0.7, 1.], "softmax", [false, true], [0, 0.1], sample(2:8), false)

    push!(ex1_v1_total_surplus, calculate_surplus(ex1_v1_sim, "total", true))
    push!(ex1_v1_producer_surplus, calculate_surplus(ex1_v1_sim, "producer", true))
    push!(ex1_v1_consumer_surplus, calculate_surplus(ex1_v1_sim, "consumer,total",true))
    push!(ex1_v1_quality, trim_first.(getfield.(ex1_v1_sim.sellers, :quality_history); trimmed = 5))
    push!(ex1_v1_durability, trim_first.(getfield.(ex1_v1_sim.sellers, :durability_history); trimmed = 5))
    push!(ex1_v1_margin, trim_first.(getfield.(ex1_v1_sim.sellers, :margin_history); trimmed = 5))
    push!(ex1_v1_price, trim_first.(calculate_price_history.(ex1_v1_sim.sellers; product_life = 5); trimmed = 5))
    push!(ex1_v1_quantity_produced, trim_first.(getfield.(ex1_v1_sim.sellers, :quantity_produced_history); trimmed = 5))
    push!(ex1_v1_quantity_sold, trim_first.(getfield.(ex1_v1_sim.sellers, :quantity_sold_history); trimmed = 5))
    push!(ex1_v1_producer_surplus_singleton, trim_first.(calculate_profit_history.(ex1_v1_sim.sellers); trimmed = 5))
    push!(ex1_v1_reselling, trim_first.(getfield.(ex1_v1_sim.sellers, :reselling_history); trimmed = 5))
    push!(ex1_v1_buying_history, trim_first.(getfield.(ex1_v1_sim.buyers, :unit_buying_selling_history); trimmed = 5))
    push!(ex1_v1_quality_exp, trim_first.(mean([b.quality_expectation_history for b in ex1_v1_sim.buyers]); trimmed = 5))
    push!(ex1_v1_durability_exp, trim_first.(mean([b.durability_expectation_history for b in ex1_v1_sim.buyers]); trimmed = 5))
    push!(ex1_v1_quality_std, [std(mean.([getindex.(x,1) for x in getfield.(ex1_v1_sim.buyers, :quality_expectation_history)])), std(mean.([getindex.(x,2) for x in getfield.(ex1_v1_sim.buyers, :quality_expectation_history)]))])
    push!(ex1_v1_durability_std, [std(mean.([getindex.(x,1) for x in getfield.(ex1_v1_sim.buyers, :durability_expectation_history)])), std(mean.([getindex.(x,2) for x in getfield.(ex1_v1_sim.buyers, :durability_expectation_history)]))])

end


ex1_p1 = plot_ecdf(true, mean.(getindex.(ex1_v1_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, mean.(getindex.(ex1_v1_producer_surplus_singleton, 2)), "Producent bada oczekiwania konsumentów")

L1 = getindex.(ex1_v12_λ, 1) 
L2 = getindex.(ex1_v12_λ, 2)

hm_rs = [mean(sum.(getindex.(ex1_v1_producer_surplus_singleton, 2))[L2 .== l2]) for l2 in sort(unique(L2))]
hm_nrs = [mean(sum.(getindex.(ex1_v1_producer_surplus_singleton, 1))[L2 .== l2]) for l2 in sort(unique(L2))]

Plots.plot(sort(unique(L2)), hm_rs)
Plots.plot!(sort(unique(L2)), hm_nrs)



[count(L2 .== l2) for l2 in sort(unique(L2))]

55 .- vec([count((L1 .== l1) .& (L2 .== l2)) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]')

#### Total market surplus

hm_eq = [mean_trim_k(ex1_v1_total_surplus[(L1 .== l1) .& (L2 .== l2)];k=2) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'

ex4_p1 = StatsPlots.heatmap(sort(unique(L1)), sort(unique(L2)), hm_eq, xlabel = "λ, produkty oceniane osobiście", ylabel = "θ, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita", titlefontsize = 8)

using Plots
using LaTeXStrings
pyplot()

calculate_posterior(λ_i, λ_w; K1, K2) = K1*(1 - λ_i - λ_w) + K2 * (λ_i + λ_w)

li = LinRange(0.01:0.01:0.99)
lw = LinRange(0.01:0.01:0.99)

ac_p1 = contour(li, lw, (li,lw)->calculate_posterior(li,lw; K1=K1, K2=K2), levels=50, xlabel = "Szybkość uczenia się, własne doświadczenie, " * L"K_{jt}", ylabel = "Trwałość dobra, " * L"D_{jt}",title = "Całkowita użyteczność dla H = 1", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6)

Plots.savefig(ex4_p1, pwd() * "\\Plots\\ex4\\total surplus equal heatmap.svg")

#%% Producer 1 surplus, worse producer

mean_trim_k(x;k=0) = mean(sort(x)[(k+1):(end-k-1)])

hm_research = [mean_trim_k(sum.(getindex.(ex1_v1_producer_surplus_singleton, 2))[(L1 .== l1) .& (L2 .== l2)]; k = 1) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'

hm_no_research = [mean_trim_k(sum.(getindex.(ex1_v1_producer_surplus_singleton, 1))[(L1 .== l1) .& (L2 .== l2)]; k =1) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'

plot_min = min(minimum(hm_research), minimum(hm_no_research))
plot_max = max(maximum(hm_research), maximum(hm_no_research))

ex4_p3 = StatsPlots.heatmap(sort(unique(L1)), sort(unique(L2)), hm_research, xlabel = "λ, produkty oceniane osobiście", ylabel = "θ, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka producenta badającego konsumentów", titlefontsize = 8, clim = (plot_min, plot_max))

Plots.savefig(ex4_p3, pwd() * "\\plots\\ex4\\prod surplus research.svg")

ex4_p4 = StatsPlots.heatmap(sort(unique(L1)), sort(unique(L2)), hm_no_research, xlabel = "λ, produkty oceniane osobiście", ylabel = "θ, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka producenta niebadającego konsumentów", titlefontsize = 8, clim = (plot_min, plot_max))

Plots.savefig(ex4_p4, pwd() * "\\plots\\ex4\\prod surplus no research.svg")

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