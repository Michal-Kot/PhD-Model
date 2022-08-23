# Experiment 1

include(pwd() * "\\methods\\methods.jl")

# Impact on social communication of simulation dynamics

ex1_v1_total_surplus = []
ex1_v1_producer_surplus = []
ex1_v1_consumer_surplus = []
ex1_v1_price = []
ex1_v1_quantity_produced = []
ex1_v1_quantity_sold = []

ex1_v2_total_surplus = []
ex1_v2_producer_surplus = []
ex1_v2_consumer_surplus = []
ex1_v2_price = []
ex1_v2_quantity_produced = []
ex1_v2_quantity_sold = []

ex1_v12_λ = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    λ_ind = rand()
    λ_wom = rand()

    push!(ex1_v12_λ, [λ_ind, λ_wom])

    ex1_v1_sim = TO_GO(150, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.3, 0.3], "random", λ_ind, λ_wom, "stochastic", ["stochastic", "stochastic"],[[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]],[[0.,2.], [0.,2.]])
    push!(ex1_v1_total_surplus, calculate_surplus(ex1_v1_sim, "total", true))
    push!(ex1_v1_producer_surplus, calculate_surplus(ex1_v1_sim, "producer", true))
    push!(ex1_v1_consumer_surplus, calculate_surplus(ex1_v1_sim, "consumer,total",true))
    push!(ex1_v1_price, calculate_price_history.(ex1_v1_sim.sellers))
    push!(ex1_v1_quantity_produced, getfield.(ex1_v1_sim.sellers, :quantity_produced_history))
    push!(ex1_v1_quantity_sold, getfield.(ex1_v1_sim.sellers, :quantity_sold_history))

    ex1_v2_sim = TO_GO(150, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.3, 0.3], "random", λ_ind, λ_wom, "stochastic", ["stochastic", "stochastic"],[[0.9, 1.5], [0.6, 1.1]], [[0.2, 0.5], [0.3, 0.6]],[[0.,2.], [0.,2.]])
    push!(ex1_v2_total_surplus, calculate_surplus(ex1_v2_sim, "total", true))
    push!(ex1_v2_producer_surplus, calculate_surplus(ex1_v2_sim, "producer", true))
    push!(ex1_v2_consumer_surplus, calculate_surplus(ex1_v2_sim, "consumer,total",true))
    push!(ex1_v2_price, calculate_price_history.(ex1_v2_sim.sellers))
    push!(ex1_v2_quantity_produced, getfield.(ex1_v2_sim.sellers, :quantity_produced_history))
    push!(ex1_v2_quantity_sold, getfield.(ex1_v2_sim.sellers, :quantity_sold_history))

end

L1 = getindex.(ex1_v12_λ, 1) 
L2 = getindex.(ex1_v12_λ, 2)

L1, L1u = cut_integer(L1, 20)
L2, L2u = cut_integer(L2, 20)


ex1_v12_λ[562]

hm_eq = [mean_c(ex1_v1_total_surplus[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'
hm_dq = [mean_c(ex1_v2_total_surplus[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'

max_plot = maximum([maximum(hm_eq), maximum(hm_dq)])
min_plot = minimum([minimum(hm_eq), minimum(hm_dq)])

ex4_p1 = heatmap(L1u, L2u, hm_eq, xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość początkowa", titlefontsize = 8)

savefig(ex4_p1, pwd() * "\\plots\\ex4_equal quality heatmap lambdas.svg")

ex4_p2 = heatmap(L1u, L2u, hm_dq, xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, różna jakość początkowa", titlefontsize = 8)

savefig(ex4_p2, pwd() * "\\plots\\ex4_different quality heatmap lambdas.svg")

#

plot(ex1_v1_quantity_produced[1])
plot!(ex1_v1_quantity_sold[1])

sum(sum(ex1_v1_quantity_produced[1])[10:end]) - sum(sum(ex1_v1_quantity_sold[1])[10:end])

scatter(getindex.(ex1_v12_λ,2),sum.(sum.(ex1_v1_quantity_sold)))
scatter!(getindex.(ex1_v12_λ,2),sum.(sum.(ex1_v2_quantity_sold)))



ex1_v1_sim = TO_GO(150, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.3, 0.3], "random", 0.25, 0.25, "stochastic")