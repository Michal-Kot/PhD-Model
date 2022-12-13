
# Experiment 4

include(pwd() * "\\methods\\methods.jl")

# Impact of secondary market on model dynamics.

ex4_v1_total_surplus = []
ex4_v1_producer_surplus = []
ex4_v1_consumer_surplus = []
ex4_v1_price = []
ex4_v1_quantity_produced = []
ex4_v1_quantity_sold = []
ex4_v1_quantity_leased = []
ex4_v1_reselling = []
ex4_v1_producer_surplus_singleton = []
ex4_v1_buying_history = []
ex4_v1_quality = []
ex4_v1_durability = []
ex4_v1_margin = []
ex4_v1_quality_exp = []
ex4_v1_durability_exp = []  

ex4_v2_total_surplus = []
ex4_v2_producer_surplus = []
ex4_v2_consumer_surplus = []
ex4_v2_price = []
ex4_v2_quantity_produced = []
ex4_v2_quantity_sold = []
ex4_v2_quantity_leased = []
ex4_v2_reselling = []
ex4_v2_producer_surplus_singleton = []
ex4_v2_buying_history = []
ex4_v2_quality = []
ex4_v2_durability = []
ex4_v2_margin = []
ex4_v2_quality_exp = []
ex4_v2_durability_exp = []  

ex4_v3_total_surplus = []
ex4_v3_producer_surplus = []
ex4_v3_consumer_surplus = []
ex4_v3_price = []
ex4_v3_quantity_produced = []
ex4_v3_quantity_sold = []
ex4_v3_quantity_leased = []
ex4_v3_reselling = []
ex4_v3_producer_surplus_singleton = []
ex4_v3_buying_history = []
ex4_v3_quality = []
ex4_v3_durability = []
ex4_v3_margin = []
ex4_v3_quality_exp = []
ex4_v3_durability_exp = []  

for i in 1:500

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    li = rand()
    lw = rand()
    h = sample(2:8)

    ex4_v1_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1., 1.], "random", li, lw, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[0.8, 2.], [0.8, 2.]], 0.1, true, 1, [0.7, 1.0], "softmax", [true, true], [0.1, 0.1], h, false)
    push!(ex4_v1_total_surplus, calculate_surplus(ex4_v1_sim, "total", true))
    push!(ex4_v1_producer_surplus, calculate_surplus(ex4_v1_sim, "producer", true))
    push!(ex4_v1_consumer_surplus, calculate_surplus(ex4_v1_sim, "consumer,total",true))
    push!(ex4_v1_price, calculate_price_history.(ex4_v1_sim.sellers; product_life = h))
    push!(ex4_v1_quantity_produced, getfield.(ex4_v1_sim.sellers, :quantity_produced_history))
    push!(ex4_v1_quantity_sold, getfield.(ex4_v1_sim.sellers, :quantity_sold_history))
    push!(ex4_v1_producer_surplus_singleton, calculate_profit_history.(ex4_v1_sim.sellers))
    push!(ex4_v1_reselling, getfield.(ex4_v1_sim.sellers, :reselling_history))
    push!(ex4_v1_buying_history, getfield.(ex4_v1_sim.buyers, :unit_buying_selling_history))
    push!(ex4_v1_quality, getfield.(ex4_v1_sim.sellers, :quality_history))
    push!(ex4_v1_durability, getfield.(ex4_v1_sim.sellers, :durability_history))
    push!(ex4_v1_margin, getfield.(ex4_v1_sim.sellers, :margin_history))
    push!(ex4_v1_quality_exp, mean([b.quality_expectation_history for b in ex4_v1_sim.buyers]))
    push!(ex4_v1_durability_exp, mean([b.durability_expectation_history for b in ex4_v1_sim.buyers]))

    ex4_v2_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1., 1.], "random", li, lw, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[0.8, 2.], [0.8, 2.]], 0.1, true, 1, [0.7, 1.0], "softmax", [false, true], [0, 0.1], h, false)

    push!(ex4_v2_total_surplus, calculate_surplus(ex4_v2_sim, "total", true))
    push!(ex4_v2_producer_surplus, calculate_surplus(ex4_v2_sim, "producer", true))
    push!(ex4_v2_consumer_surplus, calculate_surplus(ex4_v2_sim, "consumer,total",true))
    push!(ex4_v2_price, calculate_price_history.(ex4_v2_sim.sellers; product_life = 5))
    push!(ex4_v2_quantity_produced, getfield.(ex4_v2_sim.sellers, :quantity_produced_history))
    push!(ex4_v2_quantity_sold, getfield.(ex4_v2_sim.sellers, :quantity_sold_history))
    push!(ex4_v2_producer_surplus_singleton, calculate_profit_history.(ex4_v2_sim.sellers))
    push!(ex4_v2_reselling, getfield.(ex4_v2_sim.sellers, :reselling_history))
    push!(ex4_v2_buying_history, getfield.(ex4_v2_sim.buyers, :unit_buying_selling_history))
    push!(ex4_v2_quality, getfield.(ex4_v2_sim.sellers, :quality_history))
    push!(ex4_v2_durability, getfield.(ex4_v2_sim.sellers, :durability_history))
    push!(ex4_v2_margin, getfield.(ex4_v2_sim.sellers, :margin_history))
    push!(ex4_v2_quality_exp, mean([b.quality_expectation_history for b in ex4_v2_sim.buyers]))
    push!(ex4_v2_durability_exp, mean([b.durability_expectation_history for b in ex4_v2_sim.buyers]))

    ex4_v3_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1., 1.], "random", li, lw, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[0.8, 2.], [0.8, 2.]], 0.1, true, 1, [0.7, 1.0], "softmax", [false, false], [0., 0.], h, false)
    push!(ex4_v3_total_surplus, calculate_surplus(ex4_v3_sim, "total", true))
    push!(ex4_v3_producer_surplus, calculate_surplus(ex4_v3_sim, "producer", true))
    push!(ex4_v3_consumer_surplus, calculate_surplus(ex4_v3_sim, "consumer,total",true))
    push!(ex4_v3_price, calculate_price_history.(ex4_v3_sim.sellers; product_life = 5))
    push!(ex4_v3_quantity_produced, getfield.(ex4_v3_sim.sellers, :quantity_produced_history))
    push!(ex4_v3_quantity_sold, getfield.(ex4_v3_sim.sellers, :quantity_sold_history))
    push!(ex4_v3_producer_surplus_singleton, calculate_profit_history.(ex4_v3_sim.sellers))
    push!(ex4_v3_reselling, getfield.(ex4_v3_sim.sellers, :reselling_history))
    push!(ex4_v3_buying_history, getfield.(ex4_v3_sim.buyers, :unit_buying_selling_history))
    push!(ex4_v3_quality, getfield.(ex4_v3_sim.sellers, :quality_history))
    push!(ex4_v3_durability, getfield.(ex4_v3_sim.sellers, :durability_history))
    push!(ex4_v3_margin, getfield.(ex4_v3_sim.sellers, :margin_history))
    push!(ex4_v3_quality_exp, mean([b.quality_expectation_history for b in ex4_v3_sim.buyers]))
    push!(ex4_v3_durability_exp, mean([b.durability_expectation_history for b in ex4_v3_sim.buyers]))

end

function trim_outliers(x, p = 1)
    lower_bound = percentile(x, p)
    upper_bound = percentile(x,100-p)
    y = x[(x .>= lower_bound) .& (x .<= upper_bound)]
    return y
end

ex4_p1 = plot_ecdf(true, trim_outliers(ex4_v1_total_surplus,1), "Obaj gracze wykonują badania", xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false, trim_outliers(ex4_v2_total_surplus), "Jeden gracz wykonuje badania")
plot_ecdf(false, trim_outliers(ex4_v3_total_surplus), "Żaden gracz nie wykonuje badań")

Plots.savefig(ex4_p1, pwd() * "\\plots\\ex2\\total surplus research.svg")

# ... oraz nadwyżkę konsumenta ...

ex4_p2 = plot_ecdf(true, ex4_v1_consumer_surplus, "Obaj gracze wykonują badania", xlabel = "Nadwyżka konsumenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka konsumenta")
plot_ecdf(false, ex4_v2_consumer_surplus, "Jeden gracz wykonuje badania")
plot_ecdf(false, ex4_v3_consumer_surplus, "Żaden gracz nie wykonuje badań")

Plots.savefig(ex4_p2, pwd() * "\\plots\\ex2\\consumer surplus research.svg")

# ... a także nadwyżkę producenta ...

ex4_p3 = plot_ecdf(true, ex4_v1_producer_surplus, "Obaj gracze wykonują badania", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, ex4_v2_producer_surplus, "Jeden gracz wykonuje badania")
plot_ecdf(false, ex4_v3_producer_surplus, "Żaden gracz nie wykonuje badań")

Plots.savefig(ex4_p3, pwd() * "\\plots\\ex2\\producer surplus research.svg")

# Efekt istnienia rynku wtórnego nie jest istotny w przypadku rynku dóbr homogenicznych, ale jest istotny dla rynku, na którym występuje pionowe zróżnicowanie

ex4_p4 = plot_ecdf(true, mean.(getindex.(ex4_v1_producer_surplus_singleton,1)), "Obaj producenci badają rynek", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(getindex.(ex4_v1_producer_surplus_singleton,2)), "")

Plots.savefig(ex4_p4, pwd() * "\\plots\\ex2\\producer surplus research both.svg")

ex4_p5 = plot_ecdf(true, mean.(getindex.(ex4_v2_producer_surplus_singleton,1)), "Producent nie bada rynku", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(getindex.(ex4_v2_producer_surplus_singleton,2)), "Producent bada rynek")

Plots.savefig(ex4_p5, pwd() * "\\plots\\ex2\\producer surplus research one.svg")

ex4_p6 = plot_ecdf(true, mean.(getindex.(ex4_v3_producer_surplus_singleton,1)), "Żaden producent nie bada rynku", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(getindex.(ex4_v3_producer_surplus_singleton,2)), "",)

Plots.savefig(ex4_p6, pwd() * "\\plots\\ex2\\producer surplus research none.svg")

ex4_p7 = plot_ecdf(true, sum.(sum.(ex4_v1_quantity_sold)), "Obaj gracze wykonują badania", xlabel = "Liczba sprzedanych sztuk", ylabel = "F(x)", title = "Dystrybuanta empiryczna - popyt")
plot_ecdf(false,  sum.(sum.(ex4_v2_quantity_sold)), "Jeden gracz wykonuje badania")
plot_ecdf(false,  sum.(sum.(ex4_v3_quantity_sold)), "Żaden gracz nie wykonuje badań")

Plots.savefig(ex4_p7, pwd() * "\\plots\\ex2\\quantity sold.svg")

ex4_p8 = plot_ecdf(true, sum.(sum.(ex4_v1_quantity_produced)) .- sum.(sum.(ex4_v1_quantity_sold)), "Obaj gracze wykonują badania", xlabel = "Liczba niesprzedanych sztuk", ylabel = "F(x)", title = "Dystrybuanta empiryczna - liczba niesprzedanych sztuk")
plot_ecdf(false,  sum.(sum.(ex4_v2_quantity_produced)) .- sum.(sum.(ex4_v2_quantity_sold)), "Jeden gracz wykonuje badania")
plot_ecdf(false,  sum.(sum.(ex4_v3_quantity_produced)) .- sum.(sum.(ex4_v3_quantity_sold)), "Żaden gracz nie wykonuje badań")

Plots.savefig(ex4_p8, pwd() * "\\plots\\ex2\\quantity unsold.svg")