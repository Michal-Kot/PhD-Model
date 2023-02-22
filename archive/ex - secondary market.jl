ex3_v200_total_surplus = []
ex3_v200_producer_surplus = []
ex3_v200_consumer_surplus = []
ex3_v200_price = []
ex3_v200_quantity_produced = []
ex3_v200_quantity_sold = []
ex3_v200_quantity_leased = []
ex3_v200_reselling = []
ex3_v200_producer_surplus_singleton = []
ex3_v200_buying_history = []
ex3_v200_quality = []
ex3_v200_durability = []
ex3_v200_margin = []
ex3_v200_quality_exp = []
ex3_v200_durability_exp = []
ex3_v200_H = []
ex3_v200_Li = []
ex3_v200_Lw = []
ex3_v200_sv = []

ex3_v201_total_surplus = []
ex3_v201_producer_surplus = []
ex3_v201_consumer_surplus = []
ex3_v201_price = []
ex3_v201_quantity_produced = []
ex3_v201_quantity_sold = []
ex3_v201_quantity_leased = []
ex3_v201_reselling = []
ex3_v201_producer_surplus_singleton = []
ex3_v201_buying_history = []
ex3_v201_quality = []
ex3_v201_durability = []
ex3_v201_margin = []
ex3_v201_quality_exp = []
ex3_v201_durability_exp = []
ex3_v201_sv = []

ex3_v202_total_surplus = []
ex3_v202_producer_surplus = []
ex3_v202_consumer_surplus = []
ex3_v202_price = []
ex3_v202_quantity_produced = []
ex3_v202_quantity_sold = []
ex3_v202_quantity_leased = []
ex3_v202_reselling = []
ex3_v202_producer_surplus_singleton = []
ex3_v202_buying_history = []
ex3_v202_quality = []
ex3_v202_durability = []
ex3_v202_margin = []
ex3_v202_quality_exp = []
ex3_v202_durability_exp = []
ex3_v202_sv = []

ex3_v203_total_surplus = []
ex3_v203_producer_surplus = []
ex3_v203_consumer_surplus = []
ex3_v203_price = []
ex3_v203_quantity_produced = []
ex3_v203_quantity_sold = []
ex3_v203_quantity_leased = []
ex3_v203_reselling = []
ex3_v203_producer_surplus_singleton = []
ex3_v203_buying_history = []
ex3_v203_quality = []
ex3_v203_durability = []
ex3_v203_margin = []
ex3_v203_quality_exp = []
ex3_v203_durability_exp = []
ex3_v203_sv = []

ex3_v200201_H = []
ex3_v200201_Li = []
ex3_v200201_Lw = []
ex3_v200201_sd = []
ex3_v200201_cu = []
ex3_v200201_sm = []
ex3_v200201_nl = []
ex3_v200201_ρm = []

for i in 1:500

    println(i)

    li = rand()
    lw = rand()
    h = sample(2:8)
    sd = sample(1:100000)
    cu = rand(Uniform(0.0,0.5))
    sm = sample(0.025:0.025:0.150)
    nl = sample(100:100:1200)
    ρm = rand(Uniform(0.5, 0.9))

    push!(ex3_v200201_H, h)
    push!(ex3_v200201_Li, li)
    push!(ex3_v200201_Lw, lw)
    push!(ex3_v200201_sd, sd)
    push!(ex3_v200201_cu, cu)
    push!(ex3_v200201_sm, sm)
    push!(ex3_v200201_nl, nl)
    push!(ex3_v200201_ρm, ρm)

    Random.seed!(sd)

    ex3_v200_sim = TO_GO(200, 2, 400, nl, [0.4, 0.4], [1., 1.], "random", li, lw, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], cu, true, 1, [ρm, 1.], "softmax", [true, false], [sm, 0.], h, false, false)

    push!(ex3_v200_total_surplus, calculate_surplus(ex3_v200_sim, "total", false))
    push!(ex3_v200_producer_surplus, calculate_surplus(ex3_v200_sim, "producer", false))
    push!(ex3_v200_consumer_surplus, calculate_surplus(ex3_v200_sim, "consumer,total",false))
    push!(ex3_v200_quality, getfield.(ex3_v200_sim.sellers, :quality_history))
    push!(ex3_v200_durability, getfield.(ex3_v200_sim.sellers, :durability_history))
    push!(ex3_v200_margin, getfield.(ex3_v200_sim.sellers, :margin_history))
    push!(ex3_v200_price, calculate_price_history.(ex3_v200_sim.sellers; product_life = h))
    push!(ex3_v200_quantity_produced, getfield.(ex3_v200_sim.sellers, :quantity_produced_history))
    push!(ex3_v200_quantity_sold, getfield.(ex3_v200_sim.sellers, :quantity_sold_history))
    push!(ex3_v200_producer_surplus_singleton, calculate_profit_history.(ex3_v200_sim.sellers))
    push!(ex3_v200_reselling, getfield.(ex3_v200_sim.sellers, :reselling_history))
    push!(ex3_v200_buying_history, getfield.(ex3_v200_sim.buyers, :unit_buying_selling_history))
    push!(ex3_v200_quality_exp, mean([b.quality_expectation_history for b in ex3_v200_sim.buyers]))
    push!(ex3_v200_durability_exp, mean([b.durability_expectation_history for b in ex3_v200_sim.buyers]))
    push!(ex3_v200_sv, getfield.(ex3_v200_sim.buyers, :signal_volume))

    Random.seed!(sd)

    ex3_v201_sim = TO_GO(200, 2, 400, nl, [0.4, 0.4], [1., 1.], "random", li, lw, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], cu, false, 1, [ρm, 1.], "softmax", [true, false], [sm, 0.], h, false, false)

    push!(ex3_v201_total_surplus, calculate_surplus(ex3_v201_sim, "total", false))
    push!(ex3_v201_producer_surplus, calculate_surplus(ex3_v201_sim, "producer", false))
    push!(ex3_v201_consumer_surplus, calculate_surplus(ex3_v201_sim, "consumer,total",false))
    push!(ex3_v201_quality, getfield.(ex3_v201_sim.sellers, :quality_history))
    push!(ex3_v201_durability, getfield.(ex3_v201_sim.sellers, :durability_history))
    push!(ex3_v201_margin, getfield.(ex3_v201_sim.sellers, :margin_history))
    push!(ex3_v201_price, calculate_price_history.(ex3_v201_sim.sellers; product_life = h))
    push!(ex3_v201_quantity_produced, getfield.(ex3_v201_sim.sellers, :quantity_produced_history))
    push!(ex3_v201_quantity_sold, getfield.(ex3_v201_sim.sellers, :quantity_sold_history))
    push!(ex3_v201_producer_surplus_singleton, calculate_profit_history.(ex3_v201_sim.sellers))
    push!(ex3_v201_reselling, getfield.(ex3_v201_sim.sellers, :reselling_history))
    push!(ex3_v201_buying_history, getfield.(ex3_v201_sim.buyers, :unit_buying_selling_history))
    push!(ex3_v201_quality_exp, mean([b.quality_expectation_history for b in ex3_v201_sim.buyers]))
    push!(ex3_v201_durability_exp, mean([b.durability_expectation_history for b in ex3_v201_sim.buyers]))
    push!(ex3_v201_sv, getfield.(ex3_v201_sim.buyers, :signal_volume))

    Random.seed!(sd)

    ex3_v202_sim = TO_GO(200, 2, 400, nl, [0.4, 0.4], [1., 1.], "random", li, lw, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], cu, true, 1, [ρm, 1.], "softmax", [false, false], [0., 0.], h, false, false)

    push!(ex3_v202_total_surplus, calculate_surplus(ex3_v202_sim, "total", false))
    push!(ex3_v202_producer_surplus, calculate_surplus(ex3_v202_sim, "producer", false))
    push!(ex3_v202_consumer_surplus, calculate_surplus(ex3_v202_sim, "consumer,total",false))
    push!(ex3_v202_quality, getfield.(ex3_v202_sim.sellers, :quality_history))
    push!(ex3_v202_durability, getfield.(ex3_v202_sim.sellers, :durability_history))
    push!(ex3_v202_margin, getfield.(ex3_v202_sim.sellers, :margin_history))
    push!(ex3_v202_price, calculate_price_history.(ex3_v202_sim.sellers; product_life = h))
    push!(ex3_v202_quantity_produced, getfield.(ex3_v202_sim.sellers, :quantity_produced_history))
    push!(ex3_v202_quantity_sold, getfield.(ex3_v202_sim.sellers, :quantity_sold_history))
    push!(ex3_v202_producer_surplus_singleton, calculate_profit_history.(ex3_v202_sim.sellers))
    push!(ex3_v202_reselling, getfield.(ex3_v202_sim.sellers, :reselling_history))
    push!(ex3_v202_buying_history, getfield.(ex3_v202_sim.buyers, :unit_buying_selling_history))
    push!(ex3_v202_quality_exp, mean([b.quality_expectation_history for b in ex3_v202_sim.buyers]))
    push!(ex3_v202_durability_exp, mean([b.durability_expectation_history for b in ex3_v202_sim.buyers]))
    push!(ex3_v202_sv, getfield.(ex3_v202_sim.buyers, :signal_volume))

    Random.seed!(sd)

    ex3_v203_sim = TO_GO(200, 2, 400, nl, [0.4, 0.4], [1., 1.], "random", li, lw, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], cu, false, 1, [ρm, 1.], "softmax", [false, false], [0., 0.], h, false, false)

    push!(ex3_v203_total_surplus, calculate_surplus(ex3_v203_sim, "total", false))
    push!(ex3_v203_producer_surplus, calculate_surplus(ex3_v203_sim, "producer", false))
    push!(ex3_v203_consumer_surplus, calculate_surplus(ex3_v203_sim, "consumer,total",false))
    push!(ex3_v203_quality, getfield.(ex3_v203_sim.sellers, :quality_history))
    push!(ex3_v203_durability, getfield.(ex3_v203_sim.sellers, :durability_history))
    push!(ex3_v203_margin, getfield.(ex3_v203_sim.sellers, :margin_history))
    push!(ex3_v203_price, calculate_price_history.(ex3_v203_sim.sellers; product_life = h))
    push!(ex3_v203_quantity_produced, getfield.(ex3_v203_sim.sellers, :quantity_produced_history))
    push!(ex3_v203_quantity_sold, getfield.(ex3_v203_sim.sellers, :quantity_sold_history))
    push!(ex3_v203_producer_surplus_singleton, calculate_profit_history.(ex3_v203_sim.sellers))
    push!(ex3_v203_reselling, getfield.(ex3_v203_sim.sellers, :reselling_history))
    push!(ex3_v203_buying_history, getfield.(ex3_v203_sim.buyers, :unit_buying_selling_history))
    push!(ex3_v203_quality_exp, mean([b.quality_expectation_history for b in ex3_v203_sim.buyers]))
    push!(ex3_v203_durability_exp, mean([b.durability_expectation_history for b in ex3_v203_sim.buyers]))
    push!(ex3_v203_sv, getfield.(ex3_v203_sim.buyers, :signal_volume))

end

plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), "Rynek wtórny TAK", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), "Rynek wtórny NIE")

plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), "Rynek wtórny TAK, badania TAK", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), "Rynek wtórny NIE, badania TAK")
plot_ecdf(false, sum.(getindex.(ex3_v202_producer_surplus_singleton, 1)), "Rynek wtórny TAK, badania NIE")
plot_ecdf(false, sum.(getindex.(ex3_v203_producer_surplus_singleton, 1)), "Rynek wtórny NIE, badania NIE")

median(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))) - median(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)))
median(sum.(getindex.(ex3_v202_producer_surplus_singleton, 1))) - median(sum.(getindex.(ex3_v203_producer_surplus_singleton, 1)))

mean(mean.(getindex.(ex3_v200_price, 2))) # Rynek wtórny TAK
mean(mean.(getindex.(ex3_v201_price, 2)))

mean(mean.(getindex.(ex3_v200_price, 1))) # Rynek wtórny TAK
mean(mean.(getindex.(ex3_v201_price, 1)))

mean(mean.(getindex.(ex3_v200_quantity_produced, 2))) # Rynek wtórny TAK
mean(mean.(getindex.(ex3_v201_quantity_produced, 2)))

mean(mean.(getindex.(ex3_v200_quantity_produced, 1))) # Rynek wtórny TAK
mean(mean.(getindex.(ex3_v201_quantity_produced, 1)))

mean(mean.(getindex.(ex3_v200_quantity_sold, 2))) # Rynek wtórny TAK
mean(mean.(getindex.(ex3_v201_quantity_sold, 2)))

mean(mean.(getindex.(ex3_v200_quantity_sold, 1))) # Rynek wtórny TAK
mean(mean.(getindex.(ex3_v201_quantity_sold, 1)))

plot_ecdf(true, mean.(getindex.(ex3_v200_quantity_produced, 1)), "Rynek wtórny TAK", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, mean.(getindex.(ex3_v201_quantity_produced, 1)), "Rynek wtórny NIE")

plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)[ex3_v200201_H .<= 5]), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)[ex3_v200201_H .<= 5]), "Producent bada oczekiwania konsumentów")

plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)[ex3_v200201_H .> 5]), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)[ex3_v200201_H .> 5]), "Producent bada oczekiwania konsumentów")

mdf = DataFrame(
    y = vcat(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))),
    h = Int.(repeat(ex3_v200201_H, outer = 2)),
    s = vcat(fill(1, lastindex(ex3_v200_producer_surplus_singleton)), fill(0, lastindex(ex3_v200_producer_surplus_singleton))),
    nl = Int.(repeat(ex3_v200201_nl, outer = 2))
)

model = GLM.lm(@formula(y~h+s+nl), mdf)
r2(model)

Plots.plot(mdf.y)
Plots.plot!(fitted(model))