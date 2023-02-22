# Experiment 3

include(pwd() * "\\methods\\methods.jl")

# Impact on social communication of simulation dynamics

ex6_v1_total_surplus = []
ex6_v1_producer_surplus = []
ex6_v1_consumer_surplus = []
ex6_v1_price = []
ex6_v1_quantity_produced = []
ex6_v1_quantity_sold = []
ex6_v1_reselling = []
ex6_v1_producer_surplus_singleton = []
ex6_v1_buying_history = []
ex6_v1_quality = []
ex6_v1_durability = []
ex6_v1_margin = []

ex6_ss = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    ss = sample(LinRange(0.01:0.01:0.20))
    push!(ex6_ss, ss)

    ex6_v1_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1.1, 1.1], "random", rand(), rand(), "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[0.8, 2.], [0.8, 2.]], 0.1, true, 1, [0.7, 1.0], "softmax", [false, true], [0, ss], sample(3:7), false)

    push!(ex6_v1_total_surplus, calculate_surplus(ex6_v1_sim, "total", true))
    push!(ex6_v1_producer_surplus, calculate_surplus(ex6_v1_sim, "producer", true))
    push!(ex6_v1_consumer_surplus, calculate_surplus(ex6_v1_sim, "consumer,total",true))
    push!(ex6_v1_quality, getfield.(ex6_v1_sim.sellers, :quality_history))
    push!(ex6_v1_durability, getfield.(ex6_v1_sim.sellers, :durability_history))
    push!(ex6_v1_margin, getfield.(ex6_v1_sim.sellers, :margin_history))
    push!(ex6_v1_price, calculate_price_history.(ex6_v1_sim.sellers; product_life = 5))
    push!(ex6_v1_quantity_produced, getfield.(ex6_v1_sim.sellers, :quantity_produced_history))
    push!(ex6_v1_quantity_sold, getfield.(ex6_v1_sim.sellers, :quantity_sold_history))
    push!(ex6_v1_producer_surplus_singleton, calculate_profit_history.(ex6_v1_sim.sellers))
    push!(ex6_v1_reselling, getfield.(ex6_v1_sim.sellers, :reselling_history))
    push!(ex6_v1_buying_history, getfield.(ex6_v1_sim.buyers, :unit_buying_selling_history))

end


ex6_p1 = Plots.scatter(sort(unique(round.(ex6_ss, digits = 3))), [mean(sum.(getindex.(ex6_v1_producer_surplus_singleton, 2))[round.(ex6_ss, digits = 3) .== ss]) for ss in sort(unique(round.(ex6_ss, digits = 3)))], label = "Nadwyżka producenta", legend = :bottomright, title = "Nadwyżka producenta a wielkość próby", xlabel = "Wielkość próby badawczej jako % populacji", ylabel = "Zysk producenta")

ss = sort(unique(ex6_ss))

ss_df = DataFrame(y = [mean(sum.(getindex.(ex6_v1_producer_surplus_singleton, 2))[(ex6_ss .== s)]) for s in ss], xx = float.(ss), logx = log.(ss), sx = sqrt.(ss))

ss_model_l = fit(LinearModel, @formula(y~1+logx), ss_df)
adjr2(ss_model_l)
aic(ss_model_l)

ss_model_s = fit(LinearModel, @formula(y~1+sx), ss_df)
adjr2(ss_model_s)
aic(ss_model_s)

Plots.plot!(ss, fitted(ss_model), xlabel = "% populacji w badaniu", ylabel = "Zysk producenta wykorzystującego badanie", label = "Dopasowana funkcja logarytmiczna")

Plots.plot!(ss, fitted(ss_model) .- ss * 200 * 0.01 * 400, label = "Zysk producenta - koszt badania")

Plots.savefig(ex6_p1, pwd() * "\\Plots\\ex3\\research perc vs income.svg")

######

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

for i in 1:250

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    ss = sample(LinRange(0.05:0.05:0.20),2)
    push!(ex6_ss, ss)

    ex6_v1_sim = TO_GO(200, 2, 200, 300, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], 0.50, true, 1, [0.7, 1.], "softmax", [false, false], [0.0, 0.0], 5, false)

    push!(ex6_v1_total_surplus, calculate_surplus(ex6_v1_sim, "total", true))
    push!(ex6_v1_producer_surplus, calculate_surplus(ex6_v1_sim, "producer", true))
    push!(ex6_v1_consumer_surplus, calculate_surplus(ex6_v1_sim, "consumer,total",true))
    push!(ex6_v1_quality, getfield.(ex6_v1_sim.sellers, :quality_history))
    push!(ex6_v1_durability, getfield.(ex6_v1_sim.sellers, :durability_history))
    push!(ex6_v1_margin, getfield.(ex6_v1_sim.sellers, :margin_history))
    push!(ex6_v1_price, calculate_price_history.(ex6_v1_sim.sellers; product_life = 5))
    push!(ex6_v1_quantity_produced, getfield.(ex6_v1_sim.sellers, :quantity_produced_history))
    push!(ex6_v1_quantity_sold, getfield.(ex6_v1_sim.sellers, :quantity_sold_history))
    push!(ex6_v1_producer_surplus_singleton, calculate_profit_history.(ex6_v1_sim.sellers))
    push!(ex6_v1_reselling, getfield.(ex6_v1_sim.sellers, :reselling_history))
    push!(ex6_v1_buying_history, getfield.(ex6_v1_sim.buyers, :unit_buying_selling_history))

    ex6_v2_sim = TO_GO(200, 2, 200, 300, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], 0.50, true, 1, [0.7, 1.], "softmax", [true, true], ss, 5, false)

    push!(ex6_v2_total_surplus, calculate_surplus(ex6_v2_sim, "total", true))
    push!(ex6_v2_producer_surplus, calculate_surplus(ex6_v2_sim, "producer", true))
    push!(ex6_v2_consumer_surplus, calculate_surplus(ex6_v2_sim, "consumer,total",true))
    push!(ex6_v2_quality, getfield.(ex6_v2_sim.sellers, :quality_history))
    push!(ex6_v2_durability, getfield.(ex6_v2_sim.sellers, :durability_history))
    push!(ex6_v2_margin, getfield.(ex6_v2_sim.sellers, :margin_history))
    push!(ex6_v2_price, calculate_price_history.(ex6_v2_sim.sellers; product_life = 5))
    push!(ex6_v2_quantity_produced, getfield.(ex6_v2_sim.sellers, :quantity_produced_history))
    push!(ex6_v2_quantity_sold, getfield.(ex6_v2_sim.sellers, :quantity_sold_history))
    push!(ex6_v2_producer_surplus_singleton, calculate_profit_history.(ex6_v2_sim.sellers))
    push!(ex6_v2_reselling, getfield.(ex6_v2_sim.sellers, :reselling_history))
    push!(ex6_v2_buying_history, getfield.(ex6_v2_sim.buyers, :unit_buying_selling_history))

end

f(d;H) = (H-1)*d^(H) - H*d^(H-1) + 1

f(0.5; H=1)

RMSE(x,y) = sqrt(sum((x .- y) .^ 2))

Plots.scatter(getindex.(ex6_ss, 2), RMSE.(getindex.(ex6_vs2_quality, 1), getindex.(ex6_v2_quality, 2)))


Plots.scatter(mean.(getindex.(ex6_v1_quality, 1)), mean.(getindex.(ex6_v1_producer_surplus_singleton, 1)), smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0)
Plots.scatter!(mean.(getindex.(ex6_v1_quality, 2)), mean.(getindex.(ex6_v1_producer_surplus_singleton, 2)), smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0)
Plots.scatter!(mean.(getindex.(ex6_v2_quality, 1)), mean.(getindex.(ex6_v2_producer_surplus_singleton, 1)), smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0)
Plots.scatter!(mean.(getindex.(ex6_v2_quality, 2)), mean.(getindex.(ex6_v2_producer_surplus_singleton, 2)), smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0)

Plots.scatter(mean.(getindex.(ex6_v1_durability, 1)), mean.(getindex.(ex6_v1_producer_surplus_singleton, 1)), smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0)
Plots.scatter!(mean.(getindex.(ex6_v1_durability, 2)), mean.(getindex.(ex6_v1_producer_surplus_singleton, 2)), smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0)
Plots.scatter!(mean.(getindex.(ex6_v2_durability, 1)), mean.(getindex.(ex6_v2_producer_surplus_singleton, 1)), smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0)
Plots.scatter!(mean.(getindex.(ex6_v2_durability, 2)), mean.(getindex.(ex6_v2_producer_surplus_singleton, 2)), smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0)

ex2_p1 = Plots.scatter(mean.(getindex.(ex6_v1_price, 1)), mean.(getindex.(ex6_v1_producer_surplus_singleton, 1)), smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0, xlabel = "Marża %", ylabel = "Zysk firmy", label = "Producent nie bada konsumentów")
#Plots.scatter!(mean.(getindex.(ex6_v1_margin, 2))[mean.(getindex.(ex6_v1_margin, 2)) .>= 1.2], mean.(getindex.(ex6_v1_producer_surplus_singleton, 2))[mean.(getindex.(ex6_v1_margin, 2)) .>= 1.2], smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0)
Plots.scatter!(mean.(getindex.(ex6_v2_price, 1)), mean.(getindex.(ex6_v2_producer_surplus_singleton, 1)), smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0, label = "Producent bada konsumentów")
#Plots.scatter!(mean.(getindex.(ex6_v2_margin, 2))[mean.(getindex.(ex6_v2_margin, 2)) .>= 1.2], mean.(getindex.(ex6_v2_producer_surplus_singleton, 2))[mean.(getindex.(ex6_v2_margin, 2)) .>= 1.2], smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0)

Plots.savefig(ex2_p1, pwd() * "\\plots\\ex2\\Margin vs income research no research.svg")

StatsPlots.boxplot([mean.(getindex.(ex6_v2_margin, 1))[getindex.(ex6_ss, 1) .== ss] for ss in [0.05, 0.10, 0.15, 0.20]], markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0, label = "Producent bada konsumentów")

S1 = round.(getindex.(ex6_ss, 1), digits = 3)
S2 = round.(getindex.(ex6_ss, 2), digits = 3)

S1u = sort(unique(S1))
S2u = sort(unique(S2))

strategies = S1u

hm_eq = [mean(ex6_v2_total_surplus[(S1 .== s1) .& (S2 .== s2)]) for s1 in S1u, s2 in S2u]'

Plots.scatter(mean.(getindex.(ex6_v2_quality, 1)), mean.(getindex.(ex6_v2_producer_surplus_singleton, 1)))
Plots.scatter!(mean.(getindex.(ex6_v2_quality, 2)), mean.(getindex.(ex6_v2_producer_surplus_singleton, 2)))

ex6_p1 = StatsPlots.heatmap(S1u, S2u, hm_eq, xlabel = "% populacji, Producent 1", ylabel = "% populacji, Producent 2", title = "% badanej populacji a nadwyżka całkowita", titlefontsize = 8)

hm_p1 = [mean(sum.(trim_extremes.(getindex.(ex6_v2_producer_surplus_singleton, 1)[(S1 .== s1) .& (S2 .== s2)]))) for s1 in S1u, s2 in S2u]

hm_p1 = hm_p1 .- (hcat([fill(x,4) for x in [0.05, 0.10, 0.15, 0.20]]...) * 200 * 0.5)'

ex6_p1 = StatsPlots.heatmap(S1u, S2u, hm_p1', xlabel = "% populacji, Producent 1", ylabel = "% populacji, Producent 2", title = "% badanej populacji a nadwyżka  producenta 1", titlefontsize = 8)

hm_p2 = [mean(sum.(trim_extremes.(getindex.(ex6_v2_producer_surplus_singleton, 2)[(S1 .== s1) .& (S2 .== s2)]))) for s1 in S1u, s2 in S2u]

hm_p2 = hm_p2 .- (hcat([fill(x,4) for x in [0.05, 0.10, 0.15, 0.20]]...) * 200 * 0.5)

ex6_p2 = StatsPlots.heatmap(S1u, S2u, hm_p2', xlabel ="% populacji, Producent 1", ylabel = "% populacji, Producent 2", title = "% badanej populacji a nadwyżka  producenta 1", titlefontsize = 8)




CSV.write("C:/Users/User/Documents/test.csv", Tables.table(hm_p1), delim = ';', decimal = ',')
CSV.write("C:/Users/User/Documents/test2.csv", Tables.table(hm_p2), delim = ';', decimal = ',')
# producer 1 plays





function find_NE(x_init, y_init)

    chosen_strategies = []

    p1s = x_init
    p2s = y_init
    
    k = 1

    for i in 1:10

        if k == 1

            p1s = strategies[argmax(mean.([sum.(getindex.(ex6_v2_producer_surplus_singleton, 1))[(S1 .== s1) .& (S2 .== p2s)] for s1 in strategies]))]

        elseif k == 2

            p2s = strategies[argmax(mean.([sum.(getindex.(ex6_v2_producer_surplus_singleton, 2))[(S2 .== s2) .& (S1 .== p1s)] for s2 in strategies]))]

        end

        k = 3 - k

    end

    push!(chosen_strategies, [p1s, p2s])

    p1s = x_init
    p2s = y_init

    k = 2

    for i in 1:10

        if k == 1

            p1s = strategies[argmax(mean.([sum.(getindex.(ex6_v2_producer_surplus_singleton, 1))[(S1 .== s1) .& (S2 .== p2s)] for s1 in strategies]))]

        elseif k == 2

            p2s = strategies[argmax(mean.([sum.(getindex.(ex6_v2_producer_surplus_singleton, 2))[(S2 .== s2) .& (S1 .== p1s)] for s2 in strategies]))]

        end

        k = 3 - k

    end

    push!(chosen_strategies, [p1s, p2s])

    return(chosen_strategies)

end

NE_sim = [find_NE(x,y) for x in strategies, y in strategies]

[sum(x .^ 2) for x in getindex.(NE_sim, 1)]

Plots.heatmap([sum(x .^ 2) for x in getindex.(NE_sim, 1)])
Plots.heatmap([sum(x .^ 2) for x in getindex.(NE_sim, 2)])



x = getindex.(chosen_strategies, 1)
y = getindex.(chosen_strategies, 2)
d_x = diff(x)
push!(d_x, 0)
d_y = diff(y)
push!(d_y,0)

quiver(x,y,quiver=(d_x, d_y), color = :green, xticks = (round.(LinRange(0.0:0.025:0.15), digits=3), round.(LinRange(0.0:0.025:0.15), digits=3)), yticks = (round.(LinRange(0.0:0.025:0.15), digits=3), round.(LinRange(0.0:0.025:0.15), digits=3)), xlim = (0, 0.15), ylim = (0, 0.15))
Plots.scatter!(x,y, xlabel = "Strategia, Firma 1", ylabel = "Strategia, Firma 2", markershape = :none, markercolor = :white, markerstrokecolor = :white, markersize = 0)






p2s = argmax(mean.([sum.(getindex.(ex6_v2_producer_surplus_singleton, 2))[(S2 .== s2) .& (S1 .== p1s)] for s2 in sort(unique(S2))]))

p1s = argmax(mean.([sum.(getindex.(ex6_v2_producer_surplus_singleton, 1))[(S1 .== s1) .& (S2 .== p2s)] for s1 in sort(unique(S1))]))

p2s = argmax(mean.([sum.(getindex.(ex6_v2_producer_surplus_singleton, 2))[(S2 .== s2) .& (S1 .== p1s)] for s2 in sort(unique(S2))]))

p1s = argmax(mean.([sum.(getindex.(ex6_v2_producer_surplus_singleton, 1))[(S1 .== s1) .& (S2 .== p2s)] for s1 in sort(unique(S1))]))

p2s = argmax(mean.([sum.(getindex.(ex6_v2_producer_surplus_singleton, 2))[(S2 .== s2) .& (S1 .== p1s)] for s2 in sort(unique(S2))]))

Plots.scatter(sort(unique(ex6_ss)), [mean(mean.(getindex.(ex6_v1_producer_surplus_singleton, 2))[ex6_ss .== ss]) for ss in sort(unique(ex6_ss))], smooth = true)

Plots.scatter(ex6_nn, mean.(getindex.(ex6_v1_producer_surplus_singleton, 1)), smooth = true)
Plots.scatter(ex6_nn, mean.(getindex.(ex6_v1_producer_surplus_singleton, 2)), smooth = true)

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

discount_min = []

for i in 1:100

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    dm = rand(Uniform(0.5,0.95))
    push!(discount_min, dm)

    ex6_v2_sim = TO_GO(200, 2, 200, 300, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.5, 0.95], [0.5, 0.95]], [[.8, 2.], [.8, 2.]], 0.10, true, true, 1, [dm, 1.], "softmax", [true, false], [0.1, 0.0])

    push!(ex6_v2_total_surplus, calculate_surplus(ex6_v2_sim, "total", true))
    push!(ex6_v2_producer_surplus, calculate_surplus(ex6_v2_sim, "producer", true))
    push!(ex6_v2_consumer_surplus, calculate_surplus(ex6_v2_sim, "consumer,total",true))
    push!(ex6_v2_quality, getfield.(ex6_v2_sim.sellers, :quality_history))
    push!(ex6_v2_durability, getfield.(ex6_v2_sim.sellers, :durability_history))
    push!(ex6_v2_margin, getfield.(ex6_v2_sim.sellers, :margin_history))
    push!(ex6_v2_price, calculate_price_history.(ex6_v2_sim.sellers))
    push!(ex6_v2_quantity_produced, getfield.(ex6_v2_sim.sellers, :quantity_produced_history))
    push!(ex6_v2_quantity_sold, getfield.(ex6_v2_sim.sellers, :quantity_sold_history))
    push!(ex6_v2_quantity_leased, getfield.(ex6_v2_sim.sellers, :quantity_leased_history))
    push!(ex6_v2_producer_surplus_singleton, calculate_profit_history.(ex6_v2_sim.sellers))
    push!(ex6_v2_reselling, getfield.(ex6_v2_sim.sellers, :reselling_history))
    push!(ex6_v2_buying_history, getfield.(ex6_v2_sim.buyers, :unit_buying_selling_history))

end

Plots.scatter(discount_min,  mean.(getindex.(ex6_v2_quantity_sold, 1)) ./ (mean.(getindex.(ex6_v2_quantity_sold, 1)) .+  mean.(getindex.(ex6_v2_quantity_leased, 1))), smooth = true)
Plots.scatter(discount_min, mean.(getindex.(ex6_v2_quantity_sold, 2)) ./ (mean.(getindex.(ex6_v2_quantity_sold, 2)) .+  mean.(getindex.(ex6_v2_quantity_leased, 2))), smooth = true)

#####



ex6_v3_total_surplus = []
ex6_v3_producer_surplus = []
ex6_v3_consumer_surplus = []
ex6_v3_price = []
ex6_v3_quantity_produced = []
ex6_v3_quantity_sold = []
ex6_v3_quantity_leased = []
ex6_v3_reselling = []
ex6_v3_producer_surplus_singleton = []
ex6_v3_buying_history = []
ex6_v3_quality = []
ex6_v3_durability = []
ex6_v3_margin = []
ex6_v3_common = []

ex6_v4_total_surplus = []
ex6_v4_producer_surplus = []
ex6_v4_consumer_surplus = []
ex6_v4_price = []
ex6_v4_quantity_produced = []
ex6_v4_quantity_sold = []
ex6_v4_quantity_leased = []
ex6_v4_reselling = []
ex6_v4_producer_surplus_singleton = []
ex6_v4_buying_history = []
ex6_v4_quality = []
ex6_v4_durability = []
ex6_v4_margin = []

ex6_v4_ss = []
ex6_v3_H = []
ex6_v3_Li = []
ex6_v3_Lw = []

for i in 1:250

    if (mod(i,10) == 0) | (i == 1) | (i == 2)
        println(i)
    end

    ss = sample(LinRange(0.025:0.025:0.10))
    push!(ex6_v4_ss, ss)

    li = rand()
    lw = rand()
    h = sample(2:8)

    push!(ex6_v3_H, h)
    push!(ex6_v3_Li, li)
    push!(ex6_v3_Lw, lw)

    ex6_v3_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1., 1.], "random", li, lw, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], 0.25, true, 1, [0.7, 1.], "softmax", [false, true], [0., ss], h, true)

    push!(ex6_v3_total_surplus, calculate_surplus(ex6_v3_sim, "total", true))
    push!(ex6_v3_producer_surplus, calculate_surplus(ex6_v3_sim, "producer", true))
    push!(ex6_v3_consumer_surplus, calculate_surplus(ex6_v3_sim, "consumer,total",true))
    push!(ex6_v3_quality, getfield.(ex6_v3_sim.sellers, :quality_history))
    push!(ex6_v3_durability, getfield.(ex6_v3_sim.sellers, :durability_history))
    push!(ex6_v3_margin, getfield.(ex6_v3_sim.sellers, :margin_history))
    push!(ex6_v3_price, calculate_price_history.(ex6_v3_sim.sellers; product_life = h))
    push!(ex6_v3_quantity_produced, getfield.(ex6_v3_sim.sellers, :quantity_produced_history))
    push!(ex6_v3_quantity_sold, getfield.(ex6_v3_sim.sellers, :quantity_sold_history))
    push!(ex6_v3_producer_surplus_singleton, calculate_profit_history.(ex6_v3_sim.sellers))
    push!(ex6_v3_reselling, getfield.(ex6_v3_sim.sellers, :reselling_history))
    push!(ex6_v3_buying_history, getfield.(ex6_v3_sim.buyers, :unit_buying_selling_history))

    ex6_v4_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1., 1.], "random", li, lw, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], 0.25, true, 1, [0.7, 1.], "softmax", [false, true], [0., ss], h, false)

    push!(ex6_v4_total_surplus, calculate_surplus(ex6_v4_sim, "total", true))
    push!(ex6_v4_producer_surplus, calculate_surplus(ex6_v4_sim, "producer", true))
    push!(ex6_v4_consumer_surplus, calculate_surplus(ex6_v4_sim, "consumer,total",true))
    push!(ex6_v4_quality, getfield.(ex6_v4_sim.sellers, :quality_history))
    push!(ex6_v4_durability, getfield.(ex6_v4_sim.sellers, :durability_history))
    push!(ex6_v4_margin, getfield.(ex6_v4_sim.sellers, :margin_history))
    push!(ex6_v4_price, calculate_price_history.(ex6_v4_sim.sellers; product_life = h))
    push!(ex6_v4_quantity_produced, getfield.(ex6_v4_sim.sellers, :quantity_produced_history))
    push!(ex6_v4_quantity_sold, getfield.(ex6_v4_sim.sellers, :quantity_sold_history))
    push!(ex6_v4_producer_surplus_singleton, calculate_profit_history.(ex6_v4_sim.sellers))
    push!(ex6_v4_reselling, getfield.(ex6_v4_sim.sellers, :reselling_history))
    push!(ex6_v4_buying_history, getfield.(ex6_v4_sim.buyers, :unit_buying_selling_history))

end


plot_ecdf(true, mean.(getindex.(ex6_v3_producer_surplus_singleton, 1)), "RW tak, WJ nie, CR nie"; xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Trwałość")
plot_ecdf(false, mean.(getindex.(ex6_v3_producer_surplus_singleton, 2)), "RW tak, WJ nie, CR tak")
plot_ecdf(false, mean.(getindex.(ex6_v4_producer_surplus_singleton, 1)), "RW nie, WJ nie, CR nie")
plot_ecdf(false, mean.(getindex.(ex6_v4_producer_surplus_singleton, 2)), "RW nie, WJ nie, CR tak")


Plots.scatter(getindex.(ex6_v4_ss, 1), ex6_v4_producer_surplus_singleton)

function mean_per(x_mean, x_per)

    return [mean(x_mean[x_per .== xp]) for xp in sort(unique(x_per))]

end

# zmienna grupa

function get_profit(mss, css; y, ss1, ss2)
    return mean(y[(ss1 .== mss) .& (ss2 .== css)])
end

MSS = sort(unique(getindex.(ex6_v4_ss, 1)))
CSS = sort(unique(getindex.(ex6_v4_ss, 2)))

Plots.contour(CSS, MSS, (CSS, MSS) -> get_profit(CSS, MSS; y=sum.(getindex.(ex6_v4_producer_surplus_singleton, 2)), ss1 = getindex.(ex6_v4_ss, 1), ss2 = getindex.(ex6_v4_ss, 2)), levels=4, xlabel = "Wstępna jakość dobra, " * L"K_{jt}", ylabel = "Trwałość dobra, " * L"D_{jt}",title = "Całkowita użyteczność dla H = 1", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6)

ex2_p5 =Plots.scatter(sort(unique(round.(getindex.(ex6_v4_ss, 2) .- getindex.(ex6_v4_ss, 1), digits = 2))), mean_per(sum.(getindex.(ex6_v4_producer_surplus_singleton, 2)) - 400 * 0.01 * getindex.(ex6_v4_ss, 2) * 200, round.(getindex.(ex6_v4_ss, 2) .- getindex.(ex6_v4_ss, 1), digits = 2)), xlabel = "Różnica w wielkościach grup badanych konsumentów", ylabel = "Zysk", title = "Zysk a wielkość grupy badawczej", legend = nothing)

Plots.scatter!(sort(unique(round.(getindex.(ex6_v4_ss, 1) .- getindex.(ex6_v4_ss, 2), digits = 2))), mean_per(sum.(getindex.(ex6_v4_producer_surplus_singleton, 1))- 400 * 0.01 * getindex.(ex6_v4_ss, 1) * 200, round.(getindex.(ex6_v4_ss, 1) .- getindex.(ex6_v4_ss, 2), digits = 2)), xlabel = "Różnica w wielkościach grup badanych konsumentów", ylabel = "Zysk", title = "Zysk a wielkość grupy badawczej", legend = nothing)

 mean_per(sum.(getindex.(ex6_v3_producer_surplus_singleton, 2)) , round.(getindex.(ex6_v4_ss, 2) .- getindex.(ex6_v4_ss, 1), digits = 2))

 Plots.histogram(sum.(getindex.(ex6_v3_producer_surplus_singleton, 2)[round.(getindex.(ex6_v4_ss, 2) .- getindex.(ex6_v4_ss, 1), digits = 2) .== -0.15]), alpha = 0.25, bins = 0:100:1500, normalize=true)
 Plots.histogram!(sum.(getindex.(ex6_v3_producer_surplus_singleton, 2)[round.(getindex.(ex6_v4_ss, 2) .- getindex.(ex6_v4_ss, 1), digits = 2) .== 0.0]), alpha = 0.25, bins = 0:100:1500, normalize=true)
 Plots.histogram!(sum.(getindex.(ex6_v3_producer_surplus_singleton, 2)[round.(getindex.(ex6_v4_ss, 2) .- getindex.(ex6_v4_ss, 1), digits = 2) .== 0.15]), alpha = 0.25, bins = 0:100:1500, normalize=true)

 # stała grupa

 ex2_p5 =Plots.scatter!(sort(unique(round.(getindex.(ex6_v4_ss, 2) .- getindex.(ex6_v4_ss, 1), digits = 2))), mean_per(sum.(getindex.(ex6_v3_producer_surplus_singleton, 2))- 400 * 0.01 * getindex.(ex6_v4_ss, 2) * 200 , round.(getindex.(ex6_v4_ss, 2) .- getindex.(ex6_v4_ss, 1), digits = 2)), xlabel = "Różnica w wielkościach grup badanych konsumentów", ylabel = "Zysk", title = "Zysk a wielkość grupy badawczej", legend = nothing)

 Plots.scatter!(sort(unique(round.(getindex.(ex6_v4_ss, 1) .- getindex.(ex6_v4_ss, 2), digits = 2))), mean_per(sum.(getindex.(ex6_v3_producer_surplus_singleton, 1))- 400 * 0.01 * getindex.(ex6_v4_ss, 1) * 200, round.(getindex.(ex6_v4_ss, 1) .- getindex.(ex6_v4_ss, 2), digits = 2)), xlabel = "Różnica w wielkościach grup badanych konsumentów", ylabel = "Zysk", title = "Zysk a wielkość grupy badawczej", legend = nothing)

add_smoothing_spline(round.(getindex.(ex6_v4_ss, 2) .- getindex.(ex6_v4_ss, 1), digits = 2), sum.(getindex.(ex6_v4_producer_surplus_singleton, 2)) - 0.05 * getindex.(ex6_v4_ss, 2) * 200, "blue", "", 0.0001)

Plots.savefig(ex2_p5, pwd() * "\\plots\\ex2\\profit vs diff in research gorup size.svg")

function add_smoothing_spline(x,y,clr,lbl,λ=0.05)
    spl = fit(SmoothingSpline, x, y, λ)
    y_hat = predict(spl)
    plot!(sort(x), y_hat[sortperm(x)], label = lbl, color = clr)
end

ex5_p1 = plot_ecdf(true, sum.(getindex.(ex6_v3_producer_surplus_singleton, 1))[getindex.(ex6_v4_ss, 1) .<= 0.05], "Stała grupa badanych, <= 5%"; xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Nadwyżka producenta")
plot_ecdf(false, sum.(getindex.(ex6_v3_producer_surplus_singleton, 1))[getindex.(ex6_v4_ss, 1) .>= 0.20], "Stała grupa badanych, > 10%")
plot_ecdf(false, sum.(getindex.(ex6_v4_producer_surplus_singleton, 1))[getindex.(ex6_v4_ss, 1) .<= 0.10], "Zmienna grupa badanych, <= 10%")
plot_ecdf(false, sum.(getindex.(ex6_v4_producer_surplus_singleton, 1))[getindex.(ex6_v4_ss, 1) .> 0.10], "Zmienna grupa badanych, > 10%")

Plots.scatter(sort(unique(getindex.(ex6_v4_ss, 1))), [mean(sum.(getindex.(ex6_v3_producer_surplus_singleton, 1))[getindex.(ex6_v4_ss, 1) .== x]) for x in sort(unique(getindex.(ex6_v4_ss, 1)))])

Plots.scatter(sort(unique(getindex.(ex6_v4_ss, 2))), [mean(sum.(getindex.(ex6_v3_producer_surplus_singleton, 2))[getindex.(ex6_v4_ss, 2) .== x]) for x in sort(unique(getindex.(ex6_v4_ss, 2)))])

UnequalVarianceTTest(sum.(getindex.(ex6_v3_producer_surplus_singleton, 1)), sum.(getindex.(ex6_v4_producer_surplus_singleton, 1)))
UnequalVarianceTTest(sum.(getindex.(ex6_v3_producer_surplus_singleton, 1))[getindex.(ex6_v4_ss, 1) .<= 0.10], sum.(getindex.(ex6_v4_producer_surplus_singleton, 1))[getindex.(ex6_v4_ss, 1) .<= 0.10])
UnequalVarianceTTest(sum.(getindex.(ex6_v3_producer_surplus_singleton, 1))[getindex.(ex6_v4_ss, 1) .> 0.10], sum.(getindex.(ex6_v4_producer_surplus_singleton, 1))[getindex.(ex6_v4_ss, 1) .> 0.10])

Plots.scatter(sort(unique(ex6_v3_common)), [mean(sum.(getindex.(ex6_v3_producer_surplus_singleton, 1))[ex6_v3_common .== x]) for x in sort(unique(ex6_v3_common))])
Plots.scatter!(ex6_v4_common, ex6_v4_producer_surplus)

ex5_p1 = plot_ecdf(true, ex6_v3_total_surplus[[any(x .> [0.10, 0.10]) for x in ex6_v4_ss]], "Stała grupa badanych, > 10%"; xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Trwałość")
plot_ecdf(false, ex6_v4_total_surplus[[all(x .> [0.10, 0.10]) for x in ex6_v4_ss]], "Zmienna grupa badanych, > 10%")

col_dic = Dict([0.05, 0.10, 0.15, 0.20] .=> ["red", "green", "blue", "yellow"])





ex5_p1 = plot_ecdf(true, ex6_v3_total_surplus, "Rynek wtórny istnieje, producent bada"; xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Trwałość")
plot_ecdf(false, ex6_v4_total_surplus, "Rynek wtórny istnieje, producent nie bada")

S1 = round.(getindex.(ex6_v4_ss, 1), digits = 2)
S2 = round.(getindex.(ex6_v4_ss, 2), digits = 2)

[count((S1 .== s1) .& (S2 .== s2)) for s1 in sort(unique(S1)), s2 in sort(unique(S2))]

hm_p1 = [mean(sum.(getindex.(ex6_v3_producer_surplus_singleton, 1)[(S1 .== s1) .& (S2 .== s2)])) for s1 in sort(unique(S1)), s2 in sort(unique(S2))]

hm_p2 = [mean(sum.(getindex.(ex6_v3_producer_surplus_singleton, 2)[(S1 .== s1) .& (S2 .== s2)])) for s1 in sort(unique(S1)), s2 in sort(unique(S2))]


opt_p1 = []

k = 1
for i in eachcol(hm_p1)
    push!(opt_p1, CartesianIndex(argmax(i),k))
    k += 1
end

opt_p2 = []

k = 1
for j in eachrow(hm_p2)
    push!(opt_p2, CartesianIndex(k, argmax(j)))
    k += 1
end

opt_p1
opt_p2

NE = opt_p1[findall(x->x in opt_p2, opt_p1)]
NE = Tuple.(NE)
NE = [(sort(unique(S1))[ne[1]], sort(unique(S1))[ne[2]]) for ne in NE]