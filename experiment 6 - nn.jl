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

ex6_ss = []

for i in 1:100

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    ss = sample([0.05, 0.125, 0.175, 0.200])
    push!(ex6_ss, ss)

    ex6_v1_sim = TO_GO(200, 2, 200, 300, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.5, 0.95], [0.5, 0.95]], [[.8, 2.], [.8, 2.]], 0.10, true, true, 1, [0.7, 1.], "softmax", [true, false], [ss, 0.0], 5)

    push!(ex6_v1_total_surplus, calculate_surplus(ex6_v1_sim, "total", true))
    push!(ex6_v1_producer_surplus, calculate_surplus(ex6_v1_sim, "producer", true))
    push!(ex6_v1_consumer_surplus, calculate_surplus(ex6_v1_sim, "consumer,total",true))
    push!(ex6_v1_quality, getfield.(ex6_v1_sim.sellers, :quality_history))
    push!(ex6_v1_durability, getfield.(ex6_v1_sim.sellers, :durability_history))
    push!(ex6_v1_margin, getfield.(ex6_v1_sim.sellers, :margin_history))
    push!(ex6_v1_price, calculate_price_history.(ex6_v1_sim.sellers; product_life = 5))
    push!(ex6_v1_quantity_produced, getfield.(ex6_v1_sim.sellers, :quantity_produced_history))
    push!(ex6_v1_quantity_sold, getfield.(ex6_v1_sim.sellers, :quantity_sold_history))
    push!(ex6_v1_quantity_leased, getfield.(ex6_v1_sim.sellers, :quantity_leased_history))
    push!(ex6_v1_producer_surplus_singleton, calculate_profit_history.(ex6_v1_sim.sellers))
    push!(ex6_v1_reselling, getfield.(ex6_v1_sim.sellers, :reselling_history))
    push!(ex6_v1_buying_history, getfield.(ex6_v1_sim.buyers, :unit_buying_selling_history))

end


ex6_p1 = Plots.scatter(sort(unique(round.(ex6_ss, digits = 3))), [mean(sum.(getindex.(ex6_v1_producer_surplus_singleton, 1))[round.(ex6_ss, digits = 3) .== ss]) for ss in sort(unique(round.(ex6_ss, digits = 3)))], label = "Nadwyżka producenta", legend = :topleft, title = "Nadwyżka producenta a wielkość próby", xlabel = "Wielkość próby badawczej jako % populacji", ylabel = "Zysk producenta", xlim = (0, 0.16))

ss = sort(unique(ex6_ss[round.(ex6_ss, digits = 3) .<= 0.15]))

ss_df = DataFrame(y = [mean(sum.(getindex.(ex6_v1_producer_surplus_singleton, 1))[(ex6_ss .== s)]) for s in ss], logx = log.(ss))

ss_model = GLM.lm(@formula(y~logx), ss_df)

Plots.plot!(ss, fitted(ss_model), xlabel = "% populacji w badaniu", ylabel = "Zysk producenta wykorzystującego badanie", label = "Dopasowana funkcja logarytmiczna")

Plots.plot!(ss, fitted(ss_model) .- ss * 200 * 0.8, label = "Zysk producenta - koszt badania")

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

    ex6_v1_sim = TO_GO(200, 2, 200, 300, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], 0.50, true, true, 1, [0.7, 1.], "softmax", [false, false], [0.0, 0.0], 5)

    push!(ex6_v1_total_surplus, calculate_surplus(ex6_v1_sim, "total", true))
    push!(ex6_v1_producer_surplus, calculate_surplus(ex6_v1_sim, "producer", true))
    push!(ex6_v1_consumer_surplus, calculate_surplus(ex6_v1_sim, "consumer,total",true))
    push!(ex6_v1_quality, getfield.(ex6_v1_sim.sellers, :quality_history))
    push!(ex6_v1_durability, getfield.(ex6_v1_sim.sellers, :durability_history))
    push!(ex6_v1_margin, getfield.(ex6_v1_sim.sellers, :margin_history))
    push!(ex6_v1_price, calculate_price_history.(ex6_v1_sim.sellers; product_life = 5))
    push!(ex6_v1_quantity_produced, getfield.(ex6_v1_sim.sellers, :quantity_produced_history))
    push!(ex6_v1_quantity_sold, getfield.(ex6_v1_sim.sellers, :quantity_sold_history))
    push!(ex6_v1_quantity_leased, getfield.(ex6_v1_sim.sellers, :quantity_leased_history))
    push!(ex6_v1_producer_surplus_singleton, calculate_profit_history.(ex6_v1_sim.sellers))
    push!(ex6_v1_reselling, getfield.(ex6_v1_sim.sellers, :reselling_history))
    push!(ex6_v1_buying_history, getfield.(ex6_v1_sim.buyers, :unit_buying_selling_history))

    ex6_v2_sim = TO_GO(200, 2, 200, 300, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], 0.50, true, true, 1, [0.7, 1.], "softmax", [true, true], ss, 5)

    push!(ex6_v2_total_surplus, calculate_surplus(ex6_v2_sim, "total", true))
    push!(ex6_v2_producer_surplus, calculate_surplus(ex6_v2_sim, "producer", true))
    push!(ex6_v2_consumer_surplus, calculate_surplus(ex6_v2_sim, "consumer,total",true))
    push!(ex6_v2_quality, getfield.(ex6_v2_sim.sellers, :quality_history))
    push!(ex6_v2_durability, getfield.(ex6_v2_sim.sellers, :durability_history))
    push!(ex6_v2_margin, getfield.(ex6_v2_sim.sellers, :margin_history))
    push!(ex6_v2_price, calculate_price_history.(ex6_v2_sim.sellers; product_life = 5))
    push!(ex6_v2_quantity_produced, getfield.(ex6_v2_sim.sellers, :quantity_produced_history))
    push!(ex6_v2_quantity_sold, getfield.(ex6_v2_sim.sellers, :quantity_sold_history))
    push!(ex6_v2_quantity_leased, getfield.(ex6_v2_sim.sellers, :quantity_leased_history))
    push!(ex6_v2_producer_surplus_singleton, calculate_profit_history.(ex6_v2_sim.sellers))
    push!(ex6_v2_reselling, getfield.(ex6_v2_sim.sellers, :reselling_history))
    push!(ex6_v2_buying_history, getfield.(ex6_v2_sim.buyers, :unit_buying_selling_history))

end

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

ex2_p1 = Plots.scatter(mean.(getindex.(ex6_v1_margin, 1))[mean.(getindex.(ex6_v1_margin, 1)) .>= 1.2], mean.(getindex.(ex6_v1_producer_surplus_singleton, 1))[mean.(getindex.(ex6_v1_margin, 1)) .>= 1.2], smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0, xlabel = "Marża %", ylabel = "Zysk firmy", label = "Producent nie bada konsumentów")
#Plots.scatter!(mean.(getindex.(ex6_v1_margin, 2))[mean.(getindex.(ex6_v1_margin, 2)) .>= 1.2], mean.(getindex.(ex6_v1_producer_surplus_singleton, 2))[mean.(getindex.(ex6_v1_margin, 2)) .>= 1.2], smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0)
Plots.scatter!(mean.(getindex.(ex6_v2_margin, 1))[mean.(getindex.(ex6_v2_margin, 1)) .>= 1.2], mean.(getindex.(ex6_v2_producer_surplus_singleton, 1))[mean.(getindex.(ex6_v2_margin, 1)) .>= 1.2], smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0, label = "Producent bada konsumentów")
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
NE = [(strategies[ne[1]], strategies[ne[2]]) for ne in NE]

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

