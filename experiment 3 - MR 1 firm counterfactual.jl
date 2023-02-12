# Experiment 3

include(pwd() * "\\methods\\methods.jl")

println(Threads.nthreads())

# Impact on social communication of simulation dynamics

ex3_v100_total_surplus = []
ex3_v100_producer_surplus = []
ex3_v100_consumer_surplus = []
ex3_v100_price = []
ex3_v100_quantity_produced = []
ex3_v100_quantity_sold = []
ex3_v100_quantity_leased = []
ex3_v100_reselling = []
ex3_v100_producer_surplus_singleton = []
ex3_v100_buying_history = []
ex3_v100_quality = []
ex3_v100_durability = []
ex3_v100_margin = []
ex3_v100_quality_exp = []
ex3_v100_durability_exp = []
ex3_v100_H = []
ex3_v100_Li = []
ex3_v100_Lw = []
ex3_v100_sv = []

ex3_v101_total_surplus = []
ex3_v101_producer_surplus = []
ex3_v101_consumer_surplus = []
ex3_v101_price = []
ex3_v101_quantity_produced = []
ex3_v101_quantity_sold = []
ex3_v101_quantity_leased = []
ex3_v101_reselling = []
ex3_v101_producer_surplus_singleton = []
ex3_v101_buying_history = []
ex3_v101_quality = []
ex3_v101_durability = []
ex3_v101_margin = []
ex3_v101_quality_exp = []
ex3_v101_durability_exp = []
ex3_v101_sv = []

ex3_v100101_H = []
ex3_v100101_Li = []
ex3_v100101_Lw = []
ex3_v100101_sd = []
ex3_v100101_cu = []
ex3_v100101_sm = []
ex3_v100101_nl = []
ex3_v100101_ρm = []
ex3_v100101_cc = []
ex3_v100101_m = []

for i in 1:100

    println(i)

    li = sample(0.0:0.25:1.0)
    lw = sample(0.0:0.25:1.0)
    h = sample(2:8)
    sd = sample(1:100000)
    cu = sample(0.1:0.1:0.5)
    sm = sample(0.025:0.025:0.150)
    nl = sample(200:200:2400)
    ρm = rand(Uniform(0.5, 0.9))
    cc = sample(0.3:0.1:0.6)
    m = sample(1.0:0.1:1.8)

    push!(ex3_v100101_H, h)
    push!(ex3_v100101_Li, li)
    push!(ex3_v100101_Lw, lw)
    push!(ex3_v100101_sd, sd)
    push!(ex3_v100101_cu, cu)
    push!(ex3_v100101_sm, sm)
    push!(ex3_v100101_nl, nl)
    push!(ex3_v100101_ρm, ρm)
    push!(ex3_v100101_cc, cc)
    push!(ex3_v100101_m, m)

    Random.seed!(sd)

    ex3_v100_sim = TO_GO(200, 2, 400, nl, [cc, cc], [m, m], "random", li, lw, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], cu, true, 1, [ρm, 1.], "softmax", [false, false], [0., 0.], h, false, true)

    push!(ex3_v100_total_surplus, calculate_surplus(ex3_v100_sim, "total", false))
    push!(ex3_v100_producer_surplus, calculate_surplus(ex3_v100_sim, "producer", false))
    push!(ex3_v100_consumer_surplus, calculate_surplus(ex3_v100_sim, "consumer,total",false))
    push!(ex3_v100_quality, getfield.(ex3_v100_sim.sellers, :quality_history))
    push!(ex3_v100_durability, getfield.(ex3_v100_sim.sellers, :durability_history))
    push!(ex3_v100_margin, getfield.(ex3_v100_sim.sellers, :margin_history))
    push!(ex3_v100_price, calculate_price_history.(ex3_v100_sim.sellers; product_life = h))
    push!(ex3_v100_quantity_produced, getfield.(ex3_v100_sim.sellers, :quantity_produced_history))
    push!(ex3_v100_quantity_sold, getfield.(ex3_v100_sim.sellers, :quantity_sold_history))
    push!(ex3_v100_producer_surplus_singleton, calculate_profit_history.(ex3_v100_sim.sellers))
    push!(ex3_v100_reselling, getfield.(ex3_v100_sim.sellers, :reselling_history))
    push!(ex3_v100_buying_history, getfield.(ex3_v100_sim.buyers, :unit_buying_selling_history))
    push!(ex3_v100_quality_exp, mean([b.quality_expectation_history for b in ex3_v100_sim.buyers]))
    push!(ex3_v100_durability_exp, mean([b.durability_expectation_history for b in ex3_v100_sim.buyers]))
    push!(ex3_v100_sv, getfield.(ex3_v100_sim.buyers, :signal_volume))

    Random.seed!(sd)

    ex3_v101_sim = TO_GO(200, 2, 400, nl, [cc, cc], [m, m], "random", li, lw, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], cu, true, 1, [ρm, 1.], "softmax", [true, false], [sm, 0.], h, false, true)

    push!(ex3_v101_total_surplus, calculate_surplus(ex3_v101_sim, "total", false))
    push!(ex3_v101_producer_surplus, calculate_surplus(ex3_v101_sim, "producer", false))
    push!(ex3_v101_consumer_surplus, calculate_surplus(ex3_v101_sim, "consumer,total",false))
    push!(ex3_v101_quality, getfield.(ex3_v101_sim.sellers, :quality_history))
    push!(ex3_v101_durability, getfield.(ex3_v101_sim.sellers, :durability_history))
    push!(ex3_v101_margin, getfield.(ex3_v101_sim.sellers, :margin_history))
    push!(ex3_v101_price, calculate_price_history.(ex3_v101_sim.sellers; product_life = h))
    push!(ex3_v101_quantity_produced, getfield.(ex3_v101_sim.sellers, :quantity_produced_history))
    push!(ex3_v101_quantity_sold, getfield.(ex3_v101_sim.sellers, :quantity_sold_history))
    push!(ex3_v101_producer_surplus_singleton, calculate_profit_history.(ex3_v101_sim.sellers))
    push!(ex3_v101_reselling, getfield.(ex3_v101_sim.sellers, :reselling_history))
    push!(ex3_v101_buying_history, getfield.(ex3_v101_sim.buyers, :unit_buying_selling_history))
    push!(ex3_v101_quality_exp, mean([b.quality_expectation_history for b in ex3_v101_sim.buyers]))
    push!(ex3_v101_durability_exp, mean([b.durability_expectation_history for b in ex3_v101_sim.buyers]))
    push!(ex3_v101_sv, getfield.(ex3_v101_sim.buyers, :signal_volume))

end

RMSE(x,y) = sqrt(mean((x .- y).^2))
xydiff(x,y) = mean(x .- y)
cv(x) = std(x) / mean(x)

###############################

Plots.plot(mean(getindex.(ex3_v100_quality,1)))
Plots.plot!(mean(getindex.(ex3_v101_quality,1)))
Plots.plot!(mean(getindex.(ex3_v100_quality,2)))
Plots.plot!(mean(getindex.(ex3_v101_quality,2)))

Plots.plot(mean(getindex.(ex3_v100_durability,1)))
Plots.plot!(mean(getindex.(ex3_v101_durability,1)))
Plots.plot!(mean(getindex.(ex3_v100_durability,2)))
Plots.plot!(mean(getindex.(ex3_v101_durability,2)))

Plots.plot(mean(getindex.(ex3_v100_margin,1)))
Plots.plot!(mean(getindex.(ex3_v101_margin,1)))
Plots.plot!(mean(getindex.(ex3_v100_margin,2)))
Plots.plot!(mean(getindex.(ex3_v101_margin,2)))

Plots.plot(mean(getindex.(ex3_v100_producer_surplus_singleton, 1)))
Plots.plot!(mean(getindex.(ex3_v101_producer_surplus_singleton, 1)))
Plots.plot!(mean(getindex.(ex3_v100_producer_surplus_singleton, 2)))
Plots.plot!(mean(getindex.(ex3_v101_producer_surplus_singleton, 2)))

Plots.plot(mean(getindex.(ex3_v100_quality,1)))
Plots.plot!(mean([getindex.(x,1) for x in ex3_v100_quality_exp]))

Plots.plot(mean(getindex.(ex3_v100_quality,2)))
Plots.plot!(mean([getindex.(x,2) for x in ex3_v100_quality_exp]))

Plots.plot(mean(getindex.(ex3_v101_quality,2)))
Plots.plot!(mean([getindex.(x,1) for x in ex3_v101_quality_exp]))

###############################

# Gęstość, różnice między oczekiwaniami a realnymi charakterystykami

ex3_p1_pl = StatsPlots.density(xydiff.(getindex.(ex3_v100_quality, 1), [getindex.(x,1) for x in ex3_v100_quality_exp]), label = "Producent nie prowadzi badań", xlabel = "Odchylenie parametrów produktu [jakość] od oczekiwań konsumentów", ylabel = "f(x)", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, title = "Funkcja gęstości odchyleń parametrów produktu od oczekiwań konsumentów", legendfontsize = 6, legend = :topleft, normalize = true)
StatsPlots.density!(xydiff.(getindex.(ex3_v101_quality, 1), [getindex.(x,1) for x in ex3_v101_quality_exp]), label = "Producent prowadzi badania")
Plots.plot!([0,0], [0,10], color = "black")

savefigs(ex3_p1_pl, "\\plots\\ex2\\PL density diff expectations vs average")

mean(xydiff.(getindex.(ex3_v100_quality, 1), [getindex.(x,1) for x in ex3_v100_quality_exp]))

StatsPlots.density(mean.(getindex.(ex3_v100_quality, 1)))


ex3_p1_eng = StatsPlots.density(xydiff.(getindex.(ex3_v100_quality, 1), [getindex.(x,1) for x in ex3_v100_quality_exp]), label = "No research", xlabel = "Average Quality deviation from consumer expectation", ylabel = "f(x)", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, title = "Density function", legendfontsize = 6, legend = :topleft, normalize = true)
StatsPlots.density!(xydiff.(getindex.(ex3_v101_quality, 1), [getindex.(x,1) for x in ex3_v101_quality_exp]), label = "Research")
Plots.plot!([0,0], [0,10], color = "black")

savefigs(ex3_p1_eng, "\\plots\\ex2\\ENG density diff expectations vs average")

####################################################################

p_ex3_2_pl = Plots.boxplot(xydiff.(getindex.(ex3_v100_quality, 1), [getindex.(x,1) for x in ex3_v100_quality_exp]), ylabel = "Odchylenie parametrów produktu od oczekiwań", legend = nothing, ylabelfontsize = 8, dpi = 300, title = "Jakość")
Plots.boxplot!(xydiff.(getindex.(ex3_v101_quality, 1), [getindex.(x,1) for x in ex3_v101_quality_exp]))
Plots.plot!(xticks = ([1,2], ["Producent nie prowadzi badań", "Producent prowadzi badania"]))

savefigs(p_ex3_2_pl, "\\plots\\ex2\\PL boxplot diff expectations vs average quality")

p_ex3_2_eng = Plots.boxplot(xydiff.(getindex.(ex3_v100_quality, 1), [getindex.(x,1) for x in ex3_v100_quality_exp]), ylabel = "Deviation of average parameters from expectations", legend = nothing, ylabelfontsize = 8, dpi = 300, title = "Quality")
Plots.boxplot!(xydiff.(getindex.(ex3_v101_quality, 1), [getindex.(x,1) for x in ex3_v101_quality_exp]))
Plots.plot!(xticks = ([1,2], ["No research", "Research"]))

savefigs(p_ex3_2_eng, "\\plots\\ex2\\ENG boxplot diff expectations vs average quality")

####################################################################

p_ex3_3_pl = Plots.boxplot(xydiff.(getindex.(ex3_v100_durability, 1), [getindex.(x,1) for x in ex3_v100_durability_exp]), ylabel = "Odchylenie parametrów produktu od oczekiwań", legend = nothing, ylabelfontsize = 8, dpi = 300, title = "Trwałość")
Plots.boxplot!(xydiff.(getindex.(ex3_v101_durability, 1), [getindex.(x,1) for x in ex3_v101_durability_exp]))
Plots.plot!(xticks = ([1,2], ["Producent nie prowadzi badań", "Producent prowadzi badania"]))

savefigs(p_ex3_3_pl, "\\plots\\ex2\\PL boxplot diff expectations vs average durability")

p_ex3_3_eng = Plots.boxplot(xydiff.(getindex.(ex3_v100_durability, 1), [getindex.(x,1) for x in ex3_v100_durability_exp]), ylabel = "Deviation of average parameters from expectations", legend = nothing, ylabelfontsize = 8, dpi = 300, title = "Durability")
Plots.boxplot!(xydiff.(getindex.(ex3_v101_durability, 1), [getindex.(x,1) for x in ex3_v101_durability_exp]))
Plots.plot!(xticks = ([1,2], ["No research", "Research"]))

savefigs(p_ex3_3_eng, "\\plots\\ex2\\ENG boxplot diff expectations vs average durability")

#############################################################################################################

ex3_4_pl = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)), "Producent bada oczekiwania konsumentów")

savefigs(ex3_4_pl, "\\plots\\ex1\\PL ECDF producer profit")

ex3_4_eng = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)), "No research", xlabel = "Profit", ylabel = "F(x)", title = "ECDF - Profit")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)), "Research")

savefigs(ex3_4_eng, "\\plots\\ex1\\ENG ECDF producer profit")

###########################################################################################

gbp_df = DataFrame(x = repeat(ex3_v100101_H, 2), y = vcat(sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)), sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))), g = repeat( ["Firma nie prowadzi badań konsumenckich", "Firma prowadzi badania konsumenckie"],inner = lastindex(ex3_v100101_H)))

gbp_df = sort!(gbp_df, :x)

ex3_5_pl = @df gbp_df groupedboxplot(:x, :y, group = :g, xlabel = "Liczba okresów przydatności dobra", ylabel = "Zysk firmy")

percentile(sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)), 5)
percentile(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)), 5)

###########################################################################################

ex3_6_pl = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_H .< 5], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "",xlim = (-2000, 5000))
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_H .< 5], "Producent bada oczekiwania konsumentów")

ex3_7_pl = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_H .>= 5], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "", xlim = (-2000,5000))
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_H .>= 5], "Producent bada oczekiwania konsumentów")

percentile(sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_H .< 5], 5)
percentile(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_H .< 5], 5)

percentile(sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_H .>= 5], 5)
percentile(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_H .>= 5], 5)

Plots.plot(ex3_6_pl, ex3_7_pl, layout = (1,2), plot_title = "Dystybuanta empiryczna - zysk", plot_titlefontsize = 12)

############################################################################################


ex3_8_pl = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_nl .<= 300], "", xlabel = "Zysk producenta", ylabel = "F(x)", title = "<300")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_nl .<= 300], "")

ex3_9_pl = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 400) .& (ex3_v100101_nl .<= 600)], "", xlabel = "Zysk producenta", ylabel = "F(x)", title = "400-600")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 400) .& (ex3_v100101_nl .<= 600)], "")

ex3_10_pl = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 700) .& (ex3_v100101_nl .<= 900)], "", xlabel = "Zysk producenta", ylabel = "F(x)", title = "700-900")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 700) .& (ex3_v100101_nl .<= 900)], "")

ex3_11_pl = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_nl .>= 1000], "", xlabel = "Zysk producenta", ylabel = "F(x)", title = "1000-1200")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_nl .>= 1000], "")

Plots.plot(ex3_8_pl, ex3_9_pl, ex3_10_pl, ex3_11_pl)

############################################################################################

bounds = collect(-2000:1000:5000)

function cut_integer_bounds(x::Vector{Float64},bounds::Vector)

    value_counts = Int64[]

    for item in eachindex(bounds)[Not(length(bounds))]
        push!(value_counts, sum((x .>= bounds[item]) .& (x .< bounds[item+1])))
    end

    return value_counts

end

function earthmoverdistance(a::Vector, b::Vector) 
    return sum( abs, cumsum( a ) .- cumsum( b ) )
end

emd300 = earthmoverdistance(
    cut_integer_bounds(sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_nl .<= 300], bounds),
    cut_integer_bounds(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_nl .<= 300], bounds)
) / sum(ex3_v100101_nl .<= 300)

emd600 = earthmoverdistance(
    cut_integer_bounds(sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 400) .& (ex3_v100101_nl .<= 600)], bounds),
    cut_integer_bounds(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 400) .& (ex3_v100101_nl .<= 600)], bounds)
) / sum((ex3_v100101_nl .>= 400) .& (ex3_v100101_nl .<= 600))

emd900 = earthmoverdistance(
    cut_integer_bounds(sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 700) .& (ex3_v100101_nl .<= 900)], bounds),
    cut_integer_bounds(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 700) .& (ex3_v100101_nl .<= 900)], bounds)
) / sum((ex3_v100101_nl .>= 700) .& (ex3_v100101_nl .<= 900))

emd1200 = earthmoverdistance(
    cut_integer_bounds(sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_nl .>= 1000], bounds),
    cut_integer_bounds(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_nl .>= 1000], bounds)
) / sum(ex3_v100101_nl .>= 1000)

earthmoverdistance(cut_integer_bounds(sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)), bounds), 
cut_integer_bounds(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)), bounds))

median_distance(x,y) = median(x) - median(y)
mean_distnace(x,y) = mean(x) - mean(y)

emd_df = DataFrame(
    Stan = ["300", "300,700", "700, 1000", "1000"],
    EMD = round.([emd300, emd600, emd900, emd1200], digits = 2),
    Mediana = round.([median_distance(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_nl .<= 300], sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_nl .<= 300]),
    median_distance(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 400) .& (ex3_v100101_nl .<= 600)], sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 400) .& (ex3_v100101_nl .<= 600)]),
    median_distance(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 700) .& (ex3_v100101_nl .<= 900)], sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 700) .& (ex3_v100101_nl .<= 900)]),
    median_distance(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_nl .>= 1000], sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_nl .>= 1000])], digits = 2),
    Średnia = round.([mean_distnace(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_nl .<= 300], sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_nl .<= 300]),
    mean_distnace(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 400) .& (ex3_v100101_nl .<= 600)], sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 400) .& (ex3_v100101_nl .<= 600)]),
    mean_distnace(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 700) .& (ex3_v100101_nl .<= 900)], sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 700) .& (ex3_v100101_nl .<= 900)]),
    mean_distnace(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_nl .>= 1000], sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_nl .>= 1000])], digits = 2)
)

latexify(emd_df, env = :table) |> print

############################################################################################

Plots.plot

plot_ecdf(true, getindex.([var.(x) for x in ex3_v100_producer_surplus_singleton], 1), "", xlabel = "Var (" * L"\pi" * ")", ylabel = "F(x)", title = "Wariancja zysku " * L"\pi")
plot_ecdf(false, getindex.([var.(x) for x in ex3_v101_producer_surplus_singleton], 1), "")

###########################################################################################

Plots.boxplot([sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_sm .== sm] for sm in sort(unique(ex3_v100101_sm))])
Plots.boxplot([sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_sm .== sm] for sm in sort(unique(ex3_v100101_sm))])

Plots.scatter(ex3_v100101_ρm, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)))
Plots.scatter(ex3_v100101_ρm, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)))

Plots.scatter(ex3_v100101_ρm, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)) .- sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)))

Plots.scatter(ex3_v100101_cu, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)), smooth = true)
Plots.scatter(ex3_v100101_cu, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)), smooth = true)

Plots.scatter(ex3_v100101_cu, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)) .- sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)))

Plots.scatter(ex3_v100101_Li, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)), smooth = true)
Plots.scatter(ex3_v100101_Li, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)), smooth = true)

Plots.scatter(ex3_v100101_Li, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)) .- sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)))

Plots.scatter(ex3_v100101_Lw, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)), smooth = true)
Plots.scatter(ex3_v100101_Lw, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)), smooth = true)

##########################################################################################


round(iqr(xydiff.(getindex.(ex3_v100_quality, 1), [getindex.(x,1) for x in ex3_v100_quality_exp])), digits = 3)
round(cv(xydiff.(getindex.(ex3_v100_quality, 1), [getindex.(x,1) for x in ex3_v100_quality_exp])), digits = 3)
round(iqr(xydiff.(getindex.(ex3_v101_quality, 2), [getindex.(x,2) for x in ex3_v101_quality_exp])), digits = 3)
round(cv(xydiff.(getindex.(ex3_v101_quality, 2), [getindex.(x,2) for x in ex3_v101_quality_exp])), digits = 3)

p_ex3_10 = Plots.boxplot(xydiff.(getindex.(ex3_v100_quality, 1), [getindex.(x,1) for x in ex3_v100_quality_exp]), ylabel = "Odchylenie parametrów produktu od oczekiwań", legend = nothing, ylabelfontsize = 8, dpi = 300, title = "Jakość")
Plots.boxplot!(xydiff.(getindex.(ex3_v101_quality, 1), [getindex.(x,1) for x in ex3_v101_quality_exp]))
Plots.plot!(xticks = ([1,2], ["Producent nie prowadzi badań", "Producent prowadzi badania"]))

#%%%

Plots.plot(mean(getindex.(ex3_v100_quality, 1)))
Plots.plot!(mean(getindex.(ex3_v101_quality, 1)))
Plots.plot!(mean([getindex.(x,1) for x in ex3_v100_quality_exp]))
Plots.plot!(mean([getindex.(x,1) for x in ex3_v101_quality_exp]))

Plots.plot(ex3_v100_quality[1][1])

getindex.(ex3_v100_quality_exp,1)


# Trend zmian jakości i trwałości

Plots.plot(mean(getindex.(ex3_v100_quality, 1)))
Plots.plot!(mean(getindex.(ex3_v101_quality, 1)))
Plots.plot!(mean(getindex.(ex3_v100_quality, 2)))
Plots.plot!(mean(getindex.(ex3_v101_quality, 2)))

Plots.plot(mean(getindex.(ex3_v100_durability, 1)))
Plots.plot!(mean(getindex.(ex3_v101_durability, 1)))
Plots.plot!(mean(getindex.(ex3_v100_durability, 2)))
Plots.plot!(mean(getindex.(ex3_v101_durability, 2)))

#%%%%

Plots.plot((ex3_v100_quality)[1][1])
Plots.plot!((ex3_v101_quality)[1][1])

StatsPlots.density(mean.(ex3_v100_sv))
StatsPlots.density!(mean.(ex3_v101_sv))

ApproximateTwoSampleKSTest(mean.(ex3_v100_sv), mean.(ex3_v101_sv))

Plots.scatter(mean.(ex3_v100_sv), sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)), smooth = true)
Plots.scatter!(mean.(ex3_v101_sv), sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)), smooth = true)

Plots.scatter(mean.(ex3_v100_sv)[ex3_v100101_H .< 3], sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_H .< 3], smooth = true)
Plots.scatter!(mean.(ex3_v101_sv), sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)), smooth = true)

ex1_p10 = StatsPlots.boxplot([sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_H .== x] for x in sort(unique(ex3_v100101_H))], legend = false, xlabel = "Maksymalna przydatność dobra [liczba okresów]", ylabel = "Zysk firmy", title = "Zysk firmy, która nie bada konsumentów", ylim = (-1500, 3000))

ex1_p10 = StatsPlots.boxplot([sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_H .== x] for x in sort(unique(ex3_v100101_H))], legend = false, xlabel = "Maksymalna przydatność dobra [liczba okresów]", ylabel = "Zysk firmy", title = "Zysk firmy, która nie bada konsumentów", ylim = (-1500, 3000))

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v100_producer_surplus_singleton, 2)), "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 2)), "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)), "Producent bada oczekiwania konsumentów")

UnequalVarianceTTest(sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)), sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)))

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_H .< 5], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_H .< 5], "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_H .>= 5], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_H .>= 5], "Producent bada oczekiwania konsumentów")

#########################################################################################

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_nl .<= 300], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_nl .<= 300], "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 400) .& (ex3_v100101_nl .<= 600)], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 400) .& (ex3_v100101_nl .<= 600)], "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 700) .& (ex3_v100101_nl .<= 900)], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[(ex3_v100101_nl .>= 700) .& (ex3_v100101_nl .<= 900)], "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))[ex3_v100101_nl .>= 1000], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))[ex3_v100101_nl .>= 1000], "Producent bada oczekiwania konsumentów")

####

mean(sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)))
mean(sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)))

mean(sum.(getindex.(ex3_v101_quantity_sold, 1)))
mean(sum.(getindex.(ex3_v100_quantity_sold, 1)))

mean(sum.(getindex.(ex3_v101_quantity_produced, 1)))
mean(sum.(getindex.(ex3_v100_quantity_produced, 1)))

mean(mean.(getindex.(ex3_v101_price, 1)))
mean(mean.(getindex.(ex3_v100_price, 1)))

mean(mean.(getindex.(ex3_v101_price, 1))) / mean(mean.(getindex.(ex3_v101_margin, 1)))
mean(mean.(getindex.(ex3_v100_price, 1))) / mean(mean.(getindex.(ex3_v100_margin, 1)))

[sum(ex3_v101_quantity_sold[k][1] .* ex3_v101_price[k][1]) - sum(ex3_v101_quantity_produced[k][1] .* ex3_v101_price[k][1] ./ ex3_v101_margin[k][1]) .+ sum((ex3_v101_quantity_produced[k][1] .- ex3_v101_quantity_sold[k][1]) * 0.5 .* ex3_v101_price[k][1] ./ ex3_v101_margin[k][1]) for k in 1:400] .- sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))

mean()

k=3
sum(ex3_v101_quantity_sold[k][1] .* ex3_v101_price[k][1]) - sum(ex3_v101_quantity_produced[k][1] .* ex3_v101_price[k][1] ./ ex3_v101_margin[k][1]) .+ sum((ex3_v101_quantity_produced[k][1] .- ex3_v101_quantity_sold[k][1]) * 0.5 .* ex3_v101_price[k][1] ./ ex3_v101_margin[k][1])

sum(ex3_v100_quantity_sold[k][1] .* ex3_v100_price[k][1]) - sum(ex3_v100_quantity_produced[k][1] .* ex3_v100_price[k][1] ./ ex3_v100_margin[k][1]) .+ sum((ex3_v100_quantity_produced[k][1] .- ex3_v100_quantity_sold[k][1]) * 0.5 .* ex3_v100_price[k][1] ./ ex3_v100_margin[k][1])

calculate_utility(k,d,h) = k .* (1 .- d .^ h) ./ (1 .- d)

Plots.plot(calculate_utility(ex3_v101_quality[k][1], ex3_v101_durability[k][1], ex3_v100101_H[k]))
Plots.plot!(calculate_utility(ex3_v100_quality[k][1], ex3_v100_durability[k][1], ex3_v100101_H[k]))

Plots.plot!(calculate_utility(getindex.(ex3_v101_quality_exp[k], 1), getindex.(ex3_v101_durability_exp[k], 1), ex3_v100101_H[k]))
Plots.plot!(calculate_utility(getindex.(ex3_v100_quality_exp[k], 1), getindex.(ex3_v100_durability_exp[k], 1), ex3_v100101_H[k]))

Plots.plot(ex3_v101_price[k][1])
Plots.plot!(ex3_v100_price[k][1])

Plots.plot(calculate_utility(ex3_v101_quality[k][1], ex3_v101_durability[k][1], ex3_v100101_H[k]) ./ ex3_v101_price[k][1])
Plots.plot!(calculate_utility(ex3_v100_quality[k][1], ex3_v100_durability[k][1], ex3_v100101_H[k]) ./ ex3_v100_price[k][1])

Plots.plot!(calculate_utility(getindex.(ex3_v101_quality_exp[k], 1), getindex.(ex3_v101_durability_exp[k], 1), ex3_v100101_H[k]) ./ ex3_v101_price[k][1])
Plots.plot!(calculate_utility(getindex.(ex3_v100_quality_exp[k], 1), getindex.(ex3_v100_durability_exp[k], 1), ex3_v100101_H[k]) ./ ex3_v100_price[k][1])

CSV.write(pwd() * "\\data_test.csv", DataFrame(b_qs = ex3_v101_quantity_sold[k][1], b_p = ex3_v101_price[k][1], b_qp = ex3_v101_quantity_produced[k][1], b_m = ex3_v101_margin[k][1], nb_qs = ex3_v100_quantity_sold[k][1], nb_p = ex3_v100_price[k][1], nb_qp = ex3_v100_quantity_produced[k][1], nb_m = ex3_v100_margin[k][1]))

#####

function consecutive_negative(x)
    current_counter = 0
    best_counter = 0
    for xi in x
        if xi < 0
            current_counter += 1
        else
            if current_counter > best_counter
                best_counter = current_counter
            end

            current_counter = 0

        end
    end
    return best_counter
end

StatsPlots.density(consecutive_negative.(getindex.(ex3_v100_producer_surplus_singleton, 1)))
StatsPlots.density!(consecutive_negative.(getindex.(ex3_v101_producer_surplus_singleton, 1)))

####

ex3_p1 = plot_ecdf(true, mean.(getindex.(ex3_v100_quality, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Jakość producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, mean.(getindex.(ex3_v101_quality, 1)), "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, mean.(getindex.(ex3_v100_durability, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Jakość producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, mean.(getindex.(ex3_v101_durability, 1)), "Producent bada oczekiwania konsumentów")



plot_ecdf(true, sum.(ex3_v100_total_surplus), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(ex3_v101_total_surplus), "Producent bada oczekiwania konsumentów")

plot_ecdf(true, sum.(ex3_v100_consumer_surplus), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(ex3_v101_consumer_surplus), "Producent bada oczekiwania konsumentów")

#%%%

function speed_of_learing(x, xe)
    dx = diff(x)
    dxe = diff(xe)
    return sum(sign.(dx) .== sign.(dxe))
end

RMSE.([getindex.(x,1) for x in ex3_v100_quality_exp], getindex.(ex3_v100_quality, 1))

Plots.scatter(RMSE.([getindex.(x,1) for x in ex3_v100_quality_exp], getindex.(ex3_v100_quality, 1)), sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)))
Plots.scatter!(RMSE.([getindex.(x,1) for x in ex3_v101_quality_exp], getindex.(ex3_v101_quality, 1)), sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)))

Plots.scatter(speed_of_learing.([getindex.(x,1) for x in ex3_v100_quality_exp], getindex.(ex3_v100_quality, 1)), sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)))
Plots.scatter!(speed_of_learing.([getindex.(x,1) for x in ex3_v101_quality_exp], getindex.(ex3_v101_quality, 1)), sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)))

Plots.scatter(std.([getindex.(x,1) for x in ex3_v100_quality_exp]), sum.(getindex.(ex3_v100_producer_surplus_singleton, 1)))

Plots.scatter!(std.([getindex.(x,1) for x in ex3_v101_quality_exp]), sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)))

sv100 = cut_integer(mean.(ex3_v100_sv), 10)[1]
sv101 = cut_integer(mean.(ex3_v101_sv), 10)[1]

[std(mean.(getindex.(ex3_v100_producer_surplus_singleton, 1))[sv100 .== x]) for x in sort(unique(sv100))]

[std(mean.(getindex.(ex3_v101_producer_surplus_singleton, 1))[sv101 .== x]) for x in sort(unique(sv101))]

Plots.scatter(mean.(ex3_v100_sv), mean.(getindex.(ex3_v100_producer_surplus_singleton, 1)))
Plots.scatter!(mean.(ex3_v101_sv), mean.(getindex.(ex3_v101_producer_surplus_singleton, 1)))

mean.(ex3_v100_sv) .- mean.(ex3_v101_sv)

Plots.scatter(mean.(ex3_v100_sv), mean.(ex3_v101_sv), smooth = true, xlim = (0,150), ylim = (0,150))


Plots.scatter(mean.(ex3_v100101_Lw), mean.(getindex.(ex3_v100_producer_surplus_singleton, 1)) .- mean.(getindex.(ex3_v101_producer_surplus_singleton, 1)))

λ = cut_integer(float.(ex3_v100101_Li), 5)[1]
θ = cut_integer(float.(ex3_v100101_Lw), 5)[1]

hm_rs = [mean((mean.(getindex.(ex3_v101_producer_surplus_singleton, 1)) .- mean.(getindex.(ex3_v100_producer_surplus_singleton, 1)))[(λ .== l) .& (θ .== t)]) for l in sort(unique(λ)), t in sort(unique(θ))]

mean((mean.(getindex.(ex3_v101_producer_surplus_singleton, 1)) .- mean.(getindex.(ex3_v100_producer_surplus_singleton, 1)))[(λ .== 2) .& (θ .== 3)])

mean((mean.(getindex.(ex3_v101_producer_surplus_singleton, 1)) .- mean.(getindex.(ex3_v100_producer_surplus_singleton, 1)))[(λ .== 3) .& (θ .== 2)])

Plots.plot([mean(x) for x in eachrow(hm_rs)]) # by λ
Plots.plot([mean(x) for x in eachcol(hm_rs)]) # by θ

ex4_p1 = StatsPlots.heatmap(sort(unique(λ)), sort(unique(θ)), hm_rs', xlabel = "λ, produkty oceniane osobiście", ylabel = "θ, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita", titlefontsize = 8)

####

Y = sum.(getindex.(ex3_v101_producer_surplus_singleton, 1)) .- sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))

Y1 = sum.(getindex.(ex3_v100_producer_surplus_singleton, 1))
Y2 = sum.(getindex.(ex3_v101_producer_surplus_singleton, 1))

H = float.(ex3_v100101_H)
λ = float.(ex3_v100101_Li)
θ = float.(ex3_v100101_Lw)
sd = float.(ex3_v100101_sd)
cu = float.(ex3_v100101_cu)
sm = float.(ex3_v100101_sm)
nl = float.(ex3_v100101_nl)
rm = float.(ex3_v100101_ρm)


k1 = float.(mean.(getindex.(ex3_v100_quality, 1)))
k2 = float.(mean.(getindex.(ex3_v101_quality, 1)))

Y = k2.-k1

d1 = float.(mean.(getindex.(ex3_v100_durability, 1)))
d2 = float.(mean.(getindex.(ex3_v101_durability, 1)))

m1 = float.(mean.(getindex.(ex3_v100_margin, 1)))
m2 = float.(mean.(getindex.(ex3_v101_margin, 1)))

r1 = mean.(getindex.(ex3_v100_reselling, 1))
r2 = mean.(getindex.(ex3_v101_reselling, 1))

ek1 = float.(mean.([getindex.(x,1) for x in ex3_v100_quality_exp]))
ek2 = float.(mean.([getindex.(x,1) for x in ex3_v101_quality_exp]))

ed1 = float.(mean.([getindex.(x,1) for x in ex3_v100_durability_exp]))
ed2 = float.(mean.([getindex.(x,1) for x in ex3_v101_durability_exp]))

optimal_research_100 = GLM.lm(@formula(Y1 ~ λ + θ + k1 + d1 + r1 + ek1 + ed1), DataFrame(Y1 = Y1, H = H, λ = λ, θ = θ, sd = sd, k1 = k1, k2 = k2, d1 = d1, d2 = d2, m1 = m1, m2 = m2, r1 = r1, r2 = r2, ek1 = ek1, ek2 = ek2, ed1 = ed1, ed2 = ed2))

round(adjr2(optimal_research_100), digits = 2)

optimal_research_101 = GLM.lm(@formula(Y2 ~ λ + θ + k2 + d2 + r2 + ek2 + ed2), DataFrame(Y1 = Y1, Y2 = Y2, H = H, λ = λ, θ = θ, sd = sd, k1 = k1, k2 = k2, d1 = d1, d2 = d2, m1 = m1, m2 = m2, r1 = r1, r2 = r2, ek1 = ek1, ek2 = ek2, ed1 = ed1, ed2 = ed2))

round(adjr2(optimal_research_101), digits = 2)



####