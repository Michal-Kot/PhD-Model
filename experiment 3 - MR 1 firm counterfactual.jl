# Experiment 3

include(pwd() * "\\methods\\methods.jl")

# Impact on social communication of simulation dynamics

ex3_v1_total_surplus = []
ex3_v1_producer_surplus = []
ex3_v1_consumer_surplus = []
ex3_v1_price = []
ex3_v1_quantity_produced = []
ex3_v1_quantity_sold = []
ex3_v1_quantity_leased = []
ex3_v1_reselling = []
ex3_v1_producer_surplus_singleton = []
ex3_v1_buying_history = []
ex3_v1_quality = []
ex3_v1_durability = []
ex3_v1_margin = []
ex3_v1_quality_exp = []
ex3_v1_durability_exp = []
ex3_v1_H = []
ex3_v1_Li = []
ex3_v1_Lw = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    li = rand()
    lw = rand()
    h = sample(2:8)

    push!(ex3_v1_H, h)
    push!(ex3_v1_Li, li)
    push!(ex3_v1_Lw, lw)

    ex3_v1_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1.1, 1.1], "random", li, lw, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], 0.50, true, 1, [0.7, 1.], "softmax", [false, true], [0, 0.1], h, false)

    push!(ex3_v1_total_surplus, calculate_surplus(ex3_v1_sim, "total", false))
    push!(ex3_v1_producer_surplus, calculate_surplus(ex3_v1_sim, "producer", false))
    push!(ex3_v1_consumer_surplus, calculate_surplus(ex3_v1_sim, "consumer,total",false))
    push!(ex3_v1_quality, getfield.(ex3_v1_sim.sellers, :quality_history))
    push!(ex3_v1_durability, getfield.(ex3_v1_sim.sellers, :durability_history))
    push!(ex3_v1_margin, getfield.(ex3_v1_sim.sellers, :margin_history))
    push!(ex3_v1_price, calculate_price_history.(ex3_v1_sim.sellers; product_life = h))
    push!(ex3_v1_quantity_produced, getfield.(ex3_v1_sim.sellers, :quantity_produced_history))
    push!(ex3_v1_quantity_sold, getfield.(ex3_v1_sim.sellers, :quantity_sold_history))
    push!(ex3_v1_producer_surplus_singleton, calculate_profit_history.(ex3_v1_sim.sellers))
    push!(ex3_v1_reselling, getfield.(ex3_v1_sim.sellers, :reselling_history))
    push!(ex3_v1_buying_history, getfield.(ex3_v1_sim.buyers, :unit_buying_selling_history))
    push!(ex3_v1_quality_exp, mean([b.quality_expectation_history for b in ex3_v1_sim.buyers]))
    push!(ex3_v1_durability_exp, mean([b.durability_expectation_history for b in ex3_v1_sim.buyers]))
end

RMSE(x,y) = sum((x .- y).^2)
xydiff(x,y) = mean(x .- y)
cv(x) = std(x) / mean(x)

ex3_p1 = StatsPlots.density(xydiff.(getindex.(ex3_v1_quality, 1), [getindex.(x,1) for x in ex3_v1_quality_exp]), label = "Producent nie prowadzi badań", xlabel = "Odchylenie parametrów produktu [jakość] od oczekiwań konsumentów", ylabel = "f(x)", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, title = "Funkcja gęstości odchyleń parametrów produktu od oczekiwań konsumentów", legendfontsize = 6, legend = :topleft, normalize = true)
StatsPlots.density!(xydiff.(getindex.(ex3_v1_quality, 2), [getindex.(x,2) for x in ex3_v1_quality_exp]), label = "Producent prowadzi badania")

Plots.savefig(ex3_p1, pwd() * "\\plots\\ex2\\quality_research_no_research.svg")

iqr(xydiff.(getindex.(ex3_v1_quality, 1), [getindex.(x,1) for x in ex3_v1_quality_exp]))
cv(xydiff.(getindex.(ex3_v1_quality, 1), [getindex.(x,1) for x in ex3_v1_quality_exp]))
iqr(xydiff.(getindex.(ex3_v1_quality, 2), [getindex.(x,2) for x in ex3_v1_quality_exp]))
cv(xydiff.(getindex.(ex3_v1_quality, 2), [getindex.(x,2) for x in ex3_v1_quality_exp]))

ex3_p2 = StatsPlots.density(xydiff.(getindex.(ex3_v1_durability, 1), [getindex.(x,1) for x in ex3_v1_durability_exp]), label = "Producent nie prowadzi badań", xlabel = "Odchylenie parametrów produktu [trwałość] od oczekiwań konsumentów", ylabel = "f(x)", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, title = "Funkcja gęstości odchyleń parametrów produktu od oczekiwań konsumentów", legendfontsize = 6, legend = :topleft, normalize = true)
StatsPlots.density!(xydiff.(getindex.(ex3_v1_durability, 2), [getindex.(x,2) for x in ex3_v1_durability_exp]), label = "Producent prowadzi badania")

iqr(xydiff.(getindex.(ex3_v1_durability, 1), [getindex.(x,1) for x in ex3_v1_durability_exp]))
cv(xydiff.(getindex.(ex3_v1_durability, 1), [getindex.(x,1) for x in ex3_v1_durability_exp]))
iqr(xydiff.(getindex.(ex3_v1_durability, 2), [getindex.(x,2) for x in ex3_v1_durability_exp]))
cv(xydiff.(getindex.(ex3_v1_durability, 2), [getindex.(x,2) for x in ex3_v1_durability_exp]))

Plots.savefig(ex3_p2, pwd() * "\\plots\\ex2\\durability_research_no_research.svg")



ex1_p10 = StatsPlots.boxplot([sum.(getindex.(ex3_v1_producer_surplus_singleton, 1))[ex3_v1_H .== x] for x in sort(unique(ex3_v1_H))], legend = false, xlabel = "Maksymalna przydatność dobra [liczba okresów]", ylabel = "Zysk firmy", title = "Zysk firmy, która nie bada konsumentów", ylim = (-1500, 1500))

surplus_per_h_no_research = [sum.(getindex.(ex3_v1_producer_surplus_singleton, 1))[ex3_v1_H .== x] for x in sort(unique(ex3_v1_H))]
surplus_per_h_research = [sum.(getindex.(ex3_v1_producer_surplus_singleton, 2))[ex3_v1_H .== x] for x in sort(unique(ex3_v1_H))]

[count(x .<= 0) / length(x) for x in surplus_per_h_no_research]
[count(x .<= 0) / length(x) for x in surplus_per_h_research]

Plots.savefig(ex1_p10, pwd() * "\\plots\\ex1\\profit per life cycle no research.svg")

ex1_p11 = StatsPlots.boxplot([sum.(getindex.(ex3_v1_producer_surplus_singleton, 2))[ex3_v1_H .== x] for x in sort(unique(ex3_v1_H))], legend = false, xlabel = "Maksymalna przydatność dobra [liczba okresów]", ylabel = "Zysk firmy", title = "Zysk firmy, która bada konsumentów", ylim = (-1500, 1500))

Plots.savefig(ex1_p11, pwd() * "\\plots\\ex1\\profit per life cycle research.svg")

ex2_p1 = Plots.scatter(mean.(getindex.(ex3_v1_price, 1)), sum.(getindex.(ex3_v1_quantity_sold, 1)), smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0, xlabel = "Cena", ylabel = "Wielkość sprzedaży", label = "Producent nie bada konsumentów", xaxis = :log, yaxis = :log)
Plots.scatter!(mean.(getindex.(ex3_v1_price, 2)), sum.(getindex.(ex3_v1_quantity_sold, 2)), smooth = true, markeralpha = 0.25, linewidth = 2, markerstrokewidth = 0, label = "Producent bada konsumentów")

mnb = GLM.lm(@formula(ly~lx), DataFrame(ly = sum.(getindex.(ex3_v1_quantity_sold, 1)), lx = mean.(getindex.(ex3_v1_price, 1))))

mb = GLM.lm(@formula(ly~lx), DataFrame(ly = sum.(getindex.(ex3_v1_quantity_sold, 2)), lx = mean.(getindex.(ex3_v1_price, 2))))

iqr(xydiff.(getindex.(ex3_v1_durability, 1), [getindex.(x,1) for x in ex3_v1_durability_exp]))
cv(xydiff.(getindex.(ex3_v1_durability, 1), [getindex.(x,1) for x in ex3_v1_durability_exp]))
iqr(xydiff.(getindex.(ex3_v1_durability, 2), [getindex.(x,2) for x in ex3_v1_durability_exp]))
cv(xydiff.(getindex.(ex3_v1_durability, 2), [getindex.(x,2) for x in ex3_v1_durability_exp]))


ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v1_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v1_producer_surplus_singleton, 2)) .- 400 * 0.05 * 200 * 0.01, "Producent bada oczekiwania konsumentów")

Plots.savefig(ex3_p1, pwd() * "\\plots\\ex1\\surplus with and without market research.svg")

Plots.plot(trim_first(mean(getindex.(ex3_v1_producer_surplus_singleton, 1)), trimmed = 20))
Plots.plot!(trim_first(mean(getindex.(ex3_v1_producer_surplus_singleton, 2)), trimmed = 20))

mean(trim_first(mean(getindex.(ex3_v1_producer_surplus_singleton, 1)), trimmed = 10))
mean(trim_first(mean(getindex.(ex3_v1_producer_surplus_singleton, 2)), trimmed = 10))

UnequalVarianceTTest(mean.(getindex.(ex3_v1_producer_surplus_singleton, 1)), mean.(getindex.(ex3_v1_producer_surplus_singleton, 2)))

mean(trim_first(mean(getindex.(ex3_v1_producer_surplus_singleton, 1)); trimmed = 5))
mean(trim_first(mean(getindex.(ex3_v1_producer_surplus_singleton, 2)); trimmed = 5))

Plots.plot(trim_first(mean(getindex.(ex3_v1_producer_surplus_singleton, 1)); trimmed = 5))
Plots.plot!(trim_first(mean(getindex.(ex3_v1_producer_surplus_singleton, 2)); trimmed = 5))

ex3_p2 = plot_ecdf(true, mean.(getindex.(ex3_v1_quality, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Średnia jakość", ylabel = "F(x)", title = "Dystrybuanta empiryczna - średnia jakość")
plot_ecdf(false, mean.(getindex.(ex3_v1_quality, 2)), "Producent bada oczekiwania konsumentów")

Plots.savefig(ex3_p2, pwd() * "\\plots\\ex1\\quality with and without market research.svg")

ex3_p3 = plot_ecdf(true, mean.(getindex.(ex3_v1_durability, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Średnia trwałość", ylabel = "F(x)", title = "Dystrybuanta empiryczna - średnia trwałość")
plot_ecdf(false, mean.(getindex.(ex3_v1_durability, 2)), "Producent bada oczekiwania konsumentów")

Plots.savefig(ex3_p3, pwd() * "\\plots\\ex1\\durability with and without market research.svg")

ex3_p4 = plot_ecdf(true, mean.(getindex.(ex3_v1_quantity_produced, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Wielkość produkcji", ylabel = "F(x)", title = "Dystrybuanta empiryczna - wielkość produkcji")
plot_ecdf(false, mean.(getindex.(ex3_v1_quantity_produced, 2)), "Producent bada oczekiwania konsumentów")

Plots.savefig(ex3_p4, pwd() * "\\plots\\ex1\\prod quan with and without market research.svg")

ex3_p5 = plot_ecdf(true, mean.(getindex.(ex3_v1_margin, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Średnia marża", ylabel = "F(x)", title = "Dystrybuanta empiryczna - średnia marża")
plot_ecdf(false, mean.(getindex.(ex3_v1_margin, 2)), "Producent bada oczekiwania konsumentów")

Plots.savefig(ex3_p5, pwd() * "\\plots\\ex1\\prod quan with and without market research.svg")

ex3_p6 = plot_ecdf(true, sum.(getindex.(ex3_v1_quantity_sold, 1)), "
Producent nie bada oczekiwań konsumentów", xlabel = "Wielkość sprzedaży i leasingu", ylabel = "F(x)", title = "Dystrybuanta empiryczna - popyt")
plot_ecdf(false, sum.(getindex.(ex3_v1_quantity_sold, 2)), "Producent bada oczekiwania konsumentów")

Plots.savefig(ex3_p6, pwd() * "\\plots\\ex1\\sell quan with and without market research.svg")