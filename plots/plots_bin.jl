
ex3_p2 = StatsPlots.density(xydiff.(getindex.(ex3_v1_durability, 1), [getindex.(x,1) for x in ex3_v1_durability_exp]), label = "Producent nie prowadzi badań", xlabel = "Odchylenie parametrów produktu [trwałość] od oczekiwań konsumentów", ylabel = "f(x)", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, title = "Funkcja gęstości odchyleń parametrów produktu od oczekiwań konsumentów", legendfontsize = 6, legend = :topleft, normalize = true)
StatsPlots.density!(xydiff.(getindex.(ex3_v1_durability, 2), [getindex.(x,2) for x in ex3_v1_durability_exp]), label = "Producent prowadzi badania")

round(iqr(xydiff.(getindex.(ex3_v1_durability, 1), [getindex.(x,1) for x in ex3_v1_durability_exp])), digits = 3)
round(cv(xydiff.(getindex.(ex3_v1_durability, 1), [getindex.(x,1) for x in ex3_v1_durability_exp])), digits = 3)
round(iqr(xydiff.(getindex.(ex3_v1_durability, 2), [getindex.(x,2) for x in ex3_v1_durability_exp])), digits = 3)
round(cv(xydiff.(getindex.(ex3_v1_durability, 2), [getindex.(x,2) for x in ex3_v1_durability_exp])), digits = 3)

mean(xydiff.(getindex.(ex3_v1_durability, 1), [getindex.(x,1) for x in ex3_v1_durability_exp]))
mean(xydiff.(getindex.(ex3_v1_durability, 2), [getindex.(x,2) for x in ex3_v1_durability_exp]))

Plots.savefig(ex3_p2, pwd() * "\\plots\\ex2\\durability_research_no_research.svg")


StatsPlots.density(mean.(getindex.(ex3_v1_margin, 1)))
StatsPlots.density!(mean.(getindex.(ex3_v1_margin, 2)))

Plots.scatter(mean.(getindex.(ex3_v1_quality, 1)), mean.(getindex.(ex3_v1_producer_surplus_singleton, 1)))

Plots.plot(mean(getindex.(ex3_v1_quality, 1)))
Plots.plot!(mean(getindex.(ex3_v1_quality, 2)))

Plots.plot(mean(getindex.(ex3_v1_durability, 1)))
Plots.plot!(mean(getindex.(ex3_v1_durability, 2)))


Plots.scatter(mean.(getindex.(ex3_v1_quality, 2)), mean.(getindex.(ex3_v1_producer_surplus_singleton, 2)))

ex1_p10 = StatsPlots.boxplot([sum.(getindex.(ex3_v1_producer_surplus_singleton, 1))[ex3_v1_H .== x] for x in sort(unique(ex3_v1_H))], legend = false, xlabel = "Maksymalna przydatność dobra [liczba okresów]", ylabel = "Zysk firmy", title = "Zysk firmy, która nie bada konsumentów", ylim = (-1500, 1500))

surplus_per_h_no_research = [sum.(getindex.(ex3_v1_producer_surplus_singleton, 1))[ex3_v1_H .== x] for x in sort(unique(ex3_v1_H))]
surplus_per_h_research = [sum.(getindex.(ex3_v1_producer_surplus_singleton, 2))[ex3_v1_H .== x] for x in sort(unique(ex3_v1_H))]

round.([count(x .<= 0) / length(x) for x in surplus_per_h_no_research], digits = 3)
round.([count(x .<= 0) / length(x) for x in surplus_per_h_research], digits = 3)

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

close_competitors = abs.(getindex.(getindex.(ex3_v1_quality,1),1) .- getindex.(getindex.(ex3_v1_quality,2),1)) .< 0.05
distant_competitors = abs.(getindex.(getindex.(ex3_v1_quality,1),1) .- getindex.(getindex.(ex3_v1_quality,2),1)) .> 0.20
sum(distant_competitors)