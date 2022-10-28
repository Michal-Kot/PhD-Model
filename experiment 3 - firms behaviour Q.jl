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
ex3_v1_quality_std = []
ex3_v1_durability_std = []

for i in 1:1

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    ex3_v1_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], 0.50, true, 1, [0.7, 1.], "softmax", [false, true], [0, 0.1], 5, false)

    push!(ex3_v1_total_surplus, calculate_surplus(ex3_v1_sim, "total", true))
    push!(ex3_v1_producer_surplus, calculate_surplus(ex3_v1_sim, "producer", true))
    push!(ex3_v1_consumer_surplus, calculate_surplus(ex3_v1_sim, "consumer,total",true))
    push!(ex3_v1_quality, trim_first.(getfield.(ex3_v1_sim.sellers, :quality_history); trimmed = 5))
    push!(ex3_v1_durability, trim_first.(getfield.(ex3_v1_sim.sellers, :durability_history); trimmed = 5))
    push!(ex3_v1_margin, trim_first.(getfield.(ex3_v1_sim.sellers, :margin_history); trimmed = 5))
    push!(ex3_v1_price, trim_first.(calculate_price_history.(ex3_v1_sim.sellers; product_life = 5); trimmed = 5))
    push!(ex3_v1_quantity_produced, trim_first.(getfield.(ex3_v1_sim.sellers, :quantity_produced_history); trimmed = 5))
    push!(ex3_v1_quantity_sold, trim_first.(getfield.(ex3_v1_sim.sellers, :quantity_sold_history); trimmed = 5))
    push!(ex3_v1_producer_surplus_singleton, trim_first.(calculate_profit_history.(ex3_v1_sim.sellers); trimmed = 5))
    push!(ex3_v1_reselling, trim_first.(getfield.(ex3_v1_sim.sellers, :reselling_history); trimmed = 5))
    push!(ex3_v1_buying_history, trim_first.(getfield.(ex3_v1_sim.buyers, :unit_buying_selling_history); trimmed = 5))
    push!(ex3_v1_quality_exp, trim_first.(mean([b.quality_expectation_history for b in ex3_v1_sim.buyers]); trimmed = 5))
    push!(ex3_v1_durability_exp, trim_first.(mean([b.durability_expectation_history for b in ex3_v1_sim.buyers]); trimmed = 5))
    push!(ex3_v1_quality_std, [std(mean.([getindex.(x,1) for x in getfield.(ex3_v1_sim.buyers, :quality_expectation_history)])), std(mean.([getindex.(x,2) for x in getfield.(ex3_v1_sim.buyers, :quality_expectation_history)]))])
    push!(ex3_v1_durability_std, [std(mean.([getindex.(x,1) for x in getfield.(ex3_v1_sim.buyers, :durability_expectation_history)])), std(mean.([getindex.(x,2) for x in getfield.(ex3_v1_sim.buyers, :durability_expectation_history)]))])


end

Plots.scatter(getindex.(ex3_v1_quality_std, 1), mean.(getindex.(ex3_v1_producer_surplus_singleton, 1)), smooth = true)
Plots.scatter!(getindex.(ex3_v1_quality_std, 2), mean.(getindex.(ex3_v1_producer_surplus_singleton, 2)), smooth = true)

StatsPlots.histogram(getindex.(ex3_v1_quality_std, 1), alpha = 0.2, bins = 0:0.01:0.20)
StatsPlots.histogram!(getindex.(ex3_v1_quality_std, 2), alpha = 0.2, bins = 0:0.01:0.20)

Plots.plot(mean(getindex.(ex3_v1_margin,1)))
Plots.plot!(mean(getindex.(ex3_v1_margin,2)))

Plots.plot(mean(getindex.(ex3_v1_quality,1)))
Plots.plot!(mean(getindex.(ex3_v1_quality,2)))

Plots.plot(mean(getindex.(ex3_v1_durability,1)))
Plots.plot!(mean(getindex.(ex3_v1_durability,2)))

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v1_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v1_producer_surplus_singleton, 2)), "Producent bada oczekiwania konsumentów")

sort(sum.(getindex.(ex3_v1_producer_surplus_singleton, 1)))

mean(sum.(getindex.(ex3_v1_producer_surplus_singleton, 1)))
mean(sum.(getindex.(ex3_v1_producer_surplus_singleton, 2)))


ex3_v1_quality_exp[1][1]
getindex.(ex3_v1_quality_exp, 1)
getindex.(ex3_v1_producer_surplus_singleton, 1)

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

ex3_p6 = plot_ecdf(true, sum.(getindex.(ex3_v1_quantity_sold, 1)) .+ sum.(getindex.(ex3_v1_quantity_leased, 1)), "
Producent nie bada oczekiwań konsumentów", xlabel = "Wielkość sprzedaży i leasingu", ylabel = "F(x)", title = "Dystrybuanta empiryczna - popyt")
plot_ecdf(false, sum.(getindex.(ex3_v1_quantity_sold, 2)) .+ sum.(getindex.(ex3_v1_quantity_leased, 2)), "Producent bada oczekiwania konsumentów")

Plots.savefig(ex3_p6, pwd() * "\\plots\\ex1\\sell quan with and without market research.svg")

################################################################################################

ex3_p3 = plot_ecdf(true, mean.(getindex.(ex3_v1_quantity_produced, 1)) .- mean.(getindex.(ex3_v1_quantity_sold, 1)) .- mean.(getindex.(ex3_v1_quantity_leased, 1)), "
Producent bada oczekiwania konsumentów", xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false,  mean.(getindex.(ex3_v1_quantity_produced, 2)) .- mean.(getindex.(ex3_v1_quantity_sold, 2)) .- mean.(getindex.(ex3_v1_quantity_leased, 2)), "Producent nie bada oczekiwań konsumentów")

ρ_cut = cut_integer(float.(ex3_ρ_down), 10)[1]

k  =5

Plots.scatter(mean.(getindex.(ex3_v1_durability, 2)), mean.(getindex.(ex3_v1_producer_surplus_singleton, 2)))

Plots.scatter!(mean.(getindex.(ex3_v2_durability, 1)), mean.(getindex.(ex3_v2_producer_surplus_singleton, 1)))
Plots.scatter!(mean.(getindex.(ex3_v2_durability, 2)), mean.(getindex.(ex3_v2_producer_surplus_singleton, 2)))

StatsPlots.histogram(cut_integer(float.(ex3_ρ_down), 5)[1])

# Istnienie rynku wtórnego podnosi dobrobyt

ex3_p1 = plot_ecdf(true, ex3_v1_total_surplus, "Identyczna jakość, rynek wtórny istnieje", xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false, ex3_v2_total_surplus, "Różna jakość, rynek wtórny istnieje")
plot_ecdf(false, ex3_v3_total_surplus, "Identyczna jakość, rynek wtórny nie istnieje")
plot_ecdf(false, ex3_v4_total_surplus, "Różna jakość, rynek wtórny nie istnieje")

Plots.savefig(ex3_p1, pwd() * "\\plots\\ex3_secmar total surplus.svg")

# ... oraz nadwyżkę konsumenta ...

ex3_p2 = plot_ecdf(true, ex3_v1_consumer_surplus, "Identyczna jakość, rynek wtórny istnieje", xlabel = "Nadwyżka konsumenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka konsumenta")
plot_ecdf(false, ex3_v2_consumer_surplus, "Różna jakość, rynek wtórny istnieje")
plot_ecdf(false, ex3_v3_consumer_surplus, "Identyczna jakość, rynek wtórny nie istnieje")
plot_ecdf(false, ex3_v4_consumer_surplus, "Różna jakość, rynek wtórny nie istnieje")

Plots.savefig(ex3_p2, pwd() * "\\plots\\ex3_secmar cons surplus.svg")

# ... a także nadwyżkę producenta ...

ex3_p3 = plot_ecdf(true, ex3_v1_producer_surplus, "Identyczna jakość, rynek wtórny istnieje", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, ex3_v2_producer_surplus, "Różna jakość, rynek wtórny istnieje")
plot_ecdf(false, ex3_v3_producer_surplus, "Identyczna jakość, rynek wtórny nie istnieje")
plot_ecdf(false, ex3_v4_producer_surplus, "Różna jakość, rynek wtórny nie istnieje")

Plots.savefig(ex3_p3, pwd() * "\\plots\\ex3_secmar prod surplus.svg")

UnequalVarianceTTest(float.(ex3_v1_producer_surplus), float.(ex3_v3_producer_surplus))
UnequalVarianceTTest(float.(ex3_v2_producer_surplus), float.(ex3_v4_producer_surplus))

# Efekt istnienia rynku wtórnego nie jest istotny w przypadku rynku dóbr homogenicznych, ale jest istotny dla rynku, na którym występuje pionowe zróżnicowanie

ex3_p4 = plot_ecdf(true, mean.(getindex.(ex3_v1_producer_surplus_singleton,1)), "Identyczna jakość, rynek wtórny istnieje", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(getindex.(ex3_v2_producer_surplus_singleton,1)), "Różna jakość, rynek wtórny istnieje")
plot_ecdf(false, mean.(getindex.(ex3_v3_producer_surplus_singleton,1)), "Identyczna jakość, rynek wtórny nie istnieje")
plot_ecdf(false, mean.(getindex.(ex3_v4_producer_surplus_singleton,1)), "Różna jakość, rynek wtórny nie istnieje")

Plots.savefig(ex3_p4, pwd() * "\\plots\\ex3_secmar prod better surplus.svg")


UnequalVarianceTTest(float.(mean.(getindex.(ex3_v1_producer_surplus_singleton,1))), float.(mean.(getindex.(ex3_v3_producer_surplus_singleton,1))))
UnequalVarianceTTest(float.(mean.(getindex.(ex3_v2_producer_surplus_singleton,1))), float.(mean.(getindex.(ex3_v4_producer_surplus_singleton,1))))

count_x(y;x) = sum(y .== x)

# Dlaczego istnienie rynku wtórnego podnosi nadwyżkę producenta?

ex3_p5 = StatsPlots.histogram(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex3_v1_buying_history]), alpha = 0.25, normalize = true, xlabel = "Częstotliwość zakupów", ylabel = "f(x)", label = "Rynek wtórny istnieje", title = "Częstotliwość zakupów, identyczna jakość", bins = 0:5:70)
StatsPlots.histogram!(mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex3_v3_buying_history]), alpha = 0.25, normalize = true, label = "Rynek wtórny nie istnieje", bins = 0:5:70, ylim = (0, 0.08))

Plots.savefig(ex3_p5, pwd() * "\\plots\\ex3_primar freq equal.svg")


UnequalVarianceTTest(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex3_v1_buying_history]), mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex3_v3_buying_history]))

ex3_p6 = StatsPlots.histogram(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex3_v2_buying_history]), alpha = 0.25, normalize = true, xlabel = "Częstotliwość zakupów", ylabel = "f(x)", label = "Rynek wtórny istnieje", title = "Częstotliwość zakupów, różna jakość", bins = 0:5:70)
StatsPlots.histogram!(mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex3_v4_buying_history]), alpha = 0.25, normalize = true, label = "Rynek wtórny nie istnieje", bins=0:5:70, ylim = (0, 0.08))

mean(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex3_v2_buying_history]))
mean(mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex3_v4_buying_history]))
UnequalVarianceTTest(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex3_v2_buying_history]), mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex3_v4_buying_history]))

mean(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex3_v1_buying_history])) / 
mean(mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex3_v3_buying_history]))

mean(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex3_v2_buying_history])) / 
mean(mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex3_v4_buying_history]))

Plots.savefig(ex3_p6, pwd() * "\\plots\\ex3_primar freq diff.svg")

# Kto najbardziej zyskuje na rynku wtórnym?

mean(mean(getindex.(ex3_v1_producer_surplus_singleton, 1))) - mean(mean(getindex.(ex3_v3_producer_surplus_singleton, 1)))
UnequalVarianceTTest(mean.(getindex.(ex3_v1_producer_surplus_singleton, 1)), mean.(getindex.(ex3_v3_producer_surplus_singleton, 1)))

mean(mean(getindex.(ex3_v2_producer_surplus_singleton, 1))) - mean(mean(getindex.(ex3_v4_producer_surplus_singleton, 1)))
UnequalVarianceTTest(mean.(getindex.(ex3_v2_producer_surplus_singleton, 1)), mean.(getindex.(ex3_v4_producer_surplus_singleton, 1)))

mean(mean(getindex.(ex3_v1_producer_surplus_singleton, 2))) - mean(mean(getindex.(ex3_v3_producer_surplus_singleton, 2)))
UnequalVarianceTTest(mean.(getindex.(ex3_v1_producer_surplus_singleton, 2)), mean.(getindex.(ex3_v3_producer_surplus_singleton, 2)))

mean(mean(getindex.(ex3_v2_producer_surplus_singleton, 2))) - mean(mean(getindex.(ex3_v4_producer_surplus_singleton, 2)))
UnequalVarianceTTest(mean.(getindex.(ex3_v2_producer_surplus_singleton, 2)), mean.(getindex.(ex3_v4_producer_surplus_singleton, 2)))

UnequalVarianceTTest(sort(mean.(getindex.(ex3_v2_producer_surplus_singleton, 1))) .- sort(mean.(getindex.(ex3_v1_producer_surplus_singleton, 1))), sort(mean.(getindex.(ex3_v4_producer_surplus_singleton, 1))) .- sort(mean.(getindex.(ex3_v3_producer_surplus_singleton, 1))))

UnequalVarianceTTest(sort(mean.(getindex.(ex3_v2_producer_surplus_singleton, 1))) .- sort(mean.(getindex.(ex3_v1_producer_surplus_singleton, 1))), sort(mean.(getindex.(ex3_v2_producer_surplus_singleton, 2))) .- sort(mean.(getindex.(ex3_v4_producer_surplus_singleton, 2))))

plot_ecdf(true, mean.(getindex.(ex3_v3_producer_surplus_singleton, 2)), "Identyczna jakość, rynek wtórny istnieje", xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false, mean.(getindex.(ex3_v4_producer_surplus_singleton, 2)), "Identyczna jakość, rynek wtórny istnieje", xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")


Plots.plot(mean(getindex.(ex3_v3_producer_surplus_singleton, 2)))
# Bo podnosi liczbę produkowanych sztuk 

ex3_p5 = plot_ecdf(true, mean.(mean.(ex3_v1_quantity_produced)), "Identyczna jakość, rynek wtórny istnieje", xlabel = "Liczba produkowanych sztuk", ylabel = "F(x)", title = "Dystrybuanta empiryczna - liczba produkowanych sztuk")
plot_ecdf(false, mean.(mean.(ex3_v2_quantity_produced)), "Różna jakość, rynek wtórny istnieje")
plot_ecdf(false, mean.(mean.(ex3_v3_quantity_produced)), "Identyczna jakość, rynek wtórny nie istnieje")
plot_ecdf(false, mean.(mean.(ex3_v4_quantity_produced)), "Różna jakość, rynek wtórny nie istnieje")

Plots.savefig(ex3_p4, pwd() * "\\plots\\ex3_secmar prod quantity.svg")

# Produkcja jest większa, gdyż i sprzedaż jest większa -> rynek wtórny i odsprzedaż zwiększają frequency zakupów

ex3_p5 = plot_ecdf(true, mean.(mean.(ex3_v1_quantity_sold)), "Identyczna jakość, rynek wtórny istnieje", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(mean.(ex3_v2_quantity_sold)), "Różna jakość, rynek wtórny istnieje")
plot_ecdf(false, mean.(mean.(ex3_v3_quantity_sold)), "Identyczna jakość, rynek wtórny nie istnieje")
plot_ecdf(false, mean.(mean.(ex3_v4_quantity_sold)), "Różna jakość, rynek wtórny nie istnieje")

ex3_p4 = plot_ecdf(true, mean.(mean.(ex3_v1_quantity_leased)), "Identyczna jakość, rynek wtórny istnieje", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(mean.(ex3_v2_quantity_leased)), "Różna jakość, rynek wtórny istnieje")
plot_ecdf(false, mean.(mean.(ex3_v3_quantity_leased)), "Identyczna jakość, rynek wtórny nie istnieje")
plot_ecdf(false, mean.(mean.(ex3_v4_quantity_leased)), "Różna jakość, rynek wtórny nie istnieje")

Plots.scatter(sum.(sum.(ex3_v1_quantity_produced)), sum.(sum.(ex3_v1_reselling)), smooth = true, markerstrokewidth = 0, markeralpha = 0.25, linewidth = 3, xlabel = "Liczba wyprodukowanych dóbr", ylabel = "Liczba odsprzedanych dóbr", label = "Identyczna jakość", legend = :bottomright)
Plots.scatter!(sum.(sum.(ex3_v2_quantity_produced)), sum.(sum.(ex3_v2_reselling)), smooth = true, markerstrokewidth = 0, markeralpha = 0.25, linewidth = 3, label = "Różna jakość")

Plots.scatter(sum.(getindex.(ex3_v1_quantity_produced, 1)), sum.(getindex.(ex3_v1_reselling, 1)), smooth = true, markerstrokewidth = 0, markeralpha = 0.25, linewidth = 3, label = "Identyczna jakość", ylabel = "Liczba odsprzedanych dóbr", xlabel = "Liczba wyprodukowanych dóbr")
Plots.scatter!(sum.(getindex.(ex3_v2_quantity_produced, 1)), sum.(getindex.(ex3_v2_reselling, 1)), smooth = true, markerstrokewidth = 0, markeralpha = 0.25, linewidth = 3, label = "Różna jakość, lepszy produkt")
Plots.scatter!(sum.(getindex.(ex3_v2_quantity_produced, 2)), sum.(getindex.(ex3_v2_reselling, 2)), smooth = true, markerstrokewidth = 0, markeralpha = 0.25, linewidth = 3, label = "Różna jakość, gorszy produkt")

GLM.lm(@formula(y~x), DataFrame(x= sum.(getindex.(ex3_v1_quantity_produced, 1)), y=sum.(getindex.(ex3_v1_reselling, 1)) ))
GLM.lm(@formula(y~x), DataFrame(x= sum.(getindex.(ex3_v2_quantity_produced, 1)), y=sum.(getindex.(ex3_v2_reselling, 1)) ))
GLM.lm(@formula(y~x), DataFrame(x= sum.(getindex.(ex3_v2_quantity_produced, 2)), y=sum.(getindex.(ex3_v2_reselling, 2)) ))


Plots.scatter(mean.(getindex.(ex3_v1_margin, 1)) .- 1, mean.(getindex.(ex3_v1_producer_surplus_singleton, 1)), markeralpha = 0.25, markerstrokewidth = 0, label = "Identyczna jakość", xlabel = "Marżowość", ylabel = "Zysk")
add_smoothing_spline(mean.(getindex.(ex3_v1_margin, 1)) .- 1, mean.(getindex.(ex3_v1_producer_surplus_singleton, 1)), "blue")

Plots.scatter(mean.(getindex.(ex3_v2_margin, 1)) .- 1, mean.(getindex.(ex3_v2_producer_surplus_singleton, 1)), markeralpha = 0.25, markerstrokewidth = 0, label = "Różna jakość, producent 1", xlabel = "Marża %", ylabel = "Zysk")
add_smoothing_spline(mean.(getindex.(ex3_v2_margin, 1)) .- 1, mean.(getindex.(ex3_v2_producer_surplus_singleton, 1)), "blue")

Plots.scatter(mean.(getindex.(ex3_v2_margin, 2)) .- 1, mean.(getindex.(ex3_v2_producer_surplus_singleton, 2)), markeralpha = 0.25, markerstrokewidth = 0, label = "Różna jakość, producent 2", xlabel = "Marża %", ylabel = "Zysk")
add_smoothing_spline(mean.(getindex.(ex3_v2_margin, 2)) .- 1, mean.(getindex.(ex3_v2_producer_surplus_singleton, 2)), "blue")

Plots.plot(xlabel = "Marża %", ylabel = "Zysk")
add_smoothing_spline(mean.(getindex.(ex3_v1_margin, 1)) .- 1, mean.(getindex.(ex3_v1_producer_surplus_singleton, 1)), "red", "Identyczna jakość")
add_smoothing_spline(mean.(getindex.(ex3_v2_margin, 1)) .- 1, mean.(getindex.(ex3_v2_producer_surplus_singleton, 1)), "green", "Różna jakość, producent 1")
add_smoothing_spline(mean.(getindex.(ex3_v2_margin, 2)) .- 1, mean.(getindex.(ex3_v2_producer_surplus_singleton, 2)), "blue", "Różna jakość, producent 2")

Plots.scatter([getindex.([mean.(p) for p in ex3_v2_margin],x) for x in 1:2], [getindex.([mean.(p) for p in ex3_v2_producer_surplus_singleton],x) for x in 1:2], xlim = (0.8,2.0), ylim = (-5,30))

plot([mean(sum.(getindex.(ex3_v1_producer_surplus_singleton, 1))[ex3_in .== x]) for x in sort(unique(ex3_in_d))])
plot([mean(sum.(getindex.(ex3_v2_producer_surplus_singleton, 1))[ex3_in_d .== x]) for x in sort(unique(ex3_in_d))])

scatter(getindex.(ex3_in, 1), sum.(getindex.(ex3_v1_quantity_leased, 1)))
scatter(ex3_in_d, sum.(getindex.(ex3_v2_producer_surplus_singleton, 1)))

ex2_p1 = plot_ecdf(ex3_v1_total_surplus, "Identyczna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", true)
plot_ecdf(ex3_v2_total_surplus, "Różna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)

ex2_p3 = plot_ecdf(ex3_v1_producer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex3_v2_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex2_p3 = plot_ecdf(ex3_v1_consumer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex3_v2_consumer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex2_p3 = plot_ecdf(sum.(sum.(ex3_v1_quantity_produced)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_produced)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(sum.(sum.(ex3_v1_quantity_sold)) .+ sum.(sum.(ex3_v1_quantity_leased)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_sold)) .+ sum.(sum.(ex3_v2_quantity_leased)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(sum.(sum.(ex3_v1_quantity_sold)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_sold)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(sum.(sum.(ex3_v1_quantity_leased)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_leased)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(sum.(sum.(getindex.(ex3_v1_producer_surplus_singleton,1))), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(getindex.(ex3_v2_producer_surplus_singleton,1))), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

L1 = getindex.(ex3_v12_r, 1) 
L2 = getindex.(ex3_v12_r, 2)

R1 = getindex.(ex3_v12_i, 1)
R2 = getindex.(ex3_v12_i, 2)


hm_eq = [mean_c(mean.(getindex.(ex3_v1_producer_surplus_singleton, 1))[(L1 .== l1) .& (R1 .== r1)]) for l1 in sort(unique(L1)), r1 in sort(unique(R1))]'
hm_dq = [mean_c(mean.(getindex.(ex3_v2_producer_surplus_singleton, 1))[(L1 .== l1) .& (R1 .== r1)]) for l1 in sort(unique(L1)), r1 in sort(unique(R1))]'

ex3_p1 = heatmap(sort(unique(L1)), sort(unique(R1)), hm_eq, xlabel = "Decrease rate", ylabel = "Increase rate", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość początkowa", titlefontsize = 8)
ex3_p1 = heatmap(sort(unique(L1)), sort(unique(R1)), hm_dq, xlabel = "Decrease rate", ylabel = "Increase rate", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość początkowa", titlefontsize = 8)

####

hm_eq = [mean_c(mean.(getindex.(ex3_v1_producer_surplus_singleton, 1))[(R1 .== r1) .& (R2 .== r2)]) for r1 in sort(unique(R1)), r2 in sort(unique(R2))]'
hm_dq = [mean_c(mean.(getindex.(ex3_v2_producer_surplus_singleton, 1))[(R1 .== r1) .& (R2 .== r2)]) for r1 in sort(unique(R1)), r2 in sort(unique(R2))]'

max_plot = maximum([maximum(hm_eq), maximum(hm_dq)])
min_plot = minimum([minimum(hm_eq), minimum(hm_dq)])

ex3_p1 = heatmap(sort(unique(R1)), sort(unique(R2)), hm_eq[2:end,2:end], xlabel = "Decrease rate, my", ylabel = "Decrease rate, competitors", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość początkowa", titlefontsize = 8, clim = (min_plot, max_plot))

ex3_p1 = heatmap(sort(unique(R1)), sort(unique(R2)), hm_dq[2:end,2:end], xlabel = "Decrease rate, my", ylabel = "Decrease rate, competitors", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość początkowa", titlefontsize = 8, clim = (min_plot, max_plot))

####


dR = R1 .- R2 # increase
dR = round.(dR, digits=2)
plot(sort(unique(dR)), [mean(mean.(getindex.(ex3_v1_producer_surplus_singleton, 1))[dR .== dR_u]) for dR_u in sort(unique(dR))])
plot!(sort(unique(dR)), [mean(mean.(getindex.(ex3_v2_producer_surplus_singleton, 1))[dR .== dR_u]) for dR_u in sort(unique(dR))])

dL = L1 .- L2 # decrease
dL = round.(dL, digits=2)
plot(sort(unique(dL)), [mean(mean.(getindex.(ex3_v1_producer_surplus_singleton, 1))[dL .== dL_u]) for dL_u in sort(unique(dL))])
plot!(sort(unique(dL)), [mean(mean.(getindex.(ex3_v2_producer_surplus_singleton, 1))[dL .== dL_u]) for dL_u in sort(unique(dL))])


ex3_p1 = plot_ecdf(ex3_v1_total_surplus, "Identyczna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", true)
plot_ecdf(ex3_v2_total_surplus, "Różna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)

ex3_p2 = plot_ecdf(ex3_v1_consumer_surplus, "equal Q", "Consumer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex3_v2_consumer_surplus, "not equal Q", "Consumer Surplus", "Probability", "ECDF", false)

ex3_p3 = plot_ecdf(ex3_v1_producer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex3_v2_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex3_p3 = plot_ecdf(sum.(sum.(ex3_v1_quantity_produced)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_produced)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex3_p3 = plot_ecdf(sum.(sum.(ex3_v1_quantity_sold)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_sold)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)



mean.(mean.(ex3_v1_price))

ex3_p3 = plot_ecdf(mean.(mean.(ex3_v1_price)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_sold)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)