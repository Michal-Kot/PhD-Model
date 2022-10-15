# Experiment 5

include(pwd() * "\\methods\\methods.jl")

# Impact of secondary market on consumer beliefs

# Dlaczego wyższe durability powoduje spadek consumer surplus?

ex5_v1_total_surplus = []
ex5_v1_producer_surplus = []
ex5_v1_consumer_surplus = []
ex5_v1_price = []
ex5_v1_quality = []
ex5_v1_durability = []
ex5_v1_margin = []
ex5_v1_quantity_produced = []
ex5_v1_quantity_sold = []
ex5_v1_quantity_leased = []
ex5_v1_producer_surplus_singleton = []
ex5_v1_reselling = []
ex5_v1_consumer_surplus_pm = []
ex5_v1_consumer_surplus_sm_s = []
ex5_v1_consumer_surplus_sm_b = []
ex5_v1_quality_exp = []
ex5_v1_durability_exp = []

ex5_v2_total_surplus = []
ex5_v2_producer_surplus = []
ex5_v2_consumer_surplus = []
ex5_v2_price = []
ex5_v2_quality = []
ex5_v2_durability = []
ex5_v2_margin = []
ex5_v2_quantity_produced = []
ex5_v2_quantity_sold = []
ex5_v2_quantity_leased = []
ex5_v2_producer_surplus_singleton = []
ex5_v2_reselling = []
ex5_v2_consumer_surplus_pm = []
ex5_v2_consumer_surplus_sm_s = []
ex5_v2_consumer_surplus_sm_b = []
ex5_v2_quality_exp = []
ex5_v2_durability_exp = []

for i in 1:250

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    ex5_v1_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.5, 0.95], [0.5, 0.95]], [[0.8, 2.], [0.8, 2.]], 0.1, true, true, 1, [0.7, 1.0], "softmax", [true, false])

    push!(ex5_v1_total_surplus, calculate_surplus(ex5_v1_sim, "total", false))
    push!(ex5_v1_producer_surplus, calculate_surplus(ex5_v1_sim, "producer", false))
    push!(ex5_v1_consumer_surplus, calculate_surplus(ex5_v1_sim, "consumer,total", false))
    push!(ex5_v1_price, calculate_price_history.(ex5_v1_sim.sellers))
    push!(ex5_v1_quality, getfield.(ex5_v1_sim.sellers, :quality_history))
    push!(ex5_v1_durability, getfield.(ex5_v1_sim.sellers, :durability_history))
    push!(ex5_v1_margin, getfield.(ex5_v1_sim.sellers, :margin_history))
    push!(ex5_v1_quantity_produced, getfield.(ex5_v1_sim.sellers, :quantity_produced_history))
    push!(ex5_v1_quantity_sold, getfield.(ex5_v1_sim.sellers, :quantity_sold_history))
    push!(ex5_v1_quantity_leased, getfield.(ex5_v1_sim.sellers, :quantity_leased_history))
    push!(ex5_v1_producer_surplus_singleton, calculate_profit_history.(ex5_v1_sim.sellers))
    push!(ex5_v1_reselling, getfield.(ex5_v1_sim.sellers, :reselling_history))
    push!(ex5_v1_consumer_surplus_pm, calculate_surplus(ex5_v1_sim, "consumer,pm", false))
    push!(ex5_v1_consumer_surplus_sm_b, calculate_surplus(ex5_v1_sim, "consumer,sm,b", false))
    push!(ex5_v1_consumer_surplus_sm_s, calculate_surplus(ex5_v1_sim, "consumer,sm,s", false))
    push!(ex5_v1_quality_exp, mean([b.quality_expectation_history for b in ex5_v1_sim.buyers]))
    push!(ex5_v1_durability_exp, mean([b.durability_expectation_history for b in ex5_v1_sim.buyers]))

    ex5_v2_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.5, 0.95], [0.5, 0.95]], [[0.8, 2.], [0.8, 2.]], 0.1, false, true, 1, [0.7, 1.0], "softmax", [true, false])

    push!(ex5_v2_total_surplus, calculate_surplus(ex5_v2_sim, "total", false))
    push!(ex5_v2_producer_surplus, calculate_surplus(ex5_v2_sim, "producer", false))
    push!(ex5_v2_consumer_surplus, calculate_surplus(ex5_v2_sim, "consumer,total",false))
    push!(ex5_v2_price, calculate_price_history.(ex5_v2_sim.sellers))
    push!(ex5_v2_quality, getfield.(ex5_v2_sim.sellers, :quality_history))
    push!(ex5_v2_durability, getfield.(ex5_v2_sim.sellers, :durability_history))
    push!(ex5_v2_margin, getfield.(ex5_v2_sim.sellers, :margin_history))
    push!(ex5_v2_quantity_produced, getfield.(ex5_v2_sim.sellers, :quantity_produced_history))
    push!(ex5_v2_quantity_sold, getfield.(ex5_v2_sim.sellers, :quantity_sold_history))
    push!(ex5_v2_quantity_leased, getfield.(ex5_v2_sim.sellers, :quantity_leased_history))
    push!(ex5_v2_producer_surplus_singleton, calculate_profit_history.(ex5_v2_sim.sellers))
    push!(ex5_v2_reselling, getfield.(ex5_v2_sim.sellers, :reselling_history))
    push!(ex5_v2_consumer_surplus_pm, calculate_surplus(ex5_v2_sim, "consumer,pm", false))
    push!(ex5_v2_consumer_surplus_sm_b, calculate_surplus(ex5_v2_sim, "consumer,sm,b", false))
    push!(ex5_v2_consumer_surplus_sm_s, calculate_surplus(ex5_v2_sim, "consumer,sm,s", false))
    push!(ex5_v2_quality_exp, mean([b.quality_expectation_history for b in ex5_v2_sim.buyers]))
    push!(ex5_v2_durability_exp, mean([b.durability_expectation_history for b in ex5_v2_sim.buyers]))

end

# Różnice w jakości nie mają istotnego wpływu na kształtowanie się trwałości produktów. Leasing ma istotny wpływ

ex5_p1 = plot_ecdf(true, mean.(mean.(ex5_v1_durability)), "Identyczna jakość, leasing możliwy"; xlabel = "Średnia trwałość dóbr", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Trwałość")
plot_ecdf(false, mean.(mean.(ex5_v2_durability)), "Różna jakość, leasing możliwy")
plot_ecdf(false, mean.(mean.(ex5_v3_durability)), "Identyczna jakość, leasing niemożliwy")
plot_ecdf(false, mean.(mean.(ex5_v4_durability)), "Różna jakość, leasing niemożliwy")

UnequalVarianceTTest(mean.(mean.(ex5_v1_durability)), mean.(mean.(ex5_v3_durability)))

UnequalVarianceTTest(mean.(mean.(ex5_v2_durability)), mean.(mean.(ex5_v4_durability)))

Plots.savefig(ex5_p1, pwd() * "\\plots\\ex2_ecdf durability.svg")



# Ma to zastosowanie do średniej rynkowej trwałości, jak również trwałości wszystkich graczy rynkowych

plot_ecdf(true, mean.(getindex.(ex5_v1_durability, 1)), "Identyczna jakość, leasing możliwy"; xlabel = "Średnia trwałość dóbr", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false, mean.(getindex.(ex5_v2_durability, 1)), "Różna jakość, leasing możliwy")
plot_ecdf(false, mean.(getindex.(ex5_v3_durability, 1)), "Identyczna jakość, leasing niemożliwy")
plot_ecdf(false, mean.(getindex.(ex5_v4_durability, 1)), "Różna jakość, leasing niemożliwy")

# Różnica w poziomie średniej rynkowej trwałości pomiędzy rynkiem z leasingiem i bez leasingu, jest wyższa na rynku gdzie dobra są homogeniczne

UnequalVarianceTTest(sort(mean.(getindex.(ex5_v1_durability, 1))) .- sort(mean.(getindex.(ex5_v3_durability, 1))), sort(mean.(getindex.(ex5_v2_durability, 1))) .- sort(mean.(getindex.(ex5_v4_durability, 1))))

mean(sort(mean.(getindex.(ex5_v1_durability, 1))) .- sort(mean.(getindex.(ex5_v3_durability, 1))))
mean(sort(mean.(getindex.(ex5_v2_durability, 1))) .- sort(mean.(getindex.(ex5_v4_durability, 1))))

# Leasing nie ma istotnego wpływu na różnice w decyzjach odnośnie jakości produktów

ex5_p1 = plot_ecdf(true, mean.(mean.(ex5_v1_quality)), "Identyczna jakość, leasing możliwy"; xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false, mean.(mean.(ex5_v2_quality)), "Różna jakość, leasing możliwy")
plot_ecdf(false, mean.(mean.(ex5_v3_quality)), "Identyczna jakość, leasing niemożliwy")
plot_ecdf(false, mean.(mean.(ex5_v4_quality)), "Różna jakość, leasing niemożliwy")

UnequalVarianceTTest(mean.(mean.(ex5_v1_quality)), mean.(mean.(ex5_v3_quality)))
UnequalVarianceTTest(mean.(mean.(ex5_v2_quality)), mean.(mean.(ex5_v4_quality)))

# Na rynku gdzie możliwy jest leasing produkowana jest wyższa liczba produktów

ex5_p1 = plot_ecdf(true, mean.(mean.(ex5_v1_quantity_produced)), "Identyczna jakość, leasing możliwy", xlabel = "Średni wolumen produkcji", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Wielkość produkcji")
plot_ecdf(false, mean.(mean.(ex5_v2_quantity_produced)), "Różna jakość, leasing możliwy")
plot_ecdf(false, mean.(mean.(ex5_v3_quantity_produced)), "Identyczna jakość, leasing niemożliwy")
plot_ecdf(false, mean.(mean.(ex5_v4_quantity_produced)), "Różna jakość, leasing niemożliwy")

UnequalVarianceTTest(mean.(mean.(ex5_v1_quantity_produced)), mean.(mean.(ex5_v1_quantity_produced)))
UnequalVarianceTTest(mean.(mean.(ex5_v1_quantity_produced)), mean.(mean.(ex5_v1_quantity_produced)))

ex5_p1 = plot_ecdf(true, mean.(mean.(ex5_v1_quantity_produced)) .- mean.(mean.(ex5_v1_quantity_sold)) .- mean.(mean.(ex5_v1_quantity_leased)), "Identyczna jakość, leasing możliwy", xlabel = "Liczba niesprzedanych sztuk", ylabel = "F(x)", title = "Dystrybuanta empiryczna - liczba niesprzedanych sztuk")
plot_ecdf(false, mean.(mean.(ex5_v2_quantity_produced)).- mean.(mean.(ex5_v2_quantity_sold)) .- mean.(mean.(ex5_v2_quantity_leased)), "Różna jakość, leasing możliwy")
plot_ecdf(false, mean.(mean.(ex5_v3_quantity_produced)).- mean.(mean.(ex5_v3_quantity_sold)) .- mean.(mean.(ex5_v3_quantity_leased)), "Identyczna jakość, leasing niemożliwy")
plot_ecdf(false, mean.(mean.(ex5_v4_quantity_produced)).- mean.(mean.(ex5_v4_quantity_sold)) .- mean.(mean.(ex5_v4_quantity_leased)), "Różna jakość, leasing niemożliwy")

# W rezultacie, leasing prowadzi do wyższej nadwyżki producenta

ex5_p1 = plot_ecdf(true, mean.(mean.(getindex.(ex5_v1_producer_surplus_singleton, 1))), "Identyczna jakość, leasing możliwy"; xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(mean.(getindex.(ex5_v2_producer_surplus_singleton, 1))), "Różna jakość, leasing możliwy")
plot_ecdf(false, mean.(mean.(getindex.(ex5_v3_producer_surplus_singleton, 1))), "Identyczna jakość, leasing niemożliwy")
plot_ecdf(false, mean.(mean.(getindex.(ex5_v4_producer_surplus_singleton, 1))), "Różna jakość, leasing niemożliwy")

# Wyższy % leasingowanych produktów, prowadzi do wyższych trwałości
# Czy odwrotnie? Wyższe trwałości prowadzą konsumentów do chętniejszego wybierania leasingu
# Czy to nie jest zaprzeczenie teorii Coase-Swana? Albo uszczegółowienie o kierunek przyczynowości?

ex5_v1_perc_leased = mean.(getindex.(ex5_v1_quantity_leased, 1)) ./ mean.(getindex.(ex5_v1_quantity_produced, 1))
ex5_v2_perc_leased = mean.(getindex.(ex5_v2_quantity_leased, 1)) ./ mean.(getindex.(ex5_v2_quantity_produced, 1))
ex5_v3_perc_leased = mean.(getindex.(ex5_v3_quantity_leased, 1)) ./ mean.(getindex.(ex5_v3_quantity_produced, 1))
ex5_v4_perc_leased = mean.(getindex.(ex5_v4_quantity_leased, 1)) ./ mean.(getindex.(ex5_v4_quantity_produced, 1))

ex5_p2 = Plots.scatter(ex5_v1_perc_leased[mean.(getindex.(ex5_v1_durability,1)) .> 0.32], mean.(getindex.(ex5_v1_durability,1))[mean.(getindex.(ex5_v1_durability,1)) .> 0.32], smooth = true, xlabel = "Udział produktów leasingowanych w całości produkcji", ylabel = "Średnia trwałość", title = "Trwałość produktów a udział leasingu", label = "Identyczna jakość", alpha = 0.25, markerstrokewidth = 0, linewidth = 3, linecolor = "blue")
Plots.scatter!(ex5_v2_perc_leased[mean.(getindex.(ex5_v2_durability,1)) .> 0.42], mean.(getindex.(ex5_v2_durability,1))[mean.(getindex.(ex5_v2_durability,1)) .> 0.42], smooth = true, label = "Różna jakość", alpha = 0.25, markerstrokewidth = 0, linewidth = 3, linecolor = "red")

Plots.scatter(mean.(getindex.(ex5_v1_durability, 2)), mean.(getindex.(ex5_v1_producer_surplus_singleton, 2)), smooth = true)

Plots.scatter(mean.(getindex.(ex5_v2_durability, 2)), mean.(getindex.(ex5_v2_producer_surplus_singleton, 2)), smooth = true)

Plots.scatter!(mean.(getindex.(ex5_v3_durability, 2)), mean.(getindex.(ex5_v3_producer_surplus_singleton, 2)), smooth = true)

Plots.scatter!(mean.(getindex.(ex5_v4_durability, 2)), mean.(getindex.(ex5_v4_producer_surplus_singleton, 2)), smooth = true)

Plots.scatter!(mean.(getindex.(ex5_v1_durability, 1)), mean.(getindex.(ex5_v1_producer_surplus_singleton, 1)), smooth = true)


Plots.scatter!(mean.(getindex.(ex5_v2_durability, 1)), mean.(getindex.(ex5_v2_producer_surplus_singleton, 1)), smooth = true)

Plots.scatter!(mean.(getindex.(ex5_v2_durability, 1)), mean.(getindex.(ex5_v2_producer_surplus_singleton, 1)))

Plots.savefig(ex5_p2, pwd() * "\\plots\\ex2_leasing vs durability.svg")

# WAŻNE: Im wyższa jakość produktu, tym silnejszy wzrost nadwyżki producenta wraz ze wzrostem trwałości produktu. Producentom dóbr wysokiej jakości opłacalne jest dbanie o ich wysoką trwałość.

ex5_p3 = Plots.Plots.scatter(mean.(getindex.(ex5_v1_durability,1)), mean.(getindex.(ex5_v1_producer_surplus_singleton,1)), markeralpha = 0.05, label = nothing, markerstrokewidth = 0, xlabel = "Trwałość", ylabel = "Nadwyżka producenta", title = "Nadwyżka producenta a trwałość")
plot_lm(mean.(getindex.(ex5_v1_durability,1)), mean.(getindex.(ex5_v1_producer_surplus_singleton,1)), "Identyczna jakość, leasing możliwy", "blue")

Plots.scatter!(mean.(getindex.(ex5_v2_durability,1)), mean.(getindex.(ex5_v2_producer_surplus_singleton,1)), markeralpha = 0.05, label = nothing, markerstrokewidth = 0)
plot_lm(mean.(getindex.(ex5_v2_durability,1)), mean.(getindex.(ex5_v2_producer_surplus_singleton,1)), "Różna jakość, leasing możliwy", "green")

Plots.scatter!(mean.(getindex.(ex5_v3_durability,1)), mean.(getindex.(ex5_v3_producer_surplus_singleton,1)), markeralpha = 0.05, label = nothing, markerstrokewidth = 0)
plot_lm(mean.(getindex.(ex5_v3_durability,1)), mean.(getindex.(ex5_v3_producer_surplus_singleton,1)), "Identyczna jakość, leasing niemożliwy", "orange")

Plots.scatter!(mean.(getindex.(ex5_v4_durability,1)), mean.(getindex.(ex5_v4_producer_surplus_singleton,1)), markeralpha = 0.05, label = nothing, markerstrokewidth = 0)
plot_lm(mean.(getindex.(ex5_v4_durability,1)), mean.(getindex.(ex5_v4_producer_surplus_singleton,1)), "Różna jakość, leasing niemożliwy", "red")

Plots.savefig(ex5_p3, pwd() * "\\plots\\ex5 optimal dec for best quality.svg")

Plots.scatter(mean.(getindex.(ex5_v1_durability,2)), mean.(getindex.(ex5_v1_producer_surplus_singleton,2)), markeralpha = 0.25, label = nothing, markerstrokewidth = 0)


Plots.scatter!(mean.(getindex.(ex5_v2_durability,2)), mean.(getindex.(ex5_v2_producer_surplus_singleton,2)), markeralpha = 0.25, label = nothing, markerstrokewidth = 0)
plot_lm(mean.(getindex.(ex5_v2_durability,2)), mean.(getindex.(ex5_v2_producer_surplus_singleton,2)), "Różna jakość, leasing możliwy", "green")

Plots.scatter!(mean.(getindex.(ex5_v3_durability,2)), mean.(getindex.(ex5_v3_producer_surplus_singleton,2)), markeralpha = 0.25, label = nothing, markerstrokewidth = 0)
plot_lm(mean.(getindex.(ex5_v3_durability,2)), mean.(getindex.(ex5_v3_producer_surplus_singleton,2)), "Identyczna jakość, leasing niemożliwy", "orange")

Plots.scatter!(mean.(getindex.(ex5_v4_durability,2)), mean.(getindex.(ex5_v4_producer_surplus_singleton,2)), markeralpha = 0.25, label = "Różna jakość, leasing niemożliwy", markerstrokewidth = 0)
plot_lm(mean.(getindex.(ex5_v1_durability,2)), mean.(getindex.(ex5_v1_producer_surplus_singleton,2)), "Identyczna jakość, leasing możliwy", "blue")

function plot_lm(x,y,lbl,col)

    model_lm = GLM.lm(@formula(y~x), DataFrame(y = y, x = x))
    y_hat = predict(model_lm, DataFrame(y = y, x = x))
    Plots.plot!(sort(x),y_hat[sortperm(x)],label = lbl, color = col)

end

plot_lm(mean.(getindex.(ex5_v1_durability,2)), mean.(getindex.(ex5_v1_producer_surplus_singleton,2)))

Plots.scatter(mean.(getindex.(ex5_v1_durability,2)), predict(ex5_m1, DataFrame(y = mean.(getindex.(ex5_v1_producer_surplus_singleton,2)), x = mean.(getindex.(ex5_v1_durability,2)))))

Plots.plot(ex5_m1.fitted)

# Wyższa jakość dóbr podnosi poziom dobrobytu. Zgodne z wynikiem Coase.

Plots.scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_total_surplus)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Identyczna jakość, leasing możliwy")
Plots.scatter!(mean.(mean.(ex5_v2_durability)), sum.(sum.(ex5_v2_total_surplus)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Różna jakość, leasing możliwy")
Plots.scatter!(mean.(mean.(ex5_v3_durability)), sum.(sum.(ex5_v3_total_surplus)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Identyczna jakość, leasing niemożliwy")
Plots.scatter!(mean.(mean.(ex5_v4_durability)), sum.(sum.(ex5_v4_total_surplus)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Różna jakość, leasing niemożliwy")

# Nie istnieje statystycznie istotna zależność pomiędzy trwałością a nadwyżką konsumenta. 
# Do sprawdzenia - czy ma na to wpływ c? 🍎

Plots.scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_consumer_surplus)))
Plots.scatter!(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_producer_surplus)))

# Ma za to wpływ na nadwyżkę na rynku wtórnym, co wydaje się poprawne.

Plots.scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_consumer_surplus_pm)), smooth=true)
Plots.scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_consumer_surplus_sm_b)), smooth=true)
Plots.scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_consumer_surplus_sm_s)), smooth=true)

Plots.scatter(mean.(mean.(ex5_v1_durability)), mean.(mean.(ex5_v1_price)))

plot(sort(0.5 * mean.(mean.(ex5_v1_margin))))

Plots.scatter!(mean.(mean.(ex5_v1_durability)), 0.5 * sum_of_geom_series.(mean.(mean.(ex5_v1_quality)), mean.(mean.(ex5_v1_durability))))

Plots.scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_quantity_sold)), smooth = true)
Plots.scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_quantity_leased)))
