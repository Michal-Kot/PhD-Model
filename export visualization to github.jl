include(pwd() * "\\methods\\methods.jl")

using Latexify

@load "C:\\Users\\User\\Documents\\PhDWorkspace_more_nl.jld2"

iterSum = lastindex(ex3_v303_producer_surplus_singleton)

p_h0_1 = Plots.plot(cumsum(mean.(getindex.(ex3_v303_quality, 1))) ./ collect(1:iterSum), xlabel = "Liczba symulacji", ylabel = "Średnia wartość miary dla \n danej liczby ukończonych symulacji", legend = false, title = "Średnia jakość produktu, \n producent pierwszy")
p_h0_2 = Plots.plot(cumsum(sum.(getindex.(ex3_v303_durability, 1))) ./ collect(1:iterSum), xlabel = "Liczba symulacji", ylabel = "Średnia wartość miary dla \n danej liczby ukończonych symulacji", legend = false, title = "Średnia trwałość produktu, \n producent pierwszy")
p_h0_3 = Plots.plot(cumsum(sum.(getindex.(ex3_v303_quantity_sold, 1))) ./ collect(1:iterSum), xlabel = "Liczba symulacji", ylabel = "Średnia wartość miary dla \n danej liczby ukończonych symulacji", legend = false, title = "Średnia liczba sprzedanych produktów, \n producent pierwszy")
p_h0_4 = Plots.plot(cumsum(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1))) ./ collect(1:iterSum), xlabel = "Liczba symulacji", ylabel = "Średnia wartość miary dla \n danej liczby ukończonych symulacji", legend = false, title = "Średni zysk, \n producent pierwszy")

p_h0_1234 = Plots.plot(p_h0_1, p_h0_2, p_h0_3, p_h0_4, layout = (2,2), size = (1000, 1000))

Plots.savefig(p_h0_1234, pwd() * "\\plots_to_export\\results_stability_means.pdf")

############################################################################################
######################################## HISTORIA ##########################################

####################################### Hipoteza 2 #########################################

"""
Ryzyko poniesienia przez producenta straty jest uzależnione od uwzględnienia preferencji konsumentów podczas podejmowania decyzji oraz od długości okresu przydatności produktów oferowanych na rynku
"""

# Czy zyski producentów prowadzących badania i nie się różnią?

p_h2_1 = plot_ecdf(true, sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)), "Żaden producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta pierwszego", ylabel = "Dystrybuanta empiryczna", xlim = (-300,1200))
plot_ecdf(false, sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)), "Producent pierwszy bada oczekiwania konsumentów jako jedyny")
plot_ecdf(false, sum.(getindex.(ex3_v302_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów jako jedyny")
plot_ecdf(false, sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)), "Obaj producenci badają oczekiwania konsumentów")

# Barrett test for s-th order SD, H0: s1 SD s-th s2, H1: ~H0

pysdtest.test_sd(sample1 = sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)), 
    sample2 = sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)), 
    ngrid = 50, s = 1, resampling = "bootstrap").testing() 
pysdtest.test_sd(sample1 = sum.(getindex.(ex3_v302_producer_surplus_singleton, 1)), 
    sample2 = sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)), 
    ngrid = 1000, s = 1, resampling = "bootstrap").testing() 
pysdtest.test_sd_SR(sample1 = sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)), 
    sample2 = sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)), 
    ngrid = 500, s = 1, resampling = "bootstrap").testing() 

latexify(
    DataFrame(Typ = [L"Żaden producent nie bada oczekiwań konsumentów", L"Producent pierwszy bada oczekiwania konsumentów jako jedyny", L"Producent drugi bada oczekiwania konsumentów jako jedyny", L"Obaj producenci badają oczekiwania konsumentów"],
    Sredni_zysk = round.([mean(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))), 
    mean(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))),
    mean(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))),
    mean(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)))], digits = 1),
    Mediana_zysku = round.([median(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))), 
    median(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))),
    median(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))),
    median(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)))], digits = 1))
, env = :table) |> print

Plots.savefig(p_h2_1, pwd() * "\\plots\\h2\\ecdf_profit.svg")
Plots.savefig(p_h2_1, pwd() * "\\plots_to_export\\ecdf_profit.pdf")

DataFrame(Typ = ["Producenci nie badają konsumentów", "Producent bada konsumentów jako jedyny", "Producent nie bada konsumentów jako jedyny", "Producenci badają konsumentów"], Zysk = round.([
    median(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))), 
    median(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))), 
    median(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))), 
    median(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)))], digits = 1))

# Tak, dlaczego?

p_h2_q = Plots.plot(ecdf(mean.(getindex.(ex3_v300_quality, 1))), label = "Żaden producent nie bada oczekiwań konsumentów", xlabel = "Średnia jakość", ylabel = "Dystrybuanta empiryczna", legendfontsize = 6, xlim = (0,1))
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_quality, 1))), label = "Tylko producent pierwszy bada oczekiwania konsumentów")
Plots.plot!(ecdf(mean.(getindex.(ex3_v302_quality, 1))), label = "Tylko producent pierwszy nie bada oczekiwań konsumentów")
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_quality, 1))), label = "Obaj producenci badają oczekiwania konsumentów", legend = :bottomright)

p_h2_d = Plots.plot(ecdf(mean.(getindex.(ex3_v300_durability, 1))), label = "Żaden producent nie bada oczekiwań konsumentów", xlabel = "Średnia trwałość", ylabel = "Dystrybuanta empiryczna", legendfontsize = 6)
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_durability, 1))), label = "Tylko producent pierwszy bada oczekiwania konsumentów")
Plots.plot!(ecdf(mean.(getindex.(ex3_v302_durability, 1))), label = "Tylko producent pierwszy nie bada oczekiwań konsumentów")
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_durability, 1))), label = "Obaj producenci badają oczekiwania konsumentów", legend = :bottomright, xlim = (0,1))

p_h2_qd = Plots.plot(p_h2_q, p_h2_d, size = (1200, 600), margin=5Plots.mm)
Plots.savefig(p_h2_qd, pwd() *  "\\plots_to_export\\ecdf_quality_durability_w_wo_research.pdf")

p_h2_m = Plots.plot(ecdf(mean.(getindex.(ex3_v300_margin, 1))), label = "Żaden producent nie bada oczekiwań konsumentów", xlabel = "Średnia marża [%]", ylabel = "Dystrybuanta empiryczna", legendfontsize = 6, xlim = (1, 1.5))
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_margin, 1))), label = "Tylko producent pierwszy bada oczekiwania konsumentów")
Plots.plot!(ecdf(mean.(getindex.(ex3_v302_margin, 1))), label = "Tylko producent pierwszy nie bada oczekiwań konsumentów")
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_margin, 1))), label = "Obaj producenci badają oczekiwania konsumentów", legend = :bottomright)

function calculate_margin_abs(k,d,m,c,h)
    return k .* (1 .- d .^ h) ./ (1 .- d) .* c .* (m .- 1)
end

margin_abs_300 = calculate_margin_abs.(mean.(getindex.(ex3_v300_quality, 1)),
mean.(getindex.(ex3_v300_durability, 1)),
mean.(getindex.(ex3_v300_margin, 1)),
ex3_v300301302303_cc,
ex3_v300301302303_H)

margin_abs_301 = calculate_margin_abs.(mean.(getindex.(ex3_v301_quality, 1)),
mean.(getindex.(ex3_v301_durability, 1)),
mean.(getindex.(ex3_v301_margin, 1)),
ex3_v300301302303_cc,
ex3_v300301302303_H)

margin_abs_302 = calculate_margin_abs.(mean.(getindex.(ex3_v302_quality, 1)),
mean.(getindex.(ex3_v302_durability, 1)),
mean.(getindex.(ex3_v302_margin, 1)),
ex3_v300301302303_cc,
ex3_v300301302303_H)

margin_abs_303 = calculate_margin_abs.(mean.(getindex.(ex3_v303_quality, 1)),
mean.(getindex.(ex3_v303_durability, 1)),
mean.(getindex.(ex3_v303_margin, 1)),
ex3_v300301302303_cc,
ex3_v300301302303_H)

p_h2_ma = Plots.plot(ecdf(margin_abs_300), label = "Żaden producent nie bada oczekiwań konsumentów", xlabel = "Średnia marża", ylabel = "Dystrybuanta empiryczna", legendfontsize = 6, xlim = (0, 0.12))
Plots.plot!(ecdf(margin_abs_301), label = "Tylko producent pierwszy bada oczekiwania konsumentów")
Plots.plot!(ecdf(margin_abs_302), label = "Tylko producent pierwszy nie bada oczekiwań konsumentów")
Plots.plot!(ecdf(margin_abs_303), label = "Obaj producenci badają oczekiwania konsumentów", legend = :bottomright)

p_h2_p = Plots.plot(ecdf(mean.(getindex.(ex3_v300_price, 1))), label = "Żaden producent nie bada oczekiwań konsumentów", xlabel = "Średnia cena", ylabel = "Dystrybuanta empiryczna", legendfontsize = 6, xlim = (0, 1.5))
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_price, 1))), label = "Tylko producent pierwszy bada oczekiwania konsumentów")
Plots.plot!(ecdf(mean.(getindex.(ex3_v302_price, 1))), label = "Tylko producent pierwszy nie bada oczekiwań konsumentów")
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_price, 1))), label = "Obaj producenci badają oczekiwania konsumentów", legend = :bottomright)

Plots.savefig(p_h2_m, pwd() *  "\\plots_to_export\\ecdf_margin_perc_w_wo_research.pdf")
Plots.savefig(p_h2_ma, pwd() *  "\\plots_to_export\\ecdf_margin_abs_w_wo_research.pdf")
Plots.savefig(p_h2_p, pwd() *  "\\plots_to_export\\ecdf_price_w_wo_research.pdf")


p_h2_qs = Plots.plot(ecdf(mean.(getindex.(ex3_v300_quantity_sold, 1))), label = "Żaden producent nie bada oczekiwań konsumentów", xlabel = "Średni popyt", ylabel = "Dystrybuanta empiryczna", legendfontsize = 6, xlim = (0,30))
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_quantity_sold, 1))), label = "Tylko producent pierwszy bada oczekiwania konsumentów")
Plots.plot!(ecdf(mean.(getindex.(ex3_v302_quantity_sold, 1))), label = "Tylko producent pierwszy nie bada oczekiwań konsumentów")
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_quantity_sold, 1))), label = "Obaj producenci badają oczekiwania konsumentów", legend = :bottomright)

p_h2_ss = Plots.plot(ecdf(mean.(getindex.(ex3_v300_quantity_sold,1)) ./ mean.(getindex.(ex3_v300_quantity_produced,1))), xlabel = "Udział sprzedanych dóbr w całości produkcji", ylabel = "Dystrybuanta empiryczna", titlefontsize = 12, label = "Nikt nie bada oczekiwań konsumentów")
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_quantity_sold,1)) ./ mean.(getindex.(ex3_v301_quantity_produced,1))), label = "Producent bada oczekiwania konsumentów jako jedyny")
Plots.plot!(ecdf(mean.(getindex.(ex3_v302_quantity_sold,1)) ./ mean.(getindex.(ex3_v302_quantity_produced,1))), label = "Producent nie bada oczekiwań konsumentów jako jedyny")
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_quantity_sold,1)) ./ mean.(getindex.(ex3_v303_quantity_produced,1))), label = "Obaj producenci badają oczekiwania konsumentów")

Plots.savefig(p_h2_qs, pwd() *  "\\plots_to_export\\ecdf_quantity_sold.pdf")
Plots.savefig(p_h2_ss, pwd() *  "\\plots_to_export\\ecdf_quantity_sold_share.pdf")

##############################################################################################

get_percentiles(x) = [percentile(x, 1),
                    percentile(x, 5),
                    percentile(x, 10),
                    percentile(x, 15),
                    percentile(x, 20),
                    percentile(x, 25),
                    percentile(x, 30)]

get_percentile(X,k) = [round(percentile(x, k), digits = 1) for x in X]

Typ = ["FNB", "FBJ", "FNJ", "OFB"]

X = [sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)),
    sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)),
    sum.(getindex.(ex3_v302_producer_surplus_singleton, 1)),
    sum.(getindex.(ex3_v303_producer_surplus_singleton, 1))]

trans_perc = [transition_percentile(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))),
transition_percentile(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))),
transition_percentile(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))),
transition_percentile(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)))]

trans_perc = [transition_percentile(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))),
transition_percentile(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))),
transition_percentile(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))),
transition_percentile(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)))]

loss_percentiles_df = DataFrame(
    Typ = Typ,
    Trans_perc = trans_perc,
    Percentile1 = get_percentile(X, 1),
    Percentile5 = get_percentile(X, 5),
    Percentile10 = get_percentile(X, 10),
    Percentile20 = get_percentile(X, 20),
    Percentile30= get_percentile(X, 30),
    Percentile40 = get_percentile(X, 40)
    )

latexify(
    loss_percentiles_df
, env = :table) |> print



DataFrame(Typ = ["Producenci nie badają konsumentów", "Producent bada konsumentów jako jedyny", "Producent nie bada konsumentów jako jedyny", "Producenci badają konsumentów"], Zysk = round.([
    mean(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)) .<= 0]), 
    mean(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)) .<= 0]), 
    mean(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v302_producer_surplus_singleton, 1)) .<= 0]), 
    mean(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)) .<= 0])], digits = 1))

UnequalVarianceTTest(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)) .<= 0], 
                    sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)) .<= 0])
UnequalVarianceTTest(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)) .<= 0], 
                    sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v302_producer_surplus_singleton, 1)) .<= 0])
UnequalVarianceTTest(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)) .<= 0], 
                    sum.(getindex.(ex3_v303_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)) .<= 0])        

p_h2_p_h = Plots.scatter(sort(unique(ex3_v300301302303_H)), [transition_percentile(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))[ex3_v300301302303_H .== h]) for h in sort(unique(ex3_v300301302303_H))], xlabel = "Okres przydatności dobra H", ylabel = "Prawdopodobieństwo poniesienia straty, [%]", label = "Producent nie prowadzi badań", markerstrokewidth = 0, markersize = 4)
Plots.scatter!(sort(unique(ex3_v300301302303_H)), [transition_percentile(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))[ex3_v300301302303_H .== h]) for h in sort(unique(ex3_v300301302303_H))], label = "Producent prowadzi badania", markerstrokewidth = 0, markersize = 4)

Plots.savefig(p_h2_p_h, pwd() * "\\plots_to_export\\loss probability vs h.pdf")

####################################### Hipoteza 3 #########################################

"""
Istnieje zależność pomiędzy częstotliwością wymiany informacji pomiędzy konsumentami a poziomem nadzwyczajnego zysku firmy wynikającego z prowadzenia badań konsumenckich
"""

avg_quality_diff_300 = mean.(getindex.(ex3_v300_quality, 1)) .- mean.(getindex.(ex3_v300_quality_exp, 1))
avg_quality_diff_301 = mean.(getindex.(ex3_v301_quality, 1)) .- mean.(getindex.(ex3_v301_quality_exp, 1))
avg_quality_diff_302 = mean.(getindex.(ex3_v302_quality, 1)) .- mean.(getindex.(ex3_v302_quality_exp, 1))
avg_quality_diff_303 = mean.(getindex.(ex3_v303_quality, 1)) .- mean.(getindex.(ex3_v303_quality_exp, 1))

p = Plots.plot(legend = nothing)
ds300 = getindex.(ex3_v300_quality, 1) .- getindex.(ex3_v300_quality_exp, 1)
ds300b = getindex.(ex3_v300_quality, 1) .- getindex.(ex3_v300_quality_exp_buyers, 1)
ds301 = getindex.(ex3_v301_quality, 1) .- getindex.(ex3_v301_quality_exp, 1)
ds301b = getindex.(ex3_v301_quality, 1) .- getindex.(ex3_v301_quality_exp_buyers, 1)

for q in ds300
    p = Plots.plot!(q, color = "grey", alpha = 0.1)
end

id_learning = (ex3_v300301302303_Li .>= 0.9) .& (ex3_v300301302303_Lw .>= 0.9)
sum(id_learing)
id_connections = (ex3_v300301302303_nl .> 800) 

mean_missing(x) = mean(x[.!ismissing.(x)])

mean_missing([1,2,3,missing])

Plots.plot(ds300[1:100], color = "grey", alpha = 0.1, legend = false)

p = Plots.plot!(mean(ds300[(ex3_v300301302303_nl .>= 2000) .& id_learning]), color = "red", linewidth = 2)
p = Plots.plot!(mean(ds301[(ex3_v300301302303_nl .>= 2000) .& id_learning]), color = "green", linewidth = 2)
p = Plots.plot!(mean_missing(ds300b), color = "red", linewidth = 2, linestyle = :dash)
p = Plots.plot!(mean_missing(ds301b), color = "green", linewidth = 2)

Plots.plot(mean(getindex.(ex3_v300_quality, 1)))
Plots.plot!(mean(getindex.(ex3_v300_quality_exp, 1)))

Plots.plot(mean(getindex.(ex3_v301_quality, 1)))
Plots.plot!(mean(getindex.(ex3_v301_quality_exp, 1)))

p_h3_1 = Plots.scatter(sort(unique(ex3_v300301302303_nl)), [mean(avg_quality_diff_300[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], smooth = true, xlabel = "Liczba połączeń w sieci", ylabel = "Różnica między oczekiwaniami konsumentów \n a średnim poziomem jakości", label = "Brak badań konsumentów", titlefontsize = 10)
Plots.scatter!(sort(unique(ex3_v300301302303_nl)), [mean(avg_quality_diff_301[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], smooth = true, label = "Wykonywanie badań konsumentów")
Plots.scatter!(sort(unique(ex3_v300301302303_nl)), [mean(avg_quality_diff_302[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], smooth = true, label = "Wykonywanie badań konsumentów")
Plots.scatter!(sort(unique(ex3_v300301302303_nl)), [mean(avg_quality_diff_303[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], smooth = true, label = "Wykonywanie badań konsumentów")

Plots.savefig(p_h3_1, pwd() * "\\plots\\h3\\expectation bias vs connection number.svg")

# Większa liczba połączeń prowadzi do zbliżenia oczekiwań konsumentów do rzeczywistych wartości

dK300 = cut_integer(abs.(mean.(getindex.(ex3_v300_quality, 1)) .- mean.(getindex.(ex3_v300_quality_exp, 1))), 8)
dK301 = cut_integer(abs.(mean.(getindex.(ex3_v301_quality, 1)) .- mean.(getindex.(ex3_v301_quality_exp, 1))), 8)

p_h3_2 = Plots.scatter(dK300[2], [mean(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))[dK300[1] .== g]) for g in 1:8], smooth = true, xlabel = "Różnica między oczekiwaniami konsumentów \n a średnim poziomem jakości", ylabel = "Średni zysk firmy", label = "Brak badań konsumentów")
Plots.scatter!(dK301[2], [mean(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))[dK301[1] .== g]) for g in 1:8], smooth = true, label = "Wykonywanie badań konsumentów")

Plots.savefig(p_h3_2, pwd() * "\\plots\\h3\\profit vs expectation bias.svg")

Plots.scatter(sort(unique(ex3_v300301302303_nl)), [mean(mean.(getindex.(ex3_v300_producer_surplus_singleton, 1))[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], smooth = true)
Plots.scatter!(sort(unique(ex3_v300301302303_nl)), [mean(mean.(getindex.(ex3_v301_producer_surplus_singleton, 1))[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], smooth = true)
Plots.scatter!(sort(unique(ex3_v300301302303_nl)), [mean(mean.(getindex.(ex3_v302_producer_surplus_singleton, 1))[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], smooth = true)
Plots.scatter!(sort(unique(ex3_v300301302303_nl)), [mean(mean.(getindex.(ex3_v303_producer_surplus_singleton, 1))[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], smooth = true)

π300_hat1 = fitted(GLM.lm(@formula(y~x), DataFrame(x = sort(unique(ex3_v300301302303_nl)), y = [mean(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))])))

π301_hat1 = fitted(GLM.lm(@formula(y~x), DataFrame(x = sort(unique(ex3_v300301302303_nl)), y = [mean(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))])))

π302_hat1 = fitted(GLM.lm(@formula(y~x), DataFrame(x = sort(unique(ex3_v300301302303_nl)), y = [mean(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))])))

π303_hat1 = fitted(GLM.lm(@formula(y~x), DataFrame(x = sort(unique(ex3_v300301302303_nl)), y = [mean(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1))[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))])))

π300_hat2 = fitted(GLM.lm(@formula(y~x), DataFrame(x = sort(unique(ex3_v300301302303_nl)), y = [mean(sum.(getindex.(ex3_v300_producer_surplus_singleton, 2))[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))])))

π301_hat2 = fitted(GLM.lm(@formula(y~x), DataFrame(x = sort(unique(ex3_v300301302303_nl)), y = [mean(sum.(getindex.(ex3_v301_producer_surplus_singleton, 2))[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))])))

π302_hat2 = fitted(GLM.lm(@formula(y~x), DataFrame(x = sort(unique(ex3_v300301302303_nl)), y = [mean(sum.(getindex.(ex3_v302_producer_surplus_singleton, 2))[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))])))

π303_hat2 = fitted(GLM.lm(@formula(y~x), DataFrame(x = sort(unique(ex3_v300301302303_nl)), y = [mean(sum.(getindex.(ex3_v303_producer_surplus_singleton, 2))[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))])))

p_h3_3 = Plots.plot(sort(unique(ex3_v300301302303_nl)), π301_hat1 .- π300_hat1, xlabel = "Liczba połączeń w sieci", ylabel = "Dodatkowy zysk firmy z tytułu \n prowadzenia badań konsumenckich", label = "Wdrożenie badań, gdy konkurent ich nie wykonuje",marker = :circle, markerstrokewidth = 0)
Plots.plot!(sort(unique(ex3_v300301302303_nl)), π303_hat1 .- π302_hat1, label = "Wdrożenie badań, gdy konkurent też je wykonuje", widen = true, marker = :square, markerstrokewidth = 0)

Plots.savefig(p_h3_3, pwd() * "\\plots\\h3\\incremental profit vs connection number.svg")
Plots.savefig(p_h3_3, pwd() * "\\plots_to_export\\incremental profit vs connection number.pdf")

# Średnia liczba sygnałów a liczba połączeń

Plots.scatter(sort(unique(ex3_v300301302303_nl)), [mean(mean.(ex3_v300_sv)[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))])

incr_profit = mean.(getindex.(ex3_v301_producer_surplus_singleton, 1)) .- mean.(getindex.(ex3_v300_producer_surplus_singleton, 1))

incr_profit_l = [mean(incr_profit[(ex3_v300301302303_Li .== li)]) for li in sort(unique(ex3_v300301302303_Li))]

incr_profit_t = [mean(incr_profit[(ex3_v300301302303_Lw .== lw)]) for lw in sort(unique(ex3_v300301302303_Lw))]

model_l = loess(ex3_v300301302303_Li, incr_profit, span = 0.75)
model_t = loess(ex3_v300301302303_Lw, incr_profit, span = 0.75)

us = range(extrema(0:1.0)...; step = 0.01)
vs_l = Loess.predict(model_l, us)
vs_t = Loess.predict(model_t, us)

lt_df = DataFrame(x = 1 ./ (sort(unique(ex3_v300301302303_Li)) .+ 0.01), yl = incr_profit_l, yt = incr_profit_t)

lt_df_pred = DataFrame(x = 1 ./ (0.01:0.01:1.0))

model_l = GLM.lm(@formula(yl~x), lt_df)
model_t = GLM.lm(@formula(yt~x), lt_df)

fit_l = GLM.predict(model_l, lt_df_pred)
fit_t = GLM.predict(model_t, lt_df_pred)

Plots.scatter(sort(unique(ex3_v300301302303_Li)), incr_profit_l, markerstrokewidth = 0, label = "Dodatkowy zysk, sygnały indywidualne, " * L"\lambda", color = :blue)
Plots.plot!(0.00:0.01:0.99, fit_l, color = "blue", label = "")
Plots.scatter!(sort(unique(ex3_v300301302303_Lw)), incr_profit_t, legend = :topleft, markerstrokewidth = 0, label = "Dodatkowy zysk, sygnały od innych konsumentów, " * L"\theta", color = :green)
Plots.plot!(0.00:0.01:0.99, fit_t, color = "green", legend = :topright, xlabel = "Wagi sygnałów " * L"\lambda, \theta", ylabel = "Dodatkowy zysk firmy z tytułu \n prowadzenia badań konsumenckich", label = "")

λ = ex3_v300301302303_Li
θ = ex3_v300301302303_Lw

nl_id = ex3_v300301302303_nl .<= 600

hm_rs = [mean(incr_profit[(λ .== l) .& (θ .== t) .& nl_id]) for l in sort(unique(λ)), t in sort(unique(θ))]

Plots.plot([mean(x) for x in eachrow(hm_rs)]) # by λ
Plots.plot([mean(x) for x in eachcol(hm_rs)]) # by θ

ex4_p1 = StatsPlots.heatmap(sort(unique(λ))[2:end], sort(unique(θ))[2:end], hm_rs'[2:end, 2:end], xlabel = "λ, produkty oceniane osobiście", ylabel = "θ, produkty oceniane przez sąsiadów", titlefontsize = 8)

########## Hipoteza 4 ############

"""
Decyzja producenta o uwzględnianiu wyników badań preferencji konsumentów w prowadzonych przez niego działaniach zależy od ceny pomiaru preferencji oraz horyzontu czasowego planowania producenta
"""

low_nl = ex3_v300301302303_nl .<= 800

payoff1_low_nl = [mean(
        [median(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))[low_nl]), median(sum.(getindex.(ex3_v300_producer_surplus_singleton, 2))[low_nl])
]),
        mean(
            [median(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))[low_nl]), median(sum.(getindex.(ex3_v302_producer_surplus_singleton, 2))[low_nl])
]),
        mean(
            [median(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))[low_nl]), median(sum.(getindex.(ex3_v301_producer_surplus_singleton, 2))[low_nl])
]),
        mean(
            [median(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1))[low_nl]), median(sum.(getindex.(ex3_v303_producer_surplus_singleton, 2))[low_nl])
])]

payoff1_high_nl = [mean(
        [median(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))[.!low_nl]), median(sum.(getindex.(ex3_v300_producer_surplus_singleton, 2))[.!low_nl])
]),
        mean(
            [median(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))[.!low_nl]), median(sum.(getindex.(ex3_v302_producer_surplus_singleton, 2))[.!low_nl])
]),
        mean(
            [median(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))[.!low_nl]), median(sum.(getindex.(ex3_v301_producer_surplus_singleton, 2))[.!low_nl])
]),
        mean(
            [median(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1))[.!low_nl]), median(sum.(getindex.(ex3_v303_producer_surplus_singleton, 2))[.!low_nl])
])]

payoff2_low_nl = payoff1_low_nl[[1,3,2,4]]
payoff2_high_nl = payoff1_high_nl[[1,3,2,4]]

costs = collect(0:50)

#### High NL

nash_equilibriums_int_high_nl = simulate_ne_for_costs(payoff1_high_nl, payoff2_high_nl, costs, "pNE")
nash_equilibriums_PD_high_nl = simulate_ne_for_costs(payoff1_high_nl, payoff2_high_nl, costs, "PD")
mixedNE_high_nl = simulate_ne_for_costs(payoff1_high_nl, payoff2_high_nl, costs, "mNE")
valid_mixed_NE_high_nl = [all((x .>= 0) .& (x .<= 1)) for x in mixedNE_high_nl]

plot_NE = Plots.plot(costs, true_or_missing.(4 .∈ nash_equilibriums_int_high_nl, k=2), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów",  markercolor = "black", linecolor = "black", label = "Strategia równowagi, strategie czyste")
#Plots.plot!(costs, true_or_missing.(3 .∈ nash_equilibriums_int; k = 3), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
#Plots.plot!(costs, true_or_missing.(2 .∈ nash_equilibriums_int; k = 2), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
Plots.plot!(costs, true_or_missing.(1 .∈ nash_equilibriums_int_high_nl; k = 1), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
Plots.plot!(costs[nash_equilibriums_PD_high_nl .== 1], sum.(nash_equilibriums_int_high_nl)[nash_equilibriums_PD_high_nl .== 1], markercolor = "red", markershape = :square, markerstrokewidth = 0, markersize = 6, linecolor = "red", label = "Rozwiązanie będące dylematem więźnia")
Plots.plot!(costs[valid_mixed_NE_high_nl], fill(1, count(valid_mixed_NE_high_nl)), markercolor = "blue", markershape = :square, markersize = 6, markerstrokewidth = 0, label = "Dodatkowe rozwiązanie w strategiach mieszanych")
Plots.plot!(costs[valid_mixed_NE_high_nl], fill(2, count(valid_mixed_NE_high_nl)), markercolor = "blue", markershape = :square, markersize = 6, markerstrokewidth = 0, label = "", legend = :topright)
Plots.plot!(yticks = ([1,2], ["Żaden producent nie prowadzi badań", "Obaj producenci prowadzą badania"]), size = (1000, 200))
xaxis!(xlabel = "Koszt pomiaru preferencji")

p_high_nl = Plots.plot(payoff1_high_nl[[1,2,4,3,1]], payoff2_high_nl[[1,2,4,3,1]], legend = :outerright, marker = :circle, label = "C = 0", xlabel = "Wypłata firmy #1", ylabel = "Wypłata firmy #2", xlim = (0,80), ylim = (0,80))
annotate!.(payoff1_high_nl[[1,2,4,3]] .+ [-2,2,2,-2], payoff2_high_nl[[1,2,4,3]] .+ [-2,2,2,2], Plots.text.(["(N,N)", "(R,N)", "(R,R)", "(N,R)"], fill(8,4), :blue), alpha = 0.5)

c = 15; Plots.plot!(payoff1_high_nl[[1,2,4,3,1]] .- [0, c, c, 0, 0], payoff2_high_nl[[1,2,4,3,1]] .- [0, 0, c, c, 0], marker = :circle, label = "C = " * string(c))
annotate!.(payoff1_high_nl[[2,4,3]] .+ [2,2,-2] .- [c, c, 0], payoff2_high_nl[[2,4,3]] .+ [2,2,-2] .- [0, c, c], Plots.text.(["(R,N)", "(R,R)", "(N,R)"], fill(8,3), :red), linealpha = 0.5)

c = 30; Plots.plot!(payoff1_high_nl[[1,2,4,3,1]] .- [0, c, c, 0, 0], payoff2_high_nl[[1,2,4,3,1]] .- [0, 0, c, c, 0], marker = :circle, label = "C = " * string(c))
annotate!.(payoff1_high_nl[[2,4,3]] .+ [2,2,-2] .- [c, c, 0], payoff2_high_nl[[2,4,3]] .+ [2,2,-2] .- [0, c, c], Plots.text.(["(R,N)", "(R,R)", "(N,R)"], fill(8,3), :green), linealpha = 0.5)

c = 45; Plots.plot!(payoff1_high_nl[[1,2,4,3,1]] .- [0, c, c, 0, 0], payoff2_high_nl[[1,2,4,3,1]] .- [0, 0, c, c, 0], marker = :circle, label = "C = " * string(c))
annotate!.(payoff1_high_nl[[2,4,3]] .+ [2,2,-2] .- [c, c, 0], payoff2_high_nl[[2,4,3]] .+ [2,2,-2] .- [0, c, c], Plots.text.(["(R,N)", "(R,R)", "(N,R)"], fill(8,3), :purple), linealpha = 0.5)

Plots.savefig(p_high_nl, pwd() * "\\plots_to_export\\nash_equilibriums_h.pdf")
Plots.savefig(p, pwd() * "\\plots\\h4\\nash_equilibriums_h.svg")


#### Low NL

nash_equilibriums_int_low_nl = simulate_ne_for_costs(payoff1_low_nl, payoff2_low_nl, costs, "pNE")
nash_equilibriums_PD_low_nl = simulate_ne_for_costs(payoff1_low_nl, payoff2_low_nl, costs, "PD")
mixedNE_low_nl = simulate_ne_for_costs(payoff1_low_nl, payoff2_low_nl, costs, "mNE")
valid_mixed_NE_low_nl = [all((x .>= 0) .& (x .<= 1)) for x in mixedNE_low_nl]

plot_NE = Plots.plot(costs, true_or_missing.(4 .∈ nash_equilibriums_int_low_nl, k=2), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów",  markercolor = "black", linecolor = "black", label = "Strategia równowagi, strategie czyste")
#Plots.plot!(costs, true_or_missing.(3 .∈ nash_equilibriums_int; k = 3), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
#Plots.plot!(costs, true_or_missing.(2 .∈ nash_equilibriums_int; k = 2), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
Plots.plot!(costs, true_or_missing.(1 .∈ nash_equilibriums_int_low_nl; k = 1), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
Plots.plot!(costs[nash_equilibriums_PD_low_nl .== 1], sum.(nash_equilibriums_int_low_nl)[nash_equilibriums_PD_low_nl .== 1], markercolor = "red", markershape = :square, markerstrokewidth = 0, markersize = 6, linecolor = "red", label = "Rozwiązanie będące dylematem więźnia")
Plots.plot!(costs[valid_mixed_NE_low_nl], fill(1, count(valid_mixed_NE_low_nl)), markercolor = "blue", markershape = :square, markersize = 6, markerstrokewidth = 0, label = "Dodatkowe rozwiązanie w strategiach mieszanych")
Plots.plot!(costs[valid_mixed_NE_low_nl], fill(2, count(valid_mixed_NE_low_nl)), markercolor = "blue", markershape = :square, markersize = 6, markerstrokewidth = 0, label = "", legend = :topright)
Plots.plot!(yticks = ([1,2], ["Żaden producent nie prowadzi badań", "Obaj producenci prowadzą badania"]), size = (1000, 200))
xaxis!(xlabel = "Koszt pomiaru preferencji")

p_low_nl = Plots.plot(payoff1_low_nl[[1,2,4,3,1]], payoff2_low_nl[[1,2,4,3,1]], legend = :outerright, marker = :circle, label = "C = 0", xlabel = "Wypłata firmy #1", ylabel = "Wypłata firmy #2", xlim = (0,80), ylim = (0,80))
annotate!.(payoff1_low_nl[[1,2,4,3]] .+ [-2,2,2,-2], payoff2_low_nl[[1,2,4,3]] .+ [-2,2,2,2], Plots.text.(["(N,N)", "(R,N)", "(R,R)", "(N,R)"], fill(8,4), :blue), alpha = 0.5)

c = 15; Plots.plot!(payoff1_low_nl[[1,2,4,3,1]] .- [0, c, c, 0, 0], payoff2_low_nl[[1,2,4,3,1]] .- [0, 0, c, c, 0], marker = :circle, label = "C = " * string(c))
annotate!.(payoff1_low_nl[[2,4,3]] .+ [2,2,-2] .- [c, c, 0], payoff2_low_nl[[2,4,3]] .+ [2,2,-2] .- [0, c, c], Plots.text.(["(R,N)", "(R,R)", "(N,R)"], fill(8,3), :red), linealpha = 0.5)

c = 30; Plots.plot!(payoff1_low_nl[[1,2,4,3,1]] .- [0, c, c, 0, 0], payoff2_low_nl[[1,2,4,3,1]] .- [0, 0, c, c, 0], marker = :circle, label = "C = " * string(c))
annotate!.(payoff1_low_nl[[2,4,3]] .+ [2,2,-2] .- [c, c, 0], payoff2_low_nl[[2,4,3]] .+ [2,2,-2] .- [0, c, c], Plots.text.(["(R,N)", "(R,R)", "(N,R)"], fill(8,3), :green), linealpha = 0.5)

c = 45; Plots.plot!(payoff1_low_nl[[1,2,4,3,1]] .- [0, c, c, 0, 0], payoff2_low_nl[[1,2,4,3,1]] .- [0, 0, c, c, 0], marker = :circle, label = "C = " * string(c))
annotate!.(payoff1_low_nl[[2,4,3]] .+ [2,2,-2] .- [c, c, 0], payoff2_low_nl[[2,4,3]] .+ [2,2,-2] .- [0, c, c], Plots.text.(["(R,N)", "(R,R)", "(N,R)"], fill(8,3), :purple), linealpha = 0.5)

p_low_nl

Plots.savefig(p_low_nl, pwd() * "\\plots_to_export\\nash_equilibriums.pdf")
Plots.savefig(p_low_nl, pwd() * "\\plots\\h4\\nash_equilibriums.svg")



######## EXPORT PLOTS ##############

source_folder = "C:\\Users\\User\\Documents\\GitHub\\PhD-Model\\plots_to_export\\"
target_folder = "C:\\Users\\User\\Documents\\GitHub\\PhD---rozprawa\\plots\\"

for file in readdir(source_folder)
    cp(
        source_folder * file,
        target_folder * file,
        force = true
    )
end