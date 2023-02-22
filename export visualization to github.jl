include(pwd() * "\\methods\\methods.jl")



@load "C:\\Users\\User\\Documents\\PhDWorkspace_4cases_2.jld2"


############################################################################################
######################################## HISTORIA ##########################################

####################################### Hipoteza 2 #########################################

"""
Ryzyko poniesienia przez producenta straty jest uzależnione od uwzględnienia preferencji konsumentów podczas podejmowania decyzji oraz od długości okresu przydatności produktów oferowanych na rynku
"""

# Czy zyski producentów prowadzących badania i nie się różnią?

p_h2_1 = plot_ecdf(true, sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)), "Nikt nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "Dystrybuanta empiryczna", title = "Dystrybuanta empiryczna - zysk producenta", xlim = (-300,1200))
plot_ecdf(false, sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)), "Producent bada oczekiwania konsumentów jako jedyny")
#plot_ecdf(false, sum.(getindex.(ex3_v302_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów jako jedyny")
plot_ecdf(false, sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)), "Obaj producenci badają oczekiwania konsumentów")

Plots.savefig(p_h2_1, pwd() * "\\plots\\h2\\ecdf profit.svg")
Plots.savefig(p_h2_1, pwd() * "\\plots_to_export\\ecdf profit.pdf")

DataFrame(Typ = ["Producenci nie badają konsumentów", "Producent bada konsumentów jako jedyny", "Producent nie bada konsumentów jako jedyny", "Producenci badają konsumentów"], Zysk = round.([
    median(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))), 
    median(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))), 
    median(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))), 
    median(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)))], digits = 1))

# Tak, dlaczego?

Plots.plot(ecdf(mean.(getindex.(ex3_v300_quantity_sold,1))), xlim = (0,50))
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_quantity_sold,1))))

Plots.plot(ecdf(mean.(getindex.(ex3_v300_quantity_produced,1))), xlim = (0,50))
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_quantity_produced,1))))

p_h2_2 = Plots.plot(ecdf(mean.(getindex.(ex3_v300_quantity_sold,1)) ./ mean.(getindex.(ex3_v300_quantity_produced,1))), xlabel = "Udział sprzedanych dóbr w całości produkcji", ylabel = "Dystrybuanta empiryczna", title = "Dystrybuanta empiryczna - % sprzedanych produktów", titlefontsize = 12, label = "Nikt nie bada oczekiwań konsumentów")
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_quantity_sold,1)) ./ mean.(getindex.(ex3_v301_quantity_produced,1))), label = "Producent bada oczekiwania konsumentów jako jedyny")
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_quantity_sold,1)) ./ mean.(getindex.(ex3_v303_quantity_produced,1))), label = "Obaj producenci badają oczekiwania konsumentów")

Plots.savefig(p_h2_2, pwd() * "\\plots\\h2\\ecdf rotation.svg")
Plots.savefig(p_h2_2, pwd() * "\\plots_to_export\\ecdf rotation.pdf")

Plots.plot(ecdf(mean.(getindex.(ex3_v300_quantity_sold,1)) ./ mean.(getindex.(ex3_v300_quantity_produced,1))))
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_quantity_sold,1)) ./ mean.(getindex.(ex3_v301_quantity_produced,1))))

mean(mean.(getindex.(ex3_v300_quantity_sold,1))) / mean(mean.(getindex.(ex3_v300_quantity_produced,1)))
mean(mean.(getindex.(ex3_v301_quantity_sold,1)))/mean(mean.(getindex.(ex3_v301_quantity_produced,1)))

# PB ma niższe ceny, stąd wyższą wyprzedawalność

Plots.plot(ecdf(mean.(getindex.(ex3_v300_price,1))))
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_price,1))))

ATC200 = mean.(fATC.(getindex.(ex3_v300_quality,1), getindex.(ex3_v300_durability,1), ex3_v300301302303_H, ex3_v300301302303_cc))
ATC301 = mean.(fATC.(getindex.(ex3_v301_quality,1), getindex.(ex3_v301_durability,1), ex3_v300301302303_H, ex3_v300301302303_cc))

ex3_4_pl = plot_ecdf(true, ATC300, "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, ATC301, "Producent bada oczekiwania konsumentów")

ex3_4_pl = plot_ecdf(true, mean.(getindex.(ex3_v300_quality,1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, mean.(getindex.(ex3_v301_quality,1)), "Producent bada oczekiwania konsumentów")

ex3_4_pl = plot_ecdf(true, mean.(getindex.(ex3_v300_durability,1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, mean.(getindex.(ex3_v301_durability,1)), "Producent bada oczekiwania konsumentów")

# Czy ma niższe przychody czy koszty?

selling_income_300 = sum.(multi.(getindex.(ex3_v300_price, 1), getindex.(ex3_v300_quantity_sold, 1)))
utilization_cost_300 = (1 .- ex3_v300301302303_cu) .* sum.(multi.(divide.(getindex.(ex3_v300_price, 1), getindex.(ex3_v300_margin, 1)), getindex.(ex3_v300_quantity_produced, 1) .- getindex.(ex3_v300_quantity_sold, 1)))
        
TR300 = selling_income_300 .+ utilization_cost_300
        
selling_income_301 = sum.(multi.(getindex.(ex3_v301_price, 1), getindex.(ex3_v301_quantity_sold, 1)))
utilization_cost_301 = (1 .- ex3_v300301302303_cu) .* sum.(multi.(divide.(getindex.(ex3_v301_price, 1), getindex.(ex3_v301_margin, 1)), getindex.(ex3_v301_quantity_produced, 1) .- getindex.(ex3_v301_quantity_sold, 1)))
        
TR301 = selling_income_301 .+ utilization_cost_301

Plots.plot(ecdf(TR300), xlim = (0,10200))
Plots.plot!(ecdf(TR301))

# Całkowity przychód jest zazwyczaj mniejszy

Plots.plot(ecdf(selling_income_300), xlim = (0,10000))
Plots.plot!(ecdf(selling_income_301))

Plots.plot(ecdf(utilization_cost_300), xlim = (0,1000))
Plots.plot!(ecdf(utilization_cost_301))

[mean(selling_income_300), mean(utilization_cost_300)] ./ mean(TR300)
[mean(selling_income_301), mean(utilization_cost_301)] ./ mean(TR301)

mean(TR301) - mean(TR300)
[mean(selling_income_301), mean(utilization_cost_301)] .- [mean(selling_income_300), mean(utilization_cost_300)]

CS300 = sum.(fTC.(getindex.(ex3_v300_quantity_produced, 1), getindex.(ex3_v300_quality, 1), getindex.(ex3_v300_durability, 1), ex3_v300301302303_H, ex3_v300301302303_cc))
CS301 = sum.(fTC.(getindex.(ex3_v301_quantity_produced, 1), getindex.(ex3_v301_quality, 1), getindex.(ex3_v301_durability, 1), ex3_v300301302303_H, ex3_v300301302303_cc))

Plots.plot(ecdf(CS300), xlim = (0,10200))
Plots.plot!(ecdf(CS301))

##############################################################################################

get_percentiles(x) = [percentile(x, 1),
                    percentile(x, 5),
                    percentile(x, 10),
                    percentile(x, 15),
                    percentile(x, 20),
                    percentile(x, 25),
                    percentile(x, 30)]

get_percentiles(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)))
get_percentiles(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)))
get_percentiles(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1)))
get_percentiles(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)))

transition_percentile(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)))
transition_percentile(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)))
transition_percentile(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1)))
transition_percentile(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)))

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

p_h2_3 = Plots.scatter(sort(unique(ex3_v300301302303_H)), [transition_percentile(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))[ex3_v300301302303_H .== h]) for h in sort(unique(ex3_v300301302303_H))], xlabel = "Okres przydatności dobra H", ylabel = "Prawdopodobieństwo poniesienia straty, [%]", label = "Producent nie prowadzi badań", markerstrokewidth = 0, markersize = 6)
Plots.scatter!(sort(unique(ex3_v300301302303_H)), [transition_percentile(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))[ex3_v300301302303_H .== h]) for h in sort(unique(ex3_v300301302303_H))], label = "Producent prowadzi badania", markerstrokewidth = 0, markersize = 6)

Plots.savefig(p_h2_3, pwd() * "\\plots\\h2\\loss probability vs h.svg")

####################################### Hipoteza 3 #########################################

"""
Istnieje zależność pomiędzy częstotliwością wymiany informacji pomiędzy konsumentami a poziomem nadzwyczajnego zysku firmy wynikającego z prowadzenia badań konsumenckich
"""

avg_quality_diff_300 = mean.(getindex.(ex3_v300_quality, 1)) .- mean.(getindex.(ex3_v300_quality_exp, 1))
avg_quality_diff_301 = mean.(getindex.(ex3_v301_quality, 1)) .- mean.(getindex.(ex3_v301_quality_exp, 1))
avg_durability_diff_300 = mean.(getindex.(ex3_v300_durability, 1)) .- mean.(getindex.(ex3_v300_durability_exp, 1))
avg_durability_diff_301 = mean.(getindex.(ex3_v301_durability, 1)) .- mean.(getindex.(ex3_v301_durability_exp, 1))

p_h3_1 = Plots.scatter(sort(unique(ex3_v300301302303_nl)), [mean(avg_quality_diff_300[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], smooth = true, xlabel = "Liczba połączeń w sieci", ylabel = "Różnica między oczekiwaniami konsumentów \n a średnim poziomem jakości", label = "Brak badań konsumentów", title = "Różnice oczekiwań konsumentów i średnich parametrów produktów \n a liczba połączeń w sieci", titlefontsize = 10)
Plots.scatter!(sort(unique(ex3_v300301302303_nl)), [mean(avg_quality_diff_301[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], smooth = true, label = "Wykonywanie badań konsumentów")

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

Plots.scatter(sort(unique(ex3_v300301302303_Li)), [mean(mean.(getindex.(ex3_v300_producer_surplus_singleton, 1))[(ex3_v300301302303_Li .== li) .& id300]) for li in sort(unique(ex3_v300301302303_Li))], smooth = true, markerstrokewidth = 0)
Plots.scatter!(sort(unique(ex3_v300301302303_Li)), [mean(mean.(getindex.(ex3_v301_producer_surplus_singleton, 1))[(ex3_v300301302303_Li .== li) .& id301]) for li in sort(unique(ex3_v300301302303_Li))], smooth = true, markerstrokewidth = 0)

Plots.scatter(sort(unique(ex3_v300301302303_Lw)), [mean(mean.(getindex.(ex3_v300_producer_surplus_singleton, 1))[(ex3_v300301302303_Lw .== lw) .& id300]) for lw in sort(unique(ex3_v300301302303_Lw))], smooth = true, legend = :topleft, markerstrokewidth = 0)
Plots.scatter!(sort(unique(ex3_v300301302303_Lw)), [mean(mean.(getindex.(ex3_v301_producer_surplus_singleton, 1))[(ex3_v300301302303_Lw .== lw) .& id301]) for lw in sort(unique(ex3_v300301302303_Lw))], smooth = true, legend = :top, markerstrokewidth = 0)

########## Hipoteza 4 ############

"""
Decyzja producenta o uwzględnianiu wyników badań preferencji konsumentów w prowadzonych przez niego działaniach zależy od ceny pomiaru preferencji oraz horyzontu czasowego planowania producenta
"""

payoff1 = round.(median.([sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)), 
sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)), 
sum.(getindex.(ex3_v302_producer_surplus_singleton, 1)), 
sum.(getindex.(ex3_v303_producer_surplus_singleton, 1))]))

payoff2 = round.(median.([sum.(getindex.(ex3_v300_producer_surplus_singleton, 2)), 
sum.(getindex.(ex3_v301_producer_surplus_singleton, 2)), 
sum.(getindex.(ex3_v302_producer_surplus_singleton, 2)), 
sum.(getindex.(ex3_v303_producer_surplus_singleton, 2))]))

payoff1 = round.(median.([vcat(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)), sum.(getindex.(ex3_v300_producer_surplus_singleton, 2))), 
vcat(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)), sum.(getindex.(ex3_v301_producer_surplus_singleton, 2))), 
vcat(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1)), sum.(getindex.(ex3_v302_producer_surplus_singleton, 2))), 
vcat(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)), sum.(getindex.(ex3_v303_producer_surplus_singleton, 2)))]))

payoff2 = payoff1

costs = collect(0:80)

nash_equilibriums_int = simulate_ne_for_costs(payoff1, payoff2, costs, "pNE")
nash_equilibriums_PD = simulate_ne_for_costs(payoff1, payoff2, costs, "PD")
mixedNE = simulate_ne_for_costs(payoff1, payoff2, costs, "mNE")
valid_mixed_NE = [all((x .>= 0) .& (x .<= 1)) for x in mixedNE]

c = 30
construct_payoff_matrix(payoff1 .- [0, c, 0, c], payoff2 .- [0, 0, c, c])

plot_NE = Plots.plot(costs, true_or_missing.(4 .∈ nash_equilibriums_int, k=2), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów",  markercolor = "black", linecolor = "black", label = "Strategia równowagi, strategie czyste")
#Plots.plot!(costs, true_or_missing.(3 .∈ nash_equilibriums_int; k = 3), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
#Plots.plot!(costs, true_or_missing.(2 .∈ nash_equilibriums_int; k = 2), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
Plots.plot!(costs, true_or_missing.(1 .∈ nash_equilibriums_int; k = 1), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
Plots.plot!(costs[nash_equilibriums_PD .== 1], sum.(nash_equilibriums_int)[nash_equilibriums_PD .== 1], markercolor = "red", markershape = :square, markerstrokewidth = 0, markersize = 6, linecolor = "red", label = "Rozwiązanie będące dylematem więźnia")
Plots.plot!(costs[valid_mixed_NE], fill(1, count(valid_mixed_NE)), markercolor = "blue", markershape = :square, markersize = 6, markerstrokewidth = 0, label = "Dodatkowe rozwiązanie w strategiach mieszanych")
Plots.plot!(costs[valid_mixed_NE], fill(2, count(valid_mixed_NE)), markercolor = "blue", markershape = :square, markersize = 6, markerstrokewidth = 0, label = "", legend = :topright)
Plots.plot!(yticks = ([1,2,3,4], ["Żaden producent nie prowadzi badań", "Obaj producenci prowadzą badania"]), size = (1000, 200))
xaxis!(xlabel = "Koszt pomiaru preferencji")

Plots.savefig(plot_NE, pwd() * "\\plots\\h4\\Nash equilibrium cost of measurement.svg")

ex3_4_pl = plot_ecdf(true, sum.(ex3_v300_total_surplus), "Nikt nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(ex3_v301_total_surplus), "Producent bada oczekiwania konsumentów jako jedyny")
#plot_ecdf(false, sum.(getindex.(ex3_v302_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów jako jedyny")
plot_ecdf(false, sum.(ex3_v303_total_surplus), "Obaj producenci badają oczekiwania konsumentów")

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