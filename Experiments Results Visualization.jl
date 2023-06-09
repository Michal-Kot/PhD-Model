include(pwd() * "\\methods\\methods.jl")

using Latexify

@load "C:\\Users\\User\\Documents\\PhDWorkspace_more_nl.jld2"

iterSum = lastindex(ex3_v303_producer_surplus_singleton)

p_h0_1 = Plots.plot(cumsum(mean.(getindex.(ex3_v303_quality, 1))) ./ collect(1:iterSum), xlabel = "Liczba symulacji", ylabel = "Średnia wartość miary", legend = false,  xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6)
p_h0_2 = Plots.plot(cumsum(sum.(getindex.(ex3_v303_durability, 1))) ./ collect(1:iterSum), xlabel = "Liczba symulacji", ylabel = "Średnia wartość miary", legend = false,  xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6)
p_h0_3 = Plots.plot(cumsum(sum.(getindex.(ex3_v303_quantity_sold, 1))) ./ collect(1:iterSum), xlabel = "Liczba symulacji", ylabel = "Średnia wartość miary", legend = false,  xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6)
p_h0_4 = Plots.plot(cumsum(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1))) ./ collect(1:iterSum), xlabel = "Liczba symulacji", ylabel = "Średnia wartość miary", legend = false,  xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6)

p_h0_1234 = Plots.plot(p_h0_1, p_h0_2, p_h0_3, p_h0_4, layout = (2,2))

Plots.savefig(p_h0_1234, pwd() * "\\plots_to_export\\results_stability_means.pdf")

Label1 = "Żaden producent nie bada przekonań konsumentów"
Label2 = "Tylko producent pierwszy bada przekonania konsumentów"
Label3 = "Tylko producent drugi bada przekonania konsumentów"
Label4 = "Obaj producenci badają przekonania konsumentów"

############################################################################################
######################################## HISTORIA ##########################################

####################################### Hipoteza 2 #########################################

"""
Ryzyko poniesienia przez producenta straty jest uzależnione od uwzględnienia preferencji konsumentów podczas podejmowania decyzji oraz od długości okresu przydatności produktów oferowanych na rynku
"""

# Czy zyski producentów prowadzących badania i nie się różnią?

p_h2_1 = plot_ecdf(true, sum.(getindex.(ex3_v300_producer_surplus_singleton, 1)), Label1, xlabel = "Zysk producenta", ylabel = "Dystrybuanta empiryczna", xlim = (-300,1250), ylim = (0,1))
plot_ecdf(false, sum.(getindex.(ex3_v301_producer_surplus_singleton, 1)), Label2)
plot_ecdf(false, sum.(getindex.(ex3_v302_producer_surplus_singleton, 1)), Label3)
plot_ecdf(false, sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)), Label4)

"""p_h2_1 = plot_ecdf(true, sum.(ex3_v300_consumer_surplus), Label1, xlabel = "Zysk producenta pierwszego", ylabel = "Dystrybuanta empiryczna")
plot_ecdf(false, sum.(ex3_v301_consumer_surplus), Label2)
plot_ecdf(false, sum.(ex3_v302_consumer_surplus), Label3)
plot_ecdf(false, sum.(ex3_v303_consumer_surplus), Label4)

p_h2_1 = plot_ecdf(true, sum.(getindex.(ex3_v300_quantity_sold, 1)), Label1, xlabel = "Zysk producenta pierwszego", ylabel = "Dystrybuanta empiryczna")
plot_ecdf(false, sum.(getindex.(ex3_v301_quantity_sold, 1)), Label2)
plot_ecdf(false, sum.(getindex.(ex3_v302_quantity_sold, 1)), Label3)
plot_ecdf(false, sum.(getindex.(ex3_v303_quantity_sold, 1)), Label4)"""

Plots.savefig(p_h2_1, pwd() * "\\plots_to_export\\ecdf_profit.pdf")

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
    DataFrame(Typ = [LLabel1, LLabel2, L"Producent drugi bada oczekiwania konsumentów jako jedyny", LLabel4]),
    Mediana_zysku = round.([median(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))), 
    median(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))),
    median(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))),
    median(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)))], digits = 1)
, env = :table) |> print




DataFrame(Typ = ["Producenci nie badają konsumentów", "Producent bada konsumentów jako jedyny", "Producent nie bada konsumentów jako jedyny", "Producenci badają konsumentów"], Zysk = round.([
    median(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))), 
    median(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))), 
    median(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))), 
    median(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1)))], digits = 1))

# Tak, dlaczego?

p_h2_q01 = Plots.plot(ecdf(mean.(getindex.(ex3_v300_quality, 1))), label = "Żaden producent nie bada \n przekonań konsumentów", xlabel = "Średnia jakość", ylabel = "Dystrybuanta empiryczna", legendfontsize = 12, xlim = (0,1), xlabelfontsize = 14, ylabelfontsize = 14)
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_quality, 1))), label = "Tylko producent pierwszy bada \n przekonania konsumentów", legend = :bottomright)

p_h2_q23 = Plots.plot(ecdf(mean.(getindex.(ex3_v302_quality, 1))), label = "Tylko producent drugi bada \n przekonania konsumentów", xlabel = "Średnia jakość", ylabel = "Dystrybuanta empiryczna", legendfontsize = 12, xlim = (0,1), color = get_color_palette(:auto, 1)[3], xlabelfontsize = 14, ylabelfontsize = 14)
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_quality, 1))), label = "Obaj producenci badają \n przekonania konsumentów", legend = :bottomright, color = get_color_palette(:auto, 1)[4])

p_h2_d01 = Plots.plot(ecdf(mean.(getindex.(ex3_v300_durability, 1))), label = "Żaden producent nie bada \n przekonań konsumentów", xlabel = "Średnia trwałość", ylabel = "Dystrybuanta empiryczna", legendfontsize = 12, xlim = (0,1), xlabelfontsize = 14, ylabelfontsize = 14)
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_durability, 1))), label = "Tylko producent pierwszy bada \n przekonania konsumentów", legend = :bottomright)

p_h2_d23 = Plots.plot(ecdf(mean.(getindex.(ex3_v302_durability, 1))), label = "Tylko producent drugi bada \n przekonania konsumentów", legendfontsize = 12, xlim = (0,1), xlabel = "Średnia trwałość", ylabel = "Dystrybuanta empiryczna", color = get_color_palette(:auto, 1)[3], xlabelfontsize = 14, ylabelfontsize = 14)
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_durability, 1))), label = "Obaj producenci badają \n przekonania konsumentów", legend = :bottomright, xlim = (0,1), color = get_color_palette(:auto, 1)[4])

p_h2_qd = Plots.plot(p_h2_q01, p_h2_q23, p_h2_d01, p_h2_d23, margin=5Plots.mm, layout = (2,2), size = (1600, 1200))
Plots.savefig(p_h2_qd, pwd() *  "\\plots_to_export\\ecdf_quality_durability_w_wo_research.pdf")

p_h2_m = Plots.plot(ecdf(mean.(getindex.(ex3_v300_margin, 1))), label = Label1, xlabel = "Średnia marża [%]", ylabel = "Dystrybuanta empiryczna", legendfontsize = 6, xlim = (1, 1.5))
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_margin, 1))), label = Label2)
Plots.plot!(ecdf(mean.(getindex.(ex3_v302_margin, 1))), label = Label3)
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_margin, 1))), label = Label4, legend = :bottomright)

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

p_h2_ma = Plots.plot(ecdf(margin_abs_300), label = Label1, xlabel = "Średnia marża", ylabel = "Dystrybuanta empiryczna", legendfontsize = 6, xlim = (0, 0.12))
Plots.plot!(ecdf(margin_abs_301), label = Label2)
Plots.plot!(ecdf(margin_abs_302), label = Label3)
Plots.plot!(ecdf(margin_abs_303), label = Label4, legend = :bottomright)

p_h2_p = Plots.plot(ecdf(mean.(getindex.(ex3_v300_price, 1))), label = Label1, xlabel = "Średnia cena", ylabel = "Dystrybuanta empiryczna", legendfontsize = 6, xlim = (0, 1.5))
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_price, 1))), label = Label2)
Plots.plot!(ecdf(mean.(getindex.(ex3_v302_price, 1))), label = Label3)
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_price, 1))), label = Label4, legend = :bottomright)

Plots.savefig(p_h2_m, pwd() *  "\\plots_to_export\\ecdf_margin_perc_w_wo_research.pdf")
Plots.savefig(p_h2_ma, pwd() *  "\\plots_to_export\\ecdf_margin_abs_w_wo_research.pdf")
Plots.savefig(p_h2_p, pwd() *  "\\plots_to_export\\ecdf_price_w_wo_research.pdf")


p_h2_qs = Plots.plot(ecdf(mean.(getindex.(ex3_v300_quantity_sold, 1))), label = Label1, xlabel = "Średni popyt", ylabel = "Dystrybuanta empiryczna", legendfontsize = 6, xlim = (0,30))
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_quantity_sold, 1))), label = Label2)
Plots.plot!(ecdf(mean.(getindex.(ex3_v302_quantity_sold, 1))), label = Label3)
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_quantity_sold, 1))), label = Label4, legend = :bottomright)

p_h2_ss = Plots.plot(ecdf(mean.(getindex.(ex3_v300_quantity_sold,1)) ./ mean.(getindex.(ex3_v300_quantity_produced,1))), xlabel = "Udział sprzedanych dóbr w całości produkcji", ylabel = "Dystrybuanta empiryczna", titlefontsize = 12, label = Label1)
Plots.plot!(ecdf(mean.(getindex.(ex3_v301_quantity_sold,1)) ./ mean.(getindex.(ex3_v301_quantity_produced,1))), label = Label2)
Plots.plot!(ecdf(mean.(getindex.(ex3_v302_quantity_sold,1)) ./ mean.(getindex.(ex3_v302_quantity_produced,1))), label = Label3)
Plots.plot!(ecdf(mean.(getindex.(ex3_v303_quantity_sold,1)) ./ mean.(getindex.(ex3_v303_quantity_produced,1))), label = Label4)

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

DataFrame(Typ = [Label1, Label2, Label3, Label4], Zysk = round.([
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

p_h2_p_h = Plots.scatter(sort(unique(ex3_v300301302303_H)), [transition_percentile(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))[ex3_v300301302303_H .== h]) for h in sort(unique(ex3_v300301302303_H))], xlabel = "Okres przydatności dobra H", ylabel = "Prawdopodobieństwo poniesienia straty, [%]", label = "Producent nie prowadzi badań", markersize = 4)
Plots.scatter!(sort(unique(ex3_v300301302303_H)), [transition_percentile(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))[ex3_v300301302303_H .== h]) for h in sort(unique(ex3_v300301302303_H))], label = "Producent prowadzi badania", markersize = 4)

Plots.savefig(p_h2_p_h, pwd() * "\\plots_to_export\\loss probability vs h.pdf")

7/19-1

####################################### Hipoteza 3 #########################################

"""
Istnieje zależność pomiędzy częstotliwością wymiany informacji pomiędzy konsumentami a poziomem nadzwyczajnego zysku firmy wynikającego z prowadzenia badań konsumenckich
"""

p_h3_sv = Plots.scatter(sort(unique(ex3_v300301302303_nl)), [mean(mean.(ex3_v300_sv)[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], xlabel = "Liczba połączeń", ylabel = "Średnia liczba otrzymanych sygnałów", legend = false, ylim = (0,300))

Plots.savefig(p_h3_sv, pwd() * "\\plots_to_export\\average volume of signals.pdf")

ds300 = getindex.(ex3_v300_quality, 1) .- getindex.(ex3_v300_quality_exp, 1)
ds300b = getindex.(ex3_v300_quality, 1) .- getindex.(ex3_v300_quality_exp_buyers, 1)
ds301 = getindex.(ex3_v301_quality, 1) .- getindex.(ex3_v301_quality_exp, 1)
ds301b = getindex.(ex3_v301_quality, 1) .- getindex.(ex3_v301_quality_exp_buyers, 1)

p_h3_nl = Plots.scatter(sort(unique(ex3_v300301302303_nl)), [mean(mean.(ds300)[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], xlabel = "Liczba połączeń w sieci", ylabel = "Różnica między przekonaniami konsumentów \n a średnim poziomem jakości", label = Label1, titlefontsize = 10, ylim = (0.02, 0.09))
Plots.scatter!(sort(unique(ex3_v300301302303_nl)), [mean(mean.(ds301)[ex3_v300301302303_nl .== nl]) for nl in sort(unique(ex3_v300301302303_nl))], label = Label2)

Plots.savefig(p_h3_nl, pwd() * "\\plots_to_export\\nl vs quality.pdf")

p_h3_l_nr = Plots.scatter(sort(unique(ex3_v300301302303_Li))[2:end], [mean(mean.(ds300)[(ex3_v300301302303_Li .== li) .& (ex3_v300301302303_Lw .== 0.0)]) for li in sort(unique(ex3_v300301302303_Li))][2:end], ylabel = "Różnica między przekonaniami konsumentów, \n a średnim poziomem jakości", xlabel = L"\lambda", label = L"\theta = 0", xlim = (0.09, 1.01), xlabelfontsize = 14, ylabelfontsize = 14, smooth = true, ylim = (0, 0.15), legend = :topleft, legendfontsize = 12)
Plots.scatter!(sort(unique(ex3_v300301302303_Li))[2:end], [mean(mean.(ds300)[(ex3_v300301302303_Li .== li) .& (ex3_v300301302303_Lw .> 0.0)]) for li in sort(unique(ex3_v300301302303_Li))][2:end], label = L"\theta > 0", smooth = true)

p_h3_l_rr = Plots.scatter(sort(unique(ex3_v300301302303_Li))[2:end], [mean(mean.(ds301)[(ex3_v300301302303_Li .== li) .& (ex3_v300301302303_Lw .== 0.0)]) for li in sort(unique(ex3_v300301302303_Li))][2:end], label= L"\theta = 0", ylabel = "Różnica między przekonaniami konsumentów, \n a średnim poziomem jakości", xlabel = L"\lambda",xlim = (0.09, 1.01), smooth = true, ylim = (0, 0.15), legend = :topleft, legendfontsize = 12, xlabelfontsize = 14, ylabelfontsize = 14)
Plots.scatter!(sort(unique(ex3_v300301302303_Li))[2:end], [mean(mean.(ds301)[(ex3_v300301302303_Li .== li) .& (ex3_v300301302303_Lw .> 0.0)]) for li in sort(unique(ex3_v300301302303_Li))][2:end], label = L"\theta > 0", smooth = true)

p_h3_l = Plots.plot(p_h3_l_nr, p_h3_l_rr, size = (1200, 600), margin=5Plots.mm)

Plots.savefig(p_h3_l, pwd() * "\\plots_to_export\\lambda vs quality.pdf")

p_h3_w = Plots.scatter(sort(unique(ex3_v300301302303_Lw)), [mean(mean.(ds300)[ex3_v300301302303_Lw .== lw]) for lw in sort(unique(ex3_v300301302303_Lw))], ylabel = "Różnica między przekonaniami konsumentów \n a średnim poziomem jakości", xlabel = L"\theta", label = Label1, ylim = (0.02, 0.09))
Plots.scatter!(sort(unique(ex3_v300301302303_Lw)), [mean(mean.(ds301)[ex3_v300301302303_Lw .== lw]) for lw in sort(unique(ex3_v300301302303_Lw))], label = Label2)

Plots.savefig(p_h3_w, pwd() * "\\plots_to_export\\theta vs quality.pdf")

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

slope(x,y) = (y[end] - y[1]) / (x[end] - x[1])

slope(sort(unique(ex3_v300301302303_nl)), (π301_hat1 .- π300_hat1)) * 100 / (π301_hat1 .- π300_hat1)[1]
slope(sort(unique(ex3_v300301302303_nl)), (π303_hat1 .- π302_hat1)) * 100 / (π303_hat1 .- π302_hat1)[1]

p_h3_3 = Plots.plot(sort(unique(ex3_v300301302303_nl)), π301_hat1 .- π300_hat1, xlabel = "Liczba połączeń w sieci", ylabel = "Dodatkowy zysk firmy z tytułu \n prowadzenia badań konsumenckich", label = "Wdrożenie badań konsumenckich, gdy konkurent ich nie wykonuje",marker = :circle, markerstrokewidth = 0, ylim = (0, 100))
Plots.plot!(sort(unique(ex3_v300301302303_nl)), π303_hat1 .- π302_hat1, label = "Wdrożenie badań konsumenckich, gdy konkurent też je wykonuje", widen = true, marker = :square, markerstrokewidth = 0)

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

p_h3_lt = Plots.scatter(sort(unique(ex3_v300301302303_Li)), incr_profit_l, markerstrokewidth = 0, label = "Dodatkowy zysk, sygnały indywidualne, " * L"(\lambda)", color = :blue)
Plots.plot!(0.00:0.01:0.99, fit_l, color = "blue", label = "")
Plots.scatter!(sort(unique(ex3_v300301302303_Lw)), incr_profit_t, legend = :topleft, markerstrokewidth = 0, label = "Dodatkowy zysk, sygnały od innych konsumentów, " * L"(\theta)", color = :green)
Plots.plot!(0.00:0.01:0.99, fit_t, color = "green", legend = :topright, xlabel = "Wagi sygnałów " * L"\lambda, \theta", ylabel = "Dodatkowy zysk firmy z tytułu \n prowadzenia badań konsumenckich", label = "")

Plots.savefig(p_h3_lt, pwd() * "\\plots_to_export\\incremental profit vs lambda and theta.pdf")

########## Hipoteza 4 ############

"""
Decyzja producenta o uwzględnianiu wyników badań preferencji konsumentów w prowadzonych przez niego działaniach zależy od ceny pomiaru preferencji oraz horyzontu czasowego planowania producenta
"""

low_nl = ex3_v300301302303_nl .<= 1000

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

payoff1 = [mean(
        [median(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))), median(sum.(getindex.(ex3_v300_producer_surplus_singleton, 2)))
]),
        mean(
            [median(sum.(getindex.(ex3_v301_producer_surplus_singleton, 1))), median(sum.(getindex.(ex3_v302_producer_surplus_singleton, 2)))
]),
        mean(
            [median(sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))), median(sum.(getindex.(ex3_v301_producer_surplus_singleton, 2)))
]),
        mean(
            [median(sum.(getindex.(ex3_v303_producer_surplus_singleton, 1))), median(sum.(getindex.(ex3_v303_producer_surplus_singleton, 2)))
])]


KruskalWallisTest(sum.(getindex.(ex3_v300_producer_surplus_singleton, 1))[low_nl], sum.(getindex.(ex3_v302_producer_surplus_singleton, 1))[low_nl])

payoff2_low_nl = payoff1_low_nl[[1,3,2,4]]
payoff2_high_nl = payoff1_high_nl[[1,3,2,4]]

costs = collect(0:50)

#### High NL

nash_equilibriums_int_high_nl = simulate_ne_for_costs(payoff1_high_nl, payoff2_high_nl, costs, "pNE")
nash_equilibriums_PD_high_nl = simulate_ne_for_costs(payoff1_high_nl, payoff2_high_nl, costs, "PD")
mixedNE_high_nl = simulate_ne_for_costs(payoff1_high_nl, payoff2_high_nl, costs, "mNE")
valid_mixed_NE_high_nl = [all((x .>= 0) .& (x .<= 1)) for x in mixedNE_high_nl]

plot_NE = Plots.plot(costs, true_or_missing.(4 .∈ nash_equilibriums_int_high_nl, k=2), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru przekonań konsumentów",  markercolor = "black", linecolor = "black", label = "Strategia równowagi, strategie czyste")
#Plots.plot!(costs, true_or_missing.(3 .∈ nash_equilibriums_int; k = 3), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
#Plots.plot!(costs, true_or_missing.(2 .∈ nash_equilibriums_int; k = 2), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
Plots.plot!(costs, true_or_missing.(1 .∈ nash_equilibriums_int_high_nl; k = 1), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
Plots.plot!(costs[nash_equilibriums_PD_high_nl .== 1], sum.(nash_equilibriums_int_high_nl)[nash_equilibriums_PD_high_nl .== 1], markercolor = "red", markershape = :square, markerstrokewidth = 0, markersize = 6, linecolor = "red", label = "Rozwiązanie będące dylematem więźnia")
Plots.plot!(costs[valid_mixed_NE_high_nl], fill(1, count(valid_mixed_NE_high_nl)), markercolor = "blue", markershape = :square, markersize = 6, markerstrokewidth = 0, label = "Dodatkowe rozwiązanie w strategiach mieszanych")
Plots.plot!(costs[valid_mixed_NE_high_nl], fill(2, count(valid_mixed_NE_high_nl)), markercolor = "blue", markershape = :square, markersize = 6, markerstrokewidth = 0, label = "", legend = :topright)
Plots.plot!(yticks = ([1,2], ["Żaden producent nie prowadzi badań", "Obaj producenci prowadzą badania"]), size = (1000, 200))
xaxis!(xlabel = "Koszt pomiaru preferencji")

p_high_nl = Plots.plot(payoff1_high_nl[[1,2,4,3,1]], payoff2_high_nl[[1,2,4,3,1]], legend = :outerright, marker = :circle, label = "C = 0", xlabel = "Zysk firmy #1", ylabel = "Zysk firmy #2", xlim = (0,80), ylim = (0,80))
annotate!.(payoff1_high_nl[[1,2,4,3]] .+ [-2,2,2,-2], payoff2_high_nl[[1,2,4,3]] .+ [-2,2,2,2], Plots.text.(["(N,N)", "(R,N)", "(R,R)", "(N,R)"], fill(8,4), :blue), alpha = 0.5)

c = 15; Plots.plot!(payoff1_high_nl[[1,2,4,3,1]] .- [0, c, c, 0, 0], payoff2_high_nl[[1,2,4,3,1]] .- [0, 0, c, c, 0], marker = :circle, label = "C = " * string(c))
annotate!.(payoff1_high_nl[[2,4,3]] .+ [2,2,-2] .- [c, c, 0], payoff2_high_nl[[2,4,3]] .+ [2,2,-2] .- [0, c, c], Plots.text.(["(R,N)", "(R,R)", "(N,R)"], fill(8,3), :red), linealpha = 0.5)

c = 30; Plots.plot!(payoff1_high_nl[[1,2,4,3,1]] .- [0, c, c, 0, 0], payoff2_high_nl[[1,2,4,3,1]] .- [0, 0, c, c, 0], marker = :circle, label = "C = " * string(c))
annotate!.(payoff1_high_nl[[2,4,3]] .+ [2,2,-2] .- [c, c, 0], payoff2_high_nl[[2,4,3]] .+ [2,2,-2] .- [0, c, c], Plots.text.(["(R,N)", "(R,R)", "(N,R)"], fill(8,3), :green), linealpha = 0.5)

c = 45; Plots.plot!(payoff1_high_nl[[1,2,4,3,1]] .- [0, c, c, 0, 0], payoff2_high_nl[[1,2,4,3,1]] .- [0, 0, c, c, 0], marker = :circle, label = "C = " * string(c))
annotate!.(payoff1_high_nl[[2,4,3]] .+ [2,2,-2] .- [c, c, 0], payoff2_high_nl[[2,4,3]] .+ [2,2,-2] .- [0, c, c], Plots.text.(["(R,N)", "(R,R)", "(N,R)"], fill(8,3), :purple), linealpha = 0.5)

Plots.savefig(p_high_nl, pwd() * "\\plots_to_export\\nash_equilibriums_h.pdf")
Plots.savefig(p_high_nl, pwd() * "\\plots\\h4\\nash_equilibriums_h.svg")


#### Low NL

nash_equilibriums_int_low_nl = simulate_ne_for_costs(payoff1_low_nl, payoff2_low_nl, costs, "pNE")
nash_equilibriums_PD_low_nl = simulate_ne_for_costs(payoff1_low_nl, payoff2_low_nl, costs, "PD")
mixedNE_low_nl = simulate_ne_for_costs(payoff1_low_nl, payoff2_low_nl, costs, "mNE")
valid_mixed_NE_low_nl = [all((x .>= 0) .& (x .<= 1)) for x in mixedNE_low_nl]

plot_NE = Plots.plot(costs, true_or_missing.(4 .∈ nash_equilibriums_int_low_nl, k=2), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru przekonań konsumentów",  markercolor = "black", linecolor = "black", label = "Strategia równowagi, strategie czyste")
#Plots.plot!(costs, true_or_missing.(3 .∈ nash_equilibriums_int; k = 3), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
#Plots.plot!(costs, true_or_missing.(2 .∈ nash_equilibriums_int; k = 2), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
Plots.plot!(costs, true_or_missing.(1 .∈ nash_equilibriums_int_low_nl; k = 1), markershape = :square, markerstrokewidth = 0, markersize = 4, xlabel = "Koszt pomiaru oczekiwań konsumentów", label = "", markercolor = "black", linecolor = "black")
Plots.plot!(costs[nash_equilibriums_PD_low_nl .== 1], sum.(nash_equilibriums_int_low_nl)[nash_equilibriums_PD_low_nl .== 1], markercolor = "red", markershape = :square, markerstrokewidth = 0, markersize = 6, linecolor = "red", label = "Rozwiązanie będące dylematem więźnia")
Plots.plot!(costs[valid_mixed_NE_low_nl], fill(1, count(valid_mixed_NE_low_nl)), markercolor = "blue", markershape = :square, markersize = 6, markerstrokewidth = 0, label = "Dodatkowe rozwiązanie w strategiach mieszanych")
Plots.plot!(costs[valid_mixed_NE_low_nl], fill(2, count(valid_mixed_NE_low_nl)), markercolor = "blue", markershape = :square, markersize = 6, markerstrokewidth = 0, label = "", legend = :topright)
Plots.plot!(yticks = ([1,2], ["Żaden producent nie prowadzi badań", "Obaj producenci prowadzą badania"]), size = (1000, 200))
xaxis!(xlabel = "Koszt pomiaru preferencji")

p_low_nl = Plots.plot(payoff1_low_nl[[1,2,4,3,1]], payoff2_low_nl[[1,2,4,3,1]], legend = :outerright, marker = :circle, label = "C = 0", xlabel = "Zysk firmy #1", ylabel = "Zysk firmy #2", xlim = (0,80), ylim = (0,80))
annotate!.(payoff1_low_nl[[1,2,4,3]] .+ [-2,2,2,-2], payoff2_low_nl[[1,2,4,3]] .+ [-2,2,2,2], Plots.text.(["(N,N)", "(R,N)", "(R,R)", "(N,R)"], fill(8,4), :blue), alpha = 0.5)

c = 15; Plots.plot!(payoff1_low_nl[[1,2,4,3,1]] .- [0, c, c, 0, 0], payoff2_low_nl[[1,2,4,3,1]] .- [0, 0, c, c, 0], marker = :circle, label = "C = " * string(c))
annotate!.(payoff1_low_nl[[2,4,3]] .+ [2,2,-2] .- [c, c, 0], payoff2_low_nl[[2,4,3]] .+ [2,2,-2] .- [0, c, c], Plots.text.(["(R,N)", "(R,R)", "(N,R)"], fill(8,3), :red), linealpha = 0.5)

c = 30; Plots.plot!(payoff1_low_nl[[1,2,4,3,1]] .- [0, c, c, 0, 0], payoff2_low_nl[[1,2,4,3,1]] .- [0, 0, c, c, 0], marker = :circle, label = "C = " * string(c))
annotate!.(payoff1_low_nl[[2,4,3]] .+ [2,2,-2] .- [c, c, 0], payoff2_low_nl[[2,4,3]] .+ [2,2,-2] .- [0, c, c], Plots.text.(["(R,N)", "(R,R)", "(N,R)"], fill(8,3), :green), linealpha = 0.5)

c = 45; Plots.plot!(payoff1_low_nl[[1,2,4,3,1]] .- [0, c, c, 0, 0], payoff2_low_nl[[1,2,4,3,1]] .- [0, 0, c, c, 0], marker = :circle, label = "C = " * string(c))
annotate!.(payoff1_low_nl[[2,4,3]] .+ [2,2,-2] .- [c, c, 0], payoff2_low_nl[[2,4,3]] .+ [2,2,-2] .- [0, c, c], Plots.text.(["(R,N)", "(R,R)", "(N,R)"], fill(8,3), :purple), linealpha = 0.5)

p_low_nl

Plots.savefig(p_low_nl, pwd() * "\\plots_to_export\\nash_equilibriums_l.pdf")
Plots.savefig(p_low_nl, pwd() * "\\plots\\h4\\nash_equilibriums_l.svg")

############## Additional analysis for H3, new DF

@load "C:\\Users\\User\\Documents\\PhDWorkspace_more_nl_varianteps.jld2"

minimum(profit_diff_302303)
argmin(profit_diff_302303)
ex3_v300301302303_ϵq[argmin(profit_diff_302303)]

profit_diff_300301 = sum.(getindex.(ex3_v301_producer_surplus_singleton, 2)) .- sum.(getindex.(ex3_v300_producer_surplus_singleton, 2))

profit_diff_302303 = sum.(getindex.(ex3_v303_producer_surplus_singleton, 2)) .- sum.(getindex.(ex3_v302_producer_surplus_singleton, 2))

p_3_add = Plots.scatter(sort(unique(ex3_v300301302303_ϵq)), [mean(profit_diff_300301[ex3_v300301302303_ϵq .== eq]) for eq in sort(unique(ex3_v300301302303_ϵq))], xlabel = "Wariancja jakości produktów", ylabel = "Dodatkowy zysk \n z prowadzenia badań konsumenckich", label = "Firma konkurencyjna nie prowadzi badań", smooth = true, xlim = (0,0.21))
Plots.scatter!(sort(unique(ex3_v300301302303_ϵq)), [mean(profit_diff_302303[ex3_v300301302303_ϵq .== eq]) for eq in sort(unique(ex3_v300301302303_ϵq))], label = "Firma konkurencyjna prowadzi badania", smooth = true)

Plots.savefig(p_3_add, pwd() * "\\plots_to_export\\quality_variance_profit.pdf")

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