include(pwd() * "\\methods\\methods.jl")
@load "C:\\Users\\User\\Documents\\PhDWorkspace_more_nl.jld2"

# Potwierdzenie istnienia planned obsolence

xd = collect(0.10:0.10:0.90)
yd = collect(0.10:0.10:0.90)

mat = zeros(Int64, lastindex(xd) - 1, lastindex(yd)-1)
iter0 = 100

for i in 1:(lastindex(xd)-1)
    for j in 1:(lastindex(yd)-1)
        mat[i,j] = count((xd[i] .<= mean([getindex.(getindex.(ex3_no_sm_durability, 1), m) for m in iter0:500]) .< xd[i+1]) .& (yd[j] .<= mean([getindex.(getindex.(ex3_no_sm_durability, 2), m) for m in iter0:500]) .< yd[j+1]))
    end
end

val_p1 = StatsPlots.heatmap(xd, yd, mat', c = cgrad(:roma, scale = :exp, rev = true), fill_z=mat, aspect_ratio=:equal, xlabel = "Trwałość producenta pierwszego", ylabel = "Trwałość producenta drugiego", xlim = (0.10, 0.90), ylim = (0.10, 0.90), colorbar_title = "") 

nrow, ncol = size(mat)
ann = [(xd[i] + 0.05,yd[j]+0.05, Plots.text(round(mat[i,j] / sum(mat), digits=3), 8, :black, :center))
            for i in 1:nrow for j in 1:ncol]
annotate!(ann, linecolor=:white)

val_p1

Plots.savefig(val_p1, pwd() * "\\plots_to_export\\val_planned_obsolence.pdf")

##  ##

ex3_no_sm_total_surplus = Vector{Vector{Float64}}()
ex3_no_sm_producer_surplus = Vector{Vector{Float64}}()
ex3_no_sm_consumer_surplus = Vector{Vector{Float64}}()
ex3_no_sm_price = Vector{Vector{Vector{Float64}}}()
ex3_no_sm_quantity_produced = Vector{Vector{Vector{Float64}}}()
ex3_no_sm_quantity_sold = Vector{Vector{Vector{Float64}}}()
ex3_no_sm_reselling = Vector{Vector{Vector{Float64}}}()
ex3_no_sm_producer_surplus_singleton = Vector{Vector{Vector{Float64}}}()
ex3_no_sm_buying_history = []
ex3_no_sm_quality = Vector{Vector{Vector{Float64}}}()
ex3_no_sm_durability = Vector{Vector{Vector{Float64}}}()
ex3_no_sm_margin = Vector{Vector{Vector{Float64}}}()
ex3_no_sm_quality_exp = Vector{Vector{Vector{Float64}}}()
ex3_no_sm_durability_exp = Vector{Vector{Vector{Float64}}}()
ex3_no_sm_sv = Vector{Vector{Float64}}()
ex3_no_sm_saturation = Vector{Vector{Vector{Int64}}}()
ex3_no_sm_quality_exp_buyers = []
ex3_no_sm_durability_exp_buyers = []

ex3_is_sm_total_surplus = Vector{Vector{Float64}}()
ex3_is_sm_producer_surplus = Vector{Vector{Float64}}()
ex3_is_sm_consumer_surplus = Vector{Vector{Float64}}()
ex3_is_sm_price = Vector{Vector{Vector{Float64}}}()
ex3_is_sm_quantity_produced = Vector{Vector{Vector{Float64}}}()
ex3_is_sm_quantity_sold = Vector{Vector{Vector{Float64}}}()
ex3_is_sm_reselling = Vector{Vector{Vector{Float64}}}()
ex3_is_sm_producer_surplus_singleton = Vector{Vector{Vector{Float64}}}()
ex3_is_sm_buying_history = []
ex3_is_sm_quality = Vector{Vector{Vector{Float64}}}()
ex3_is_sm_durability = Vector{Vector{Vector{Float64}}}()
ex3_is_sm_margin = Vector{Vector{Vector{Float64}}}()
ex3_is_sm_quality_exp = Vector{Vector{Vector{Float64}}}()
ex3_is_sm_durability_exp = Vector{Vector{Vector{Float64}}}()
ex3_is_sm_sv = Vector{Vector{Float64}}()
ex3_is_sm_saturation = Vector{Vector{Vector{Int64}}}()
ex3_is_sm_quality_exp_buyers = []
ex3_is_sm_durability_exp_buyers = []

ex3_no_sm301302303_H = Vector{Int64}()
ex3_no_sm301302303_Li = Vector{Float64}()
ex3_no_sm301302303_Lw = Vector{Float64}()
ex3_no_sm301302303_sd = Vector{Int64}()
ex3_no_sm301302303_cu = Vector{Float64}()
ex3_no_sm301302303_sm = Vector{Float64}()
ex3_no_sm301302303_nl = Vector{Int64}()
ex3_no_sm301302303_ρm = Vector{Float64}()
ex3_no_sm301302303_cc = Vector{Float64}()
ex3_no_sm301302303_m = Vector{Float64}()

for i in 1:500

    println(i)

    N = 200

    li = sample(0.0:0.1:1.0)
    lw = sample(0.0:0.1:1.0)
    h = sample(3:8)
    sd = sample(1:100000)
    cu = sample(0.0:0.1:0.5)
    sm = sample(0.10:0.05:0.30)
    nl = sample(N:N:(15*N))
    ρm = rand(Uniform(0.5, 0.95))
    cc = sample(0.4:0.1:0.6)
    m = sample(1.1:0.1:1.3)

    Random.seed!(sd)

    println([li, lw, h, cu, sm, nl, ρm, cc, m])

    ex3_no_sm_sim = TO_GO(500, 2, N, nl, [cc, cc], [m, m], "random", li, lw, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[1., 2.], [1., 2.]], cu, false, 0, [ρm, 1.], "softmax", ["internal knowledge", "internal knowledge"], [0., 0.], h, "dist", TriangularDist(0,1,0.5))

    Random.seed!(sd)

    ex3_is_sm_sim = TO_GO(500, 2, N, nl, [cc, cc], [m, m], "random", li, lw, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[1., 2.], [1., 2.]], cu, true, 0, [ρm, 1.], "softmax", ["internal knowledge", "internal knowledge"], [0., 0.], h, "dist", TriangularDist(0,1,0.5))

    push!(ex3_no_sm301302303_H, h)
    push!(ex3_no_sm301302303_Li, li)
    push!(ex3_no_sm301302303_Lw, lw)
    push!(ex3_no_sm301302303_sd, sd)
    push!(ex3_no_sm301302303_cu, cu)
    push!(ex3_no_sm301302303_sm, sm)
    push!(ex3_no_sm301302303_nl, nl)
    push!(ex3_no_sm301302303_ρm, ρm)
    push!(ex3_no_sm301302303_cc, cc)
    push!(ex3_no_sm301302303_m, m)

    push!(ex3_no_sm_total_surplus, calculate_surplus(ex3_no_sm_sim, "total", false))
    push!(ex3_no_sm_producer_surplus, calculate_surplus(ex3_no_sm_sim, "producer", false))
    push!(ex3_no_sm_consumer_surplus, calculate_surplus(ex3_no_sm_sim, "consumer,total",false))
    push!(ex3_no_sm_quality, getfield.(ex3_no_sm_sim.sellers, :quality_history))
    push!(ex3_no_sm_durability, getfield.(ex3_no_sm_sim.sellers, :durability_history))
    push!(ex3_no_sm_margin, getfield.(ex3_no_sm_sim.sellers, :margin_history))
    push!(ex3_no_sm_price, calculate_price_history.(ex3_no_sm_sim.sellers; product_life = h))
    push!(ex3_no_sm_quantity_produced, getfield.(ex3_no_sm_sim.sellers, :quantity_produced_history))
    push!(ex3_no_sm_quantity_sold, getfield.(ex3_no_sm_sim.sellers, :quantity_sold_history))
    push!(ex3_no_sm_producer_surplus_singleton, calculate_profit_history.(ex3_no_sm_sim.sellers))
    push!(ex3_no_sm_reselling, getfield.(ex3_no_sm_sim.sellers, :reselling_history))
    push!(ex3_no_sm_quality_exp, [mean([getindex.(x,y) for x in getfield.(ex3_no_sm_sim.buyers, :quality_expectation_history)]) for y in 1:2])
    push!(ex3_no_sm_durability_exp, [mean([getindex.(x,y) for x in getfield.(ex3_no_sm_sim.buyers, :durability_expectation_history)]) for y in 1:2])
    push!(ex3_no_sm_sv, getfield.(ex3_no_sm_sim.buyers, :signal_volume))
    push!(ex3_no_sm_saturation, [sum([any_vec(x,s) for x in getfield.(ex3_no_sm_sim.buyers, :unit_possessed_history)]) for s in 1:2])
    push!(ex3_no_sm_quality_exp_buyers, get_expectation_buyers(ex3_no_sm_sim.buyers, :quality_expectation_history; s = 2, T = 500))
    push!(ex3_no_sm_durability_exp_buyers, get_expectation_buyers(ex3_no_sm_sim.buyers, :durability_expectation_history; s = 2, T = 500))

    push!(ex3_is_sm_total_surplus, calculate_surplus(ex3_is_sm_sim, "total", false))
    push!(ex3_is_sm_producer_surplus, calculate_surplus(ex3_is_sm_sim, "producer", false))
    push!(ex3_is_sm_consumer_surplus, calculate_surplus(ex3_is_sm_sim, "consumer,total",false))
    push!(ex3_is_sm_quality, getfield.(ex3_is_sm_sim.sellers, :quality_history))
    push!(ex3_is_sm_durability, getfield.(ex3_is_sm_sim.sellers, :durability_history))
    push!(ex3_is_sm_margin, getfield.(ex3_is_sm_sim.sellers, :margin_history))
    push!(ex3_is_sm_price, calculate_price_history.(ex3_is_sm_sim.sellers; product_life = h))
    push!(ex3_is_sm_quantity_produced, getfield.(ex3_is_sm_sim.sellers, :quantity_produced_history))
    push!(ex3_is_sm_quantity_sold, getfield.(ex3_is_sm_sim.sellers, :quantity_sold_history))
    push!(ex3_is_sm_producer_surplus_singleton, calculate_profit_history.(ex3_is_sm_sim.sellers))
    push!(ex3_is_sm_reselling, getfield.(ex3_is_sm_sim.sellers, :reselling_history))
    push!(ex3_is_sm_quality_exp, [mean([getindex.(x,y) for x in getfield.(ex3_is_sm_sim.buyers, :quality_expectation_history)]) for y in 1:2])
    push!(ex3_is_sm_durability_exp, [mean([getindex.(x,y) for x in getfield.(ex3_is_sm_sim.buyers, :durability_expectation_history)]) for y in 1:2])
    push!(ex3_is_sm_sv, getfield.(ex3_is_sm_sim.buyers, :signal_volume))
    push!(ex3_is_sm_saturation, [sum([any_vec(x,s) for x in getfield.(ex3_is_sm_sim.buyers, :unit_possessed_history)]) for s in 1:2])
    push!(ex3_is_sm_quality_exp_buyers, get_expectation_buyers(ex3_is_sm_sim.buyers, :quality_expectation_history; s = 2, T = 500))
    push!(ex3_is_sm_durability_exp_buyers, get_expectation_buyers(ex3_is_sm_sim.buyers, :durability_expectation_history; s = 2, T = 500))

end

"""jldsave("C:\\Users\\User\\Documents\\PhDWorkspace_validation_secondary_market.jld2"; ex3_no_sm_total_surplus,
ex3_no_sm_producer_surplus,
ex3_no_sm_consumer_surplus,
ex3_no_sm_price,
ex3_no_sm_quantity_produced,
ex3_no_sm_quantity_sold,
ex3_no_sm_reselling,
ex3_no_sm_producer_surplus_singleton,
ex3_no_sm_quality,
ex3_no_sm_durability,
ex3_no_sm_margin,
ex3_no_sm_quality_exp,
ex3_no_sm_durability_exp,
ex3_no_sm_sv, 
ex3_no_sm_saturation, 
ex3_no_sm_quality_exp_buyers,
ex3_no_sm_durability_exp_buyers,

ex3_is_sm_total_surplus,
ex3_is_sm_producer_surplus,
ex3_is_sm_consumer_surplus,
ex3_is_sm_price,
ex3_is_sm_quantity_produced,
ex3_is_sm_quantity_sold,
ex3_is_sm_reselling,
ex3_is_sm_producer_surplus_singleton,
ex3_is_sm_quality,
ex3_is_sm_durability,
ex3_is_sm_margin,
ex3_is_sm_quality_exp,
ex3_is_sm_durability_exp,
ex3_is_sm_sv,
ex3_is_sm_saturation, 
ex3_is_sm_quality_exp_buyers,
ex3_is_sm_durability_exp_buyers)"""


# rynek wtórny zwiększa całkowitą nadwyżkę Ghose 2010

val_p2 = Plots.plot(ecdf(sum.(ex3_no_sm_total_surplus)), xlabel = "Nadwyżka całkowita", ylabel = "Dystrybuanta empiryczna", label = "Rynek wtórny nie istnieje")
Plots.plot!(ecdf(sum.(ex3_is_sm_total_surplus)), label = "Rynek wtórny  istnieje", legend = :bottomright)

Plots.savefig(val_p2, pwd() * "\\plots_to_export\\val_sec_mar_total_surplus.pdf")

# mało iteracji, niekonkluzywne - jeśli y2 < y1 to potwierdzenie wyników Yin oraz Chen

Plots.plot(ecdf(sum.(ex3_no_sm_producer_surplus)))
Plots.plot!(ecdf(sum.(ex3_is_sm_producer_surplus)))

# zmniejsza liczbę sprzedawanych sztuk Esteban, Shum

Plots.plot(ecdf(mean.(getindex.(ex3_no_sm_quantity_sold, 1))))
Plots.plot!(ecdf(mean.(getindex.(ex3_is_sm_quantity_sold, 1))))

# pozwala podwyższyć cenę Hendel, Lizzeri

val_p3 = Plots.plot(ecdf(mean.(getindex.(ex3_no_sm_price, 1))), xlabel = "Średnia cena", ylabel = "Dystrybuanta empiryczna", label = "Rynek wtórny nie istnieje")
Plots.plot!(ecdf(mean.(getindex.(ex3_is_sm_price, 1))), label = "Rynek wtórny  istnieje", legend = :bottomright)

Plots.savefig(val_p3, pwd() * "\\plots_to_export\\val_sec_mar_prices.pdf")

# adresuje potrzeby osób ze środka rozkładu Ghose 2010

Random.seed!(2)

@time sim_single = TO_GO(200, 2, 500, 1000, [0.6, 0.6], [1.1, 1.1], "random", .5, .5, "deterministic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[1.,2.], [1., 2.]], 0.1, true, 1, [0.7, 1.], "softmax", ["internal knowledge", "internal knowledge"], [0., 0.], 6, "dist", Uniform(0,1))

colors = fill("grey", lastindex(sim_single.buyers))
alphas = fill(0.25, lastindex(sim_single.buyers))
colors[[any(getfield.(x, :decision) .== "buy, primary market") for x in getfield.(sim_single.buyers, :unit_buying_selling_history)]] .= "red"
alphas[[any(getfield.(x, :decision) .== "buy, primary market") for x in getfield.(sim_single.buyers, :unit_buying_selling_history)]] .= 1
colors[[any(getfield.(x, :decision) .== "buy, secondary market") & all(getfield.(x, :decision) .!= "buy, primary market") for x in getfield.(sim_single.buyers, :unit_buying_selling_history)]] .= "green"
alphas[[any(getfield.(x, :decision) .== "buy, secondary market") & all(getfield.(x, :decision) .!= "buy, primary market") for x in getfield.(sim_single.buyers, :unit_buying_selling_history)]] .= 1

val_p4 = Plots.scatter((1:lastindex(sim_single.buyers)) ./ lastindex(sim_single.buyers), sort(getfield.(sim_single.buyers, :std_reservation_price), rev = true), color = colors[sortperm(getfield.(sim_single.buyers, :std_reservation_price), rev = true)], markerstrokewidth = 0, markeralpha = alphas[sortperm(getfield.(sim_single.buyers, :std_reservation_price), rev = true)], xlabel = "% populacji konsumentów", ylabel = "Cena rezerwacji " * L"\beta_i", label = "")
Plots.scatter!(1, 1, color = "red", label = "Kupno na rynku pierwotnym", markerstrokewidth = 0)
Plots.scatter!(1, 1, color = "green", label = "Kupno tylko na rynku wtórnym", markerstrokewidth = 0)
Plots.scatter!(1, 1, color = "grey", alpha = 0.25, label = "Brak zakupów", markerstrokewidth = 0)

Plots.savefig(val_p4, pwd() * "\\plots_to_export\\val_sec_mar_std_res_prc.pdf")

# odniesienie do modelu Izquierdo
# działa jeśli ϵ = 0.20

@time sim_single = TO_GO(10000, 2, 200, 0, [0.6, 0.6], [1., 1.], "random", .85, .0, "deterministic", [[0.4, 0.6], [0.4, 0.6]], [[0.4, 0.6], [0.4, 0.6]], [[.975, 1.025], [.975, 1.025]], 0.00, false, 1, [0.99, 1.], "softmax", ["internal knowledge", "internal knowledge"], [0., 0.], 4, "dist", Uniform(0,1))

val_p5 = Plots.plot(sum(getfield.(sim_single.sellers, :quantity_sold_history)), xscale = :log, xlabel = "Czas", ylabel = "Liczba sprzedanych sztuk", label = "")

Plots.savefig(val_p5, pwd() * "\\plots_to_export\\val_izquierdo_individual_learning.pdf")


