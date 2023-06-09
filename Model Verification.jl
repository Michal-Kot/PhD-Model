##################################### Unit testing ####################################################

include(pwd() * "\\methods\\methods.jl")

########################################### Functions 

function create_buying(_buyer, product_id)
    b = []
    for item in _buyer.unit_buying_selling_history
        if item.decision == "buy, primary market"
            if item.product == product_id
                push!(b, item.t)
            end
        end
    end
    return b
end

function get_signal_value(_buyer, product_id, qs)
    b = Vector{Union{Missing, Float64}}(undef, lastindex(_buyer.quality_expectation_history))
    for item in _buyer.unit_buying_selling_history
        if item.decision == "buy, primary market"
            if item.product == product_id
                b[item.t] = item.real_product_features[qs]
            end
        end
    end
    return b
end

function create_destroy(_buyer, product_id)
    b = []
    for item in _buyer.unit_buying_selling_history
        if item.decision == "destroyed"
            if item.product == product_id
                push!(b, item.t)
            end
        end
    end
    return b
end

sim_single.buyers[193].unit_buying_selling_history

function create_resell(_buyer, product_id)
    b = []
    for item in _buyer.unit_buying_selling_history
        if item.decision == "sell, secondary market"
            if item.product == product_id
                push!(b, item.t)
            end
        end
    end
    return b
end

function create_failed_resell(_buyer, product_id)
    b = []
    for item in _buyer.unit_buying_selling_history
        if item.decision == "failed to resell"
            if item.product == product_id
                push!(b, item.t)
            end
        end
    end
    return b
end

function check_impact_of_signals(sim_single, buyer_id, product_id)

    _buyer = sim_single.buyers[buyer_id]

    own_signals = []
    outer_signals = []

    for item in _buyer.unit_buying_selling_history
        if item.decision == "buy, primary market"
            if item.product == product_id
                push!(own_signals, _buyer.quality_expectation_history[item.t][product_id] - _buyer.quality_expectation_history[item.t - 1][product_id])
            end
        end
    end 

    for _buyer in sim_single.buyers[Not(buyer_id)]
        for item in _buyer.unit_buying_selling_history
            if item.decision == "buy, primary market"
                if item.product == product_id
                    push!(outer_signals, _buyer.quality_expectation_history[item.t][product_id] - _buyer.quality_expectation_history[item.t - 1][product_id])
                end
            end
        end 
    end

    return (own_signals, outer_signals)

end

function quality_sampled_products(buyers, product_id, qd)

    quality_of_unit_possessed_history = [[] for _ in 1:lastindex(buyers[1].quality_expectation_history)]

    for _buyer in buyers
        for item in _buyer.unit_buying_selling_history
            if item.decision == "buy, primary market"
                if item.product == product_id
                    push!(quality_of_unit_possessed_history[item.t], item.real_product_features[qd])
                end
            end
        end
    end

    return quality_of_unit_possessed_history

end

function product_age(buyers)

    # t = buy, t+1 = destroy

    buy_to_destroy = []

    # t = buy, t+1 = resell

    buy_to_resell = []

    # t = buy, sec mar, t+1 = destroy

    rebuy_to_destroy = []

    # t = buy, sec mar, t+1 = buy, pri mar

    rebuy_to_buy = []

    for _buyer in buyers
        for idx in 2:lastindex(_buyer.unit_buying_selling_history)
            if (_buyer.unit_buying_selling_history[idx-1].decision == "buy, primary market") & (_buyer.unit_buying_selling_history[idx].decision == "destroyed")
                push!(buy_to_destroy, _buyer.unit_buying_selling_history[idx].t - _buyer.unit_buying_selling_history[idx-1].t)
            end

            if (_buyer.unit_buying_selling_history[idx-1].decision == "buy, primary market") & (_buyer.unit_buying_selling_history[idx].decision == "buy, primary market")
                push!(buy_to_resell, _buyer.unit_buying_selling_history[idx].t - _buyer.unit_buying_selling_history[idx-1].t)
            end

            if (_buyer.unit_buying_selling_history[idx-1].decision == "buy, secondary market") & (_buyer.unit_buying_selling_history[idx].decision == "destroyed")
                push!(rebuy_to_destroy, _buyer.unit_buying_selling_history[idx].t - _buyer.unit_buying_selling_history[idx-1].t)
            end

            if (_buyer.unit_buying_selling_history[idx-1].decision == "buy, secondary market") & (_buyer.unit_buying_selling_history[idx].decision == "buy, primary market")
                push!(rebuy_to_buy, _buyer.unit_buying_selling_history[idx].t - _buyer.unit_buying_selling_history[idx-1].t)
            end

        end
    end

    return [buy_to_destroy, buy_to_resell, rebuy_to_destroy, rebuy_to_buy]

end

############################## Individual level ###############################################

########################################### If signals change expectation

Random.seed!(12345678)

sim_single = TO_GO(100, 2, 2, 1, [0.4, 0.4], [1.0, 1.0], "random", 0.5, 0.5, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.5], 4, "dist", TriangularDist(0,1,0.5))

## Number of signals

sim_single.buyers[2].unit_buying_selling_history
sim_single.buyers[1].signal_volume == countmap(getfield.(sim_single.buyers[2].unit_buying_selling_history, :decision))["buy, primary market"]

sim_single.buyers[1].unit_buying_selling_history
sim_single.buyers[2].signal_volume == countmap(getfield.(sim_single.buyers[1].unit_buying_selling_history, :decision))["buy, primary market"]

p3 = Plots.plot(sim_single.sellers[1].quality_history, label = "Średnia jakość", xlabel = "Czas", ylabel = "Jakość", legend = :outerright)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 1), label = "Przekonanie o jakości " * L"K_{ijt}")
Plots.scatter!(get_signal_value(sim_single.buyers[1], 1, 1), label = "Rzeczywista jakość " * L"\hat K_{ijt}" )
Plots.scatter!(get_signal_value(sim_single.buyers[2], 1, 1), label = "Rzeczywista jakość " * L"\bar K_{ijt}")

Plots.savefig(p3, pwd() * "\\plots_to_export\\verification_Ehatbar.pdf")

Plots.plot!(create_buying(sim_single.buyers[1], 1) .- 0.9, seriestype = :vline, linestyle = :dash)
Plots.plot!(create_buying(sim_single.buyers[2], 1) .- 1.1, seriestype = :vline)

Plots.plot(sim_single.sellers[2].quality_history)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 2))
Plots.plot!(create_buying(sim_single.buyers[1], 2) .+ 0.1, seriestype = :vline, 
linestyle = :dash)
Plots.plot!(create_buying(sim_single.buyers[2], 2) .- 0.1, seriestype = :vline, 
linestyle = :dash)

########################################### Check if 0 weight of comms results in no difference in expectations

## EXP = 0, WOM = 1

Random.seed!(1425)

sim_single = TO_GO(50, 2, 2, 1, [0.4, 0.4], [1.0, 1.0], "random", 0.0, 1.0, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.5], 4, "dist", TriangularDist(0,1,0.5))

prod(x) = any(x .> 0) ? argmax(x) : 0

Plots.scatter(1:50, prod.(sim_single.buyers[1].unit_possessed_history), color = :black)
Plots.plot!(create_buying(sim_single.buyers[1], 1), seriestype = :vline, color = :red)
Plots.plot!(create_buying(sim_single.buyers[1], 2), seriestype = :vline, color = :green)
Plots.plot!(create_destroy(sim_single.buyers[1], 1), seriestype = :vline, color = :red, linestyle = :dash)
Plots.plot!(create_destroy(sim_single.buyers[1], 2), seriestype = :vline, color = :green, linestyle = :dash)

[any(x .== 0) for x in check_impact_of_signals(sim_single, 1, 1)]

p1 = Plots.plot(sim_single.sellers[1].quality_history, label = "Średnia jakość dobra", xlabel = "Czas", ylabel = "Jakość", legend = :outerright)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 1), label = "Przekonanie o jakości " * L"K_{ijt}")
Plots.plot!(create_buying(sim_single.buyers[1], 1) .+ 0.1, seriestype = :vline, linestyle = :dash, label = "Własne doświadczenie")
Plots.plot!(create_buying(sim_single.buyers[2], 1) .- 0.1, seriestype = :vline, linestyle = :dash, label = "Wymiana informacji")

Plots.savefig(p1, pwd() * "\\plots_to_export\\verification_lambda0.pdf")

Plots.plot(sim_single.sellers[2].quality_history)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 2))
Plots.plot!(create_buying(sim_single.buyers[1], 2) .+ 0.1, seriestype = :vline, linestyle = :dash)
Plots.plot!(create_buying(sim_single.buyers[2], 2) .- 0.1, seriestype = :vline, 
linestyle = :dash)

[any(x .== 0) for x in check_impact_of_signals(sim_single, 1, 2)]

## EXP = 1, WOM = 0

Random.seed!(15512)

sim_single = TO_GO(50, 2, 2, 1, [0.4, 0.4], [1.0, 1.0], "random", 1.0, 0.0, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.5], 4, "dist", TriangularDist(0,1,0.5))

p2 = Plots.plot(sim_single.sellers[1].quality_history, label = "Średnia jakość dobra", xlabel = "Czas", ylabel = "Jakość", legend = :outerright)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 1), label = "Przekonanie o jakości " * L"K_{ijt}")
Plots.plot!(create_buying(sim_single.buyers[1], 1) .+ 0.1, seriestype = :vline, linestyle = :dash, label = "Własne doświadczenie")
Plots.plot!(create_buying(sim_single.buyers[2], 1) .- 0.1, seriestype = :vline, linestyle = :dash, label = "Wymiana informacji")


Plots.savefig(p2, pwd() * "\\plots_to_export\\verification_theta0.pdf")

Plots.scatter(1:50, prod.(sim_single.buyers[1].unit_possessed_history), color = :black)
Plots.plot!(create_buying(sim_single.buyers[1], 1), seriestype = :vline, color = :red)
Plots.plot!(create_buying(sim_single.buyers[1], 2), seriestype = :vline, color = :green)
Plots.plot!(create_destroy(sim_single.buyers[1], 1), seriestype = :vline, color = :red, linestyle = :dash)
Plots.plot!(create_destroy(sim_single.buyers[1], 2), seriestype = :vline, color = :green, linestyle = :dash)

## EXP = 0, WOM = 0

Random.seed!(20)

sim_single = TO_GO(50, 2, 2, 1, [0.4, 0.4], [1.0, 1.0], "random", 0.0, 0.0, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.5], 4, "dist", TriangularDist(0,1,0.5))

p3 = Plots.plot(sim_single.sellers[1].quality_history, label = "Średnia jakość dobra", xlabel = "Czas", ylabel = "Jakość", legend = :outerright)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 1), label = "Przekonanie o jakości " * L"K_{ijt}")
Plots.plot!(create_buying(sim_single.buyers[1], 1) .+ 0.1, seriestype = :vline, linestyle = :dash, label = "Własne doświadczenie")
Plots.plot!(create_buying(sim_single.buyers[2], 1) .- 0.1, seriestype = :vline, linestyle = :dash, label = "Wymiana informacji")


Plots.savefig(p3, pwd() * "\\plots_to_export\\verification_lambda0theta0.pdf")

Random.seed!(21)

sim_single = TO_GO(50, 2, 2, 0, [0.4, 0.4], [1.0, 1.0], "random", 1.0, 1.0, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.5], 4, "dist", TriangularDist(0,1,0.5))

p4 = Plots.plot(sim_single.sellers[1].quality_history, label = "Średnia jakość dobra", xlabel = "Czas", ylabel = "Jakość", legend = :outerright)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 1), label = "Przekonanie o jakości " * L"K_{ijt}")
Plots.plot!(create_buying(sim_single.buyers[1], 1) .+ 0.1, seriestype = :vline, linestyle = :dash, label = "Własne doświadczenie")

Plots.savefig(p4, pwd() * "\\plots_to_export\\verification_m0.pdf")

########################################### Reselling probability changes due to secondary market

Random.seed!(8)

sim_single = TO_GO(50, 2, 200, 400, [0.4, 0.4], [1.0, 1.0], "random", 0.0, 1.0, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.5], 4, "dist", TriangularDist(0,1,0.5))

function check_changes(x)
    dx = sign.(diff(x))
    lowered = any(dx .< 0) ? true : false
    uppered = any(dx .> 0) ? true : false
    if lowered & uppered
        println("Check this buyer")
    end
end

for i in 1:200
    println(i)
    check_changes(getfield.(sim_single.buyers, :reselling_probability_history)[i])
end

k=193
p2 = Plots.plot(getfield.(sim_single.buyers, :reselling_probability_history)[k], color = :black, linewidth = 2, label = L"p_{it}^R", xlabel = "Czas", ylabel = "Prawdopodobieństwo odsprzedaży", legend = :outerright)
Plots.plot!(create_failed_resell(sim_single.buyers[k], 1), linetype = :vline, color = "red", linestyle = :dash, label = "")
Plots.plot!(create_failed_resell(sim_single.buyers[k], 2), linetype = :vline, color = "red", linestyle = :dash, label = "Nieudana transakcja odsprzedaży")
Plots.plot!(create_resell(sim_single.buyers[k], 1), linetype = :vline, color = "green", linestyle = :dash, label = "Udana transakcja odsprzedaży")
Plots.plot!(create_resell(sim_single.buyers[k], 2), linetype = :vline, color = "green", linestyle = :dash, label = "Udana transakcja odsprzedaży")

Plots.savefig(p2, pwd() * "\\plots_to_export\\verification_pir.pdf")

## No secondary market -> flat probability

sim_single = TO_GO(50, 2, 200, 400, [0.4, 0.4], [1.0, 1.0], "random", 0.0, 1.0, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, false, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.5], 4, "dist", TriangularDist(0,1,0.5))

Plots.plot(getfield.(sim_single.buyers, :reselling_probability_history)[1])
Plots.plot!(create_buying(sim_single.buyers[1], 1), linetype = :vline, color = "grey")
Plots.plot!(create_buying(sim_single.buyers[1], 2), linetype = :vline, color = "grey")
Plots.plot!(create_resell(sim_single.buyers[1], 1), linetype = :vline, color = "red")
Plots.plot!(create_resell(sim_single.buyers[1], 2), linetype = :vline, color = "red")

########################################### Max breakage time

sim_single = TO_GO(1000, 2, 1000, 2000, [0.4, 0.4], [1.0, 1.0], "random", 0.0, 1.0, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, false, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.5], 4, "dist", TriangularDist(0,1,0.5))

history_of_trading = getfield.(sim_single.buyers, :unit_buying_selling_history)

max_h = 0

for b in history_of_trading
    for item in eachindex(b)
        if item > 1
            if b[item - 1].decision == "buy, primary market"
                if b[item].decision == "destroyed"
                    if b[item].t - b[item-1].t > max_h
                        max_h = b[item].t - b[item-1].t
                    end
                end
            end
        end
    end
end

max_h

###########################################

########################################### Meso level #############################################

########################################### No links
## -> no signals shall be received

sim_single = TO_GO(200, 2, 400, 0, [0.4, 0.4], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.1], 8, "dist", TriangularDist(0,1,0.5))

maximum(getfield.(sim_single.buyers, :signal_volume))

########################################### No learning
## -> expectation shall be flat

sim_single = TO_GO(200, 2, 400, 600, [0.4, 0.4], [1.0, 1.0], "random", 0., 0., "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.1], 8, "dist", TriangularDist(0,1,0.5))

allequal([getindex.(x,1) for x in getfield.(sim_single.buyers, :quality_expectation_history)])
allequal([getindex.(x,2) for x in getfield.(sim_single.buyers, :quality_expectation_history)])

########################################### Sampled product quality
## -> in the mean shall be 'close' to producers' decisions

sim_single = TO_GO(200, 2, 1000, 2000, [0.4, 0.4], [1.0, 1.0], "random", 0., 0., "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.1], 4, "dist", TriangularDist(0,1,0.5), 0.1)

plot_quantity(sim_single.sellers, 1)

p = Plots.plot(sim_single.sellers[1].quality_history, label = L"K_{jt}", color="black", linewidth = 2, legend = :outerright, xlabel = "Czas", ylabel = "Jakość")
Plots.plot!(sim_single.sellers[1].quality_history .+ 0.1, color="black", linewidth = 0.5, label = "")
Plots.plot!(sim_single.sellers[1].quality_history .- 0.1, color="black", linewidth = 0.5, label = L"K_{jt} \pm 0.1")
qs = quality_sampled_products(sim_single.buyers, 1, 1)
Plots.plot!(mean_nothing.(qs),color = "blue", label = L"E(\hat K_{ijt})", linewidth = 2)
for t in 1:200
    if length(qs[t]) > 0
        if t == 5
            Plots.plot!([t,t],[mean(qs[t])-std(qs[t]), mean(qs[t])+std(qs[t])],color="blue", label = L"E(\hat K_{ijt}) \pm S(\hat K_{ijt})", linealpha = 0.5)
        else
            Plots.plot!([t,t],[mean(qs[t])-std(qs[t]), mean(qs[t])+std(qs[t])],color="blue", label = "", linealpha = 0.5)
        end
    end
end
Plots.plot!(maximum_nothing.(qs), linewidth = 1.5, color = "blue", label = L"max(\hat K_{ijt})", linestyle = :dash)
Plots.plot!(minimum_nothing.(qs), linewidth = 1.5, color = "blue", label = L"min(\hat K_{ijt})", linestyle = :dash)
p

Plots.savefig(p, pwd() * "\\plots_to_export\\verification_hatK.pdf")

########################################### System level ###########################################

########################################### No secondary market
## -> empty reselling

sim_single = TO_GO(200, 2, 400, 600, [0.4, 0.4], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[1., 2.], [1., 2.]], 0.50, false, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.1], 8, "dist", TriangularDist(0,1,0.5))

Plots.plot(getfield.(sim_single.sellers, :reselling_history))

########################################### Secondary market
## -> market effectiveness?

sim_single = TO_GO(200, 2, 400, 600, [0.4, 0.4], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[1., 2.], [1., 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.1], 8, "dist", TriangularDist(0,1,0.5))

Plots.plot(getfield.(sim_single.sellers, :reselling_history))

########################################### Initial margin 2.0
## -> no initial trading, sellers to decrease

sim_single = TO_GO(200, 2, 400, 600, [1.0, 1.0], [2.0, 2.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[1., 2.], [1., 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.1], 8, "dist", TriangularDist(0,1,0.5))

Plots.plot(getfield.(sim_single.sellers, :quantity_sold_history))
Plots.plot(getfield.(sim_single.sellers, :margin_history))

########################################### Cost coefficient == 1
## -> no trading?

sim_single = TO_GO(200, 2, 400, 600, [1.0, 1.0], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[1., 2.], [1., 2.]], 0.50, false, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.1], 8, "dist", TriangularDist(0,1,0.5))

Plots.plot(getfield.(sim_single.sellers, :quantity_sold_history))
Plots.plot(getfield.(sim_single.sellers, :selling_income_history))

ke = [getindex.(x,1) for x in getfield.(sim_single.buyers, :quality_expectation_history)]
de = [getindex.(x,1) for x in getfield.(sim_single.buyers, :durability_expectation_history)]
srp = getfield.(sim_single.buyers, :std_reservation_price)
fd = getfield.(sim_single.buyers, :future_discount)

u(β,ρ,k,d;H) = β .* k .* (1 .- (ρ.*d).^H) ./ (1 .- (ρ.*d))

ue = u.(srp, fd, ke, de; H = 8)

pr = calculate_price_history(sim_single.sellers[1]; product_life = 8)
pr = fill(pr, 400)

_surplus = [x .> 0 for x in ue .- pr]

any(_surplus[1])

srp[5]
fd[5]

getindex.(calculate_price_history.(sim_single.sellers; product_life= 8), 68)

sim_single.sellers[1].cost_coefficient
sim_single.sellers[1].margin_history

dec = sim_single.buyers[5].unit_buying_selling_history[1]
dec.surplus .== dec.utility .- dec.prices
sim_single.buyers[5].quality_expectation_history[4]
Plots.plot(ke[5])
Plots.plot!(sim_single.sellers[1].quality_history)
Plots.plot(sim_single.sellers[1].margin_history)

utility_expectation_buyers = [[] for _ = 1:200] # to T

for _buyer in sim_single.buyers
    for item in _buyer.unit_buying_selling_history
        if item.decision == "buy, primary market"
            if item.product == 1
                push!(utility_expectation_buyers[item.t], item.utility[1])
            end
        end
    end
end

mean_nothing(x) = length(x) == 0 ? missing : mean(x)



kd(k,d;H) = k .* (1 .- d .^ H) ./ (1 .- d)

Plots.plot(kd.(sim_single.sellers[1].quality_history, sim_single.sellers[1].durability_history; H = 8))
Plots.plot!(mean_nothing.(utility_expectation_buyers))

########################################### Producer profit expectations
##

Random.seed!(12345)

sim_single = TO_GO(200, 2, 400, 600, [1.0, 1.0], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[1., 2.], [1., 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.1], 8, "dist", TriangularDist(0,1,0.5))

Plots.plot(calculate_profit_history(sim_single.sellers[1]))
Plots.plot!(vcat(0,getindex.(sim_single.profit_expected[(getindex.(sim_single.profit_expected, 1) .== 1) .& (getindex.(sim_single.profit_expected, 2) .== "p")], 3)))

Plots.plot(calculate_profit_history(sim_single.sellers[2]))
Plots.plot!(vcat(0, getindex.(sim_single.profit_expected[(getindex.(sim_single.profit_expected, 1) .== 2) .& (getindex.(sim_single.profit_expected, 2) .== "p")], 3)))

########################################### Producer profit expectations
##

Random.seed!(12345)

sim_single = TO_GO(1000, 2, 400, 600, [1.0, 1.0], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[1., 2.], [1., 2.]], 0.50, true, 0, [0.7, 1.], "softmax", ["internal knowledge", "market research"], [0, 0.1], 7, "dist", TriangularDist(0,1,0.5))

max_periods = product_age(sim_single.buyers)

using Latexify

df = DataFrame(Wartość_oczekiwana = round.(mean.(max_periods), digits = 2), 
        Mediana = round.(median.(max_periods), digits = 2),
        Odchylenie_standardowe = round.(std.(max_periods), digits = 2),
        Rozstęp_międzykwartylowy = round.(iqr.(max_periods), digits = 2))

latexify(df, env = :table) |> print

StatsPlots.density(max_periods[1])
StatsPlots.density!(max_periods[2])
StatsPlots.density!(max_periods[3])

Plots.plot(calculate_profit_history(sim_single.sellers[1]))
Plots.plot!(vcat(0,getindex.(sim_single.profit_expected[(getindex.(sim_single.profit_expected, 1) .== 1) .& (getindex.(sim_single.profit_expected, 2) .== "p")], 3)))

Plots.plot(calculate_profit_history(sim_single.sellers[2]))
Plots.plot!(vcat(0, getindex.(sim_single.profit_expected[(getindex.(sim_single.profit_expected, 1) .== 2) .& (getindex.(sim_single.profit_expected, 2) .== "p")], 3)))

countmap(getfield.(vcat(getfield.(sim_single.buyers, :unit_buying_selling_history)...), :decision))


########################################### Producer profit expectations
##

Random.seed!(12345)

sim_single = TO_GO(1000, 2, 400, 600, [0.4, 0.4], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.10, true, 0, [0.7, 1.], "softmax", [true, true], [0.1, 0.1], 10, "dist", TriangularDist(0,1,0.5))

Plots.plot(sum([[any(xi .> 0) for xi in x] for x in getfield.(sim_single.buyers, :unit_possessed_history)]))

Plots.plot([[any(xi .> 0) for xi in x] for x in getfield.(sim_single.buyers, :unit_possessed_history)], color = "grey", linealpha = 0.1)
Plots.plot!(mean([[any(xi .> 0) for xi in x] for x in getfield.(sim_single.buyers, :unit_possessed_history)]), color = "red")

plot_quantity(sim_single.sellers,1)
plot_quantity(sim_single.sellers,2)

Plots.plot(getfield.(sim_single.sellers, :quality_history))
Plots.plot(getfield.(sim_single.sellers, :durability_history))
Plots.plot(getfield.(sim_single.sellers, :margin_history))

max_periods = product_age(sim_single.buyers)

using Latexify

df = DataFrame(Wartość_oczekiwana = round.(mean.(max_periods), digits = 2), 
        Mediana = round.(median.(max_periods), digits = 2),
        Odchylenie_standardowe = round.(std.(max_periods), digits = 2),
        Rozstęp_międzykwartylowy = round.(iqr.(max_periods), digits = 2))

latexify(df, env = :table) |> print

StatsPlots.density(max_periods[1])
StatsPlots.density!(max_periods[2])
StatsPlots.density!(max_periods[3])

Plots.plot(calculate_profit_history(sim_single.sellers[1]))
Plots.plot!(vcat(0,getindex.(sim_single.profit_expected[(getindex.(sim_single.profit_expected, 1) .== 1) .& (getindex.(sim_single.profit_expected, 2) .== "p")], 3)))

Plots.plot(calculate_profit_history(sim_single.sellers[2]))
Plots.plot!(vcat(0, getindex.(sim_single.profit_expected[(getindex.(sim_single.profit_expected, 1) .== 2) .& (getindex.(sim_single.profit_expected, 2) .== "p")], 3)))

countmap(getfield.(vcat(getfield.(sim_single.buyers, :unit_buying_selling_history)...), :decision))


######

"maxIter", "num_sellers", "num_buyers", "num_links", "c", "m", "network_type", "λ_ind", "λ_wom", "buyer_behaviour", "kr", "dr", "mr", "μ_c", "secondary_market_exists", "rand_period", "future_discount_range", "method_weight", "consumer_research", "sample_size", "product_life", "samples_fixed", "do_posterior"

maxIter
num_sellers
num_buyers
num_links
network_type
λ
θ
k_range
d_range
m_range
ρ_range
product_life
secondary_market_exists