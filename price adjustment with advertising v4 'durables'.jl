using Distributions, StatsBase
using OrderedCollections
using Plots
using StatsPlots
using LinearAlgebra
using InvertedIndices
using Graphs
using GraphPlot, Compose, Cairo, Fontconfig
using HypothesisTests
using SmoothingSplines
using GLM
using DataFrames


################################# STRUCTS ####################################################################

"""

Sellers form a supply side of the model

Sellers' main goal is to maximize profit π
π = pq - c(q)p
p - price of product, q - quantity of sales, c - function of cost, linear on q, a - advertising intensity

"""

mutable struct seller
    id::Int64
    quality::Float64
    durability::Int64
    average_cost::Float64
    margin::Float64
    quantity_history::Vector{Float64}
    margin_history::Vector{Float64}
    price_return::Bool
    advertising::Float64
    advertising_history::Vector{Float64}
    risk_decisions::Float64
    persuasiveness::Float64
end

function create_sellers(num_sellers::Int64,q::Vector{Float64},c::Vector{Float64},m::Vector{Float64},a::Vector{Float64},r::Vector{Float64},p::Vector{Float64},d::Vector{Int64})::Vector{seller}
    @assert length(q) == num_sellers
    @assert length(c) == num_sellers
    @assert length(m) == num_sellers
    @assert length(a) == num_sellers
    @assert length(r) == num_sellers
    @assert length(p) == num_sellers
    sellers_vector = []
    for s in 1:num_sellers
        new_seller = seller(s, q[s], d[s], c[s], m[s], [], [], false, a[s], [], r[s], p[s])
        push!(sellers_vector, new_seller)
    end
    return sellers_vector
end

"""

Buyers form a demand side of the model

Buyers main goal is to maximize utility of consumption υ(q) = s*eq - p, subject to budget constraint ρ*eq >= p
s - seeking for quality parameter, eq = expected quality of the product, p - price of product, ρ - willingness to pay for unit of quality

"""

mutable struct buyer
    id::Int64
    neighbours::Vector{Int64} 
    std_reservation_price::Float64
    quality_seeking::Float64
    quality_expectation::Vector{Float64}
    quality_expectation_history::Vector{Vector{Float64}}
    durability_expectation::Vector{Float64}
    durability_expectation_history::Vector{Vector{Float64}}
    unit_bought::Vector{Bool}
    unit_bought_history::Vector{Vector{Bool}}
    unit_bought_time::Int64
    quality_of_unit_bought::Float64
    quality_of_unit_bought_history::Vector{Vector{Float64}}
    ad_received::Vector{Bool}
    ad_received_history::Vector{Vector{Bool}}
    broken_product::Bool
    surplus::Vector{Float64}
end

function create_buyers(num_buyers::Int64, num_sellers::Int64, network::SimpleGraph{Int64}, initial_quality_expectation::Vector{Float64}=fill(1.0, num_sellers))::Vector{buyer}
    buyers_vector = []
    for b in 1:num_buyers
        my_neighbours = neighbors(network, b)
        new_buyer = buyer(b, my_neighbours, 2*rand(), rand(Uniform(0,1)), initial_quality_expectation, [], fill(false, num_sellers), [], 0, 0.0, [], [], [], true, [])
        push!(buyers_vector, new_buyer)
    end
    return buyers_vector
end

"""

Create network of a given type - vertices are buyers

if random type, then use Erdos Renyi algorithm, specify number of vertices and edges
if preferential attachment type , then use Albert Barabasi algorithm, specify number of vertices and number of links for a single vertex (constant)

"""

function create_network(type::String; num_buyers::Int64, num_links::Int64=0, pref_attachment_links::Int64=0)

    if type == "random"
        return Graphs.erdos_renyi(num_buyers, num_links)
    elseif type == "preferential_attachment"
        return Graphs.barabasi_albert(num_buyers, pref_attachment_links)
    end
end

#################################### AUX FUNCTIONS ##############################################################

function in_boundaries(x::Float64,lb::Float64,ub::Float64)::Float64
    return max(min(x,ub),lb)
end

function calculate_price(_seller::seller)::Float64
    price = _seller.quality * _seller.average_cost + _seller.margin
    return price
end

function calculate_price_history(_seller::seller)::Vector{Float64}
    price = _seller.quality * _seller.average_cost .+ _seller.margin_history
    return price
end

function u2w(u::Vector{Float64}, p_min::Float64 = 0.1)::Vector{Float64}
    @assert p_min < 1.0
    ut = u .- minimum(u)
    wgt = ut ./ sum(ut) .* (1 - p_min)
    wgt[argmin(u)] = p_min
    return wgt
end

function calculate_profit_history(_seller::seller, γ::Float64, num_buyers::Int64)::Vector{Float64}
    _seller.quantity_history .* calculate_price_history(_seller) .- _seller.average_cost .* _seller.quality .* _seller.quantity_history .- γ * num_buyers * _seller.advertising_history
end

function plot_phasediagram(_seller::seller, function_args::Vector)
    x = calculate_price_history(_seller)
    y = calculate_profit_history(_seller, function_args.γ, function_args.num_buyers)
    d_x = diff(x)
    push!(d_x, 0)
    d_y = diff(y)
    push!(d_y,0)
    labs = string.(1:length(x))
    points_labels = ([mod(x,100) == 0 ? labs[x] : "" for x in 1:length(labs)])
    points_labels[1] = "1"
    quiver(x,y,quiver=(d_x, d_y), color = :green)
    scatter!(x,y, xlabel = "Price", ylabel = "Profit", markershape = :none, markercolor = :white, markerstrokecolor = :white, series_annotations = points_labels, markersize = 0)
end

function create_bool_purchase(n::Int64,k::Int64)::Vector{Bool}
    x = fill(false, n)
    x[k] = true
    return x
end

function calculate_total_surplus(sim_res, return_type::String = "total")

    producer_surplus = sum(sum(calculate_profit_history.(sim_res.sellers, sim_res.function_args.γ, sim_res.function_args.num_buyers)))
    consumer_surplus = sum(sum(getfield.(sim_res.buyers, :surplus)))
    total_surplus = producer_surplus + consumer_surplus

    if return_type == "total"

        return total_surplus

    elseif return_type == "producer"

        return producer_surplus

    elseif return_type == "consumer"

        return consumer_surplus

    end

end

function plot_quantity(sim_res)
    groupedbar(reduce(hcat,[getfield.(sim_res.sellers, :quantity_history)[x] for x in 1:length(sim_res.sellers)]),linecolor=nothing, bar_position = :stack, xlabel = "Time", ylabel = "Quantity of sales", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_quality_expectation(sim_res)
    plot([getindex.(mean(getfield.(sim_res.buyers, :quality_expectation_history)),x) for x in 1:sim_res.function_args.num_sellers], label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_profit_history(sim_res)
    plot(calculate_profit_history.(sim_res.sellers, sim_res.function_args.γ, sim_res.function_args.num_buyers))
end

function plot_advertising(sim_res)
    plot(getfield.(sim_res.sellers, :advertising_history), xlabel = "Time", ylabel = "Advertising intensity [share of population targeted]", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_margin(sim_res)
    Plots.plot([getfield.(sim_res.sellers, :margin_history)[x] for x in 1:length(sim_res.sellers)], xlabel = "Time", ylabel = "Margin", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end


############################################ MAIN SIM FUNCTION #####################################################

function TO_GO(num_sellers::Int64, num_buyers::Int64, max_iter::Int64, λ_ind::Float64, λ_wom::Float64, consumer_behaviour::String; q::Vector{Float64} = fill(1.0, num_sellers), c::Vector{Float64} = fill(0.8, num_sellers), m::Vector{Float64} = fill(0.2, num_sellers), a::Vector{Float64} = fill(0.05, num_sellers), r::Vector{Float64} = fill(0.05, num_sellers), ϵ::Vector{Float64} = fill(1/3, num_sellers), q_init::Vector{Float64} = fill(1.0, num_sellers), p::Vector{Float64} = fill(0.05, num_sellers), num_links::Int64 = 200, δ::Float64 = 0.05, γ::Float64 = 0.050, α::Float64 = 0.005, variant_advertising::Bool = true, allow_negative_margin::Bool = true)

    function_args = (num_sellers = num_sellers, num_buyers = num_buyers, max_iter = max_iter, λ_ind = λ_ind, λ_wom = λ_wom, q = q, c = c, m = m, a = a, r = r, ϵ = ϵ, q_init = q_init, num_links = num_links, δ = δ, γ = γ, α=α)

    # create sellers and buyers

    sellers = create_sellers(num_sellers, q, c, m, a, r, p)

    buyers_network = create_network("random", num_buyers = num_buyers, num_links = num_links)
    buyers = create_buyers(num_buyers, num_sellers, buyers_network, q_init)

    # utlity function definition

    u(k,p;β) = (0.5 + β) .* k .- p

    # TEMP 🍎

    utility_history = []

    # sellers are choosing new price

    #margin_lb = -2
    #margin_ub = 2

    for iter in 1:max_iter

        if iter == 1

            # no change to price in iter 1

            for _seller in sellers

                push!(_seller.margin_history, _seller.margin)
                push!(_seller.advertising_history, _seller.advertising)

            end

        elseif iter == 2

            # random change to price in iter 2
            # price changes by δ
            # advertising changes by α

            for _seller in sellers

                change_margin = sample([-1,0,1], Weights(fill(1/3,3)))
                new_margin = _seller.margin + change_margin * rand() * δ
                _seller.margin = new_margin
                push!(_seller.margin_history, _seller.margin)

                if variant_advertising
                    change_advertising = sample([0,1], Weights(fill(1/2,2)))
                    new_advertising = _seller.advertising + change_advertising * rand() * α
                    new_advertising = in_boundaries(new_advertising, 0.0, 0.10)
                    _seller.advertising = new_advertising
                end

                push!(_seller.advertising_history, _seller.advertising)

            end

        else

            # from iter 3 to the end, firm decides about the price, by comparing price change and resulting profit change

            for _seller in sellers

                # margin has to be in tha range between cost of production & 2
                #margin_ub = 2 - _seller.quality * _seller.average_cost

                # how did margin change?
                margin_change = sign(_seller.margin_history[end] - _seller.margin_history[end-1])

                # how did advertising change?

                advertising_change = sign(_seller.advertising_history[end] - _seller.advertising_history[end-1])

                # how did profit change?

                profit_history = calculate_profit_history(_seller, γ, num_buyers)
                profit_change = profit_history[end] - profit_history[end-1]
                
                if profit_change > 0

                    if rand() <= _seller.risk_decisions
                        new_margin = _seller.margin + margin_change * rand() * δ
                    else
                        new_margin = _seller.margin
                    end

                    if rand() <= _seller.risk_decisions
                        new_advertising = _seller.advertising + advertising_change * rand() * α
                    else
                        new_advertising = _seller.advertising
                    end

                elseif profit_change < 0

                    if rand() <= _seller.risk_decisions
                        new_margin = _seller.margin_history[end-1] + (-1) * margin_change * rand() * δ
                    else
                        new_margin = _seller.margin_history[end-1]
                    end

                    if rand() <= _seller.risk_decisions
                        new_advertising = _seller.advertising_history[end-1] + (-1) * advertising_change * rand() * α
                    else
                        new_advertising = _seller.advertising_history[end-1]
                    end

                elseif profit_change == 0

                    if rand() <= _seller.risk_decisions
                        margin_change = sample([-1,0,1], Weights(fill(1/3,3)))
                        new_margin = _seller.margin + margin_change * rand() * δ
                    else
                        new_margin = _seller.margin
                    end

                    if rand() <= _seller.risk_decisions
                        advertising_change = sample([-1,0,1], Weights(fill(1/3,3)))
                        new_advertising = _seller.advertising + advertising_change * rand() * α
                    else
                        new_advertising = _seller.advertising
                    end

                end

                if iter >= 6

                    if all((_seller.margin_history[(end-4):end] .* _seller.quantity_history[(end-4):end]) .<= 0) #& (_seller.margin > 0)
                        new_margin = _seller.margin - rand() * δ
                    end

                    if (all(calculate_profit_history(_seller, γ, num_buyers)[(end-4):end] .< 0)) & (_seller.margin < 0)
                        new_margin = _seller.margin + δ
                    end

                end

                if !allow_negative_margin
                    new_margin = in_boundaries(new_margin, 0, 2)
                end

                _seller.margin = new_margin
                push!(_seller.margin_history, _seller.margin)

                if variant_advertising
                    new_advertising = in_boundaries(new_advertising, 0.0, 0.10)
                    _seller.advertising = new_advertising
                end

                push!(_seller.advertising_history, _seller.advertising)
            
            end

        end

        #### RYNEK ####

        demand = Int64[]

        # REKLAMY

        advertising_intensity = getfield.(sellers, :advertising)
        prices = calculate_price.(sellers)

        for _buyer in buyers

            # Advertising

            if rand() < sum(advertising_intensity)
                firm_ad_seen = sample(1:num_sellers, Weights(advertising_intensity))
                advertising_bool = fill(false, num_sellers)
                advertising_bool[firm_ad_seen] = true
                _buyer.ad_received = advertising_bool
                push!(_buyer.ad_received_history, _buyer.ad_received)
            else
                _buyer.ad_received = fill(false, num_sellers)
                push!(_buyer.ad_received_history, _buyer.ad_received)
            end

            # Product usage -> if broken

            if iter >= 2

                product_durability = sellers[_buyer.unit_bought].durability
                product_amortization = iter - _buyer.unit_bought_time
                probability_to_break = 0.5 / product_durability * product_amortization

                if rand() < probability_to_break
                    _buyer.broken_product = true
                end

            end



            wtp_advertising = 1 .+ _buyer.ad_received .* p

            available_products = prices .<= _buyer.std_reservation_price .* _buyer.quality_expectation .* wtp_advertising

            if any(available_products)

                utility = u.(_buyer.quality_expectation, prices; β = _buyer.quality_seeking)

                push!(utility_history, (id=_buyer.id, qs=_buyer.quality_seeking, qe=_buyer.quality_expectation, p=prices, ap=available_products, u=utility))

                if consumer_behaviour == "deterministic"

                    # consumer choses product that maximizes utility

                    chosen_product = argmax(utility .- 1e10 .* .!available_products)

                elseif consumer_behaviour == "stochastic"

                    # consumer randomly choses products, accordingly to utility (may choose not the best one)

                    weight = u2w(utility) .* available_products

                    chosen_product = sample(1:num_sellers, Weights(weight))

                end
                
                _buyer.unit_bought = create_bool_purchase(num_sellers, chosen_product)

                push!(_buyer.unit_bought_history, _buyer.unit_bought)

                _buyer.quality_of_unit_bought = rand(Uniform(sellers[chosen_product].quality - ϵ[chosen_product], sellers[chosen_product].quality + ϵ[chosen_product]))
                
                push!(demand, chosen_product)

                push!(_buyer.quality_of_unit_bought_history, _buyer.unit_bought .*  _buyer.quality_of_unit_bought)

                buyer_surplus = _buyer.quality_expectation[chosen_product] * _buyer.std_reservation_price * wtp_advertising[chosen_product] - prices[chosen_product]

                push!(_buyer.surplus, buyer_surplus)
            
            else

                _buyer.unit_bought = fill(false, length(sellers))
                push!(_buyer.unit_bought_history, _buyer.unit_bought)
                _buyer.quality_of_unit_bought = 0
                push!(_buyer.surplus, 0.0)

            end

        end

        
        # kupujący weryfikują oczekiwaną jakość produktów

        for _buyer in buyers

            if length(_buyer.neighbours) > 0
                average_quality_neighbours = reduce(+,getfield.(buyers[_buyer.neighbours], :quality_of_unit_bought) .* getfield.(buyers[_buyer.neighbours], :unit_bought)) ./ reduce(+, getfield.(buyers[_buyer.neighbours], :unit_bought))
                bought_neighbours = .!isnan.(average_quality_neighbours)
                average_quality_neighbours[.!bought_neighbours] .= 0
            else
                bought_neighbours = fill(0, num_sellers)
                average_quality_neighbours = fill(0, num_sellers)
            end

            new_quality_expectation = _buyer.quality_expectation .+ λ_ind .* _buyer.unit_bought .* (_buyer.quality_of_unit_bought .- _buyer.quality_expectation) .+ λ_wom .* bought_neighbours .* (average_quality_neighbours .- _buyer.quality_expectation)

            # forgetting
            # to be decided if to use it 🍓

            #new_quality_expectation = new_quality_expectation .+ (_buyer.quality_expectation .>= 1.0) .* λ_fgt .* (1 .- _buyer.quality_expectation)

            _buyer.quality_expectation = new_quality_expectation

            push!(_buyer.quality_expectation_history, new_quality_expectation)

        end

        # firma weryfikuje Q sprzedaży

        for _seller in sellers

            push!(_seller.quantity_history, count(demand .== _seller.id))
        end

    end

    return (sellers = sellers, buyers = buyers, function_args = function_args, ut_his = utility_history)

end

__precompile__()

# DEAD MARKET, no comms, no p2p comms

# Cel: pokazać zbieżność proponowanego modelu z modelem Izquierdo, który w identycznej sytuacji prowadzi do tej samej konkluzji. Brak komunikacji pomiędzy agentami prowadzi, w długim terminie do zaniku rynku, czyli braku popytu na produkty.

# Przykład bez zapominania, rynek przestaje istnieć.

# Adjustment ceny jest zbyt wolny, żeby przeciwstawić się negatywnemu wpływowi braku komunikacji p2p.

sim_with_obs_11 = TO_GO(4, 500, 10000, 1.0, 0.0, "deterministic"; q = [1.0, 1.0, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], num_links = 0, variant_advertising = true, allow_negative_margin = true, δ = 0.005)

plot_margin(sim_with_obs_11)
plot_advertising(sim_with_obs_11)
plot_quantity(sim_with_obs_11)
plot_quality_expectation(sim_with_obs_11)

####################### EXPERIMENT 1, EQUAL QUALITIES ############################################

# Cel 1: Weryfikacja hipotezy, że rynki gdzie występuje zróżnicowanie jakości produktów są bardziej efektywne niż pozostałe. 
# Cel 2: Sprawdzić, czy komunikacja reklamowa zwiększa efektywność rynku w przypadku, jeśli jakość produktów jest "równa" (EX1) oraz jakość produktów jest różna - istnieje firma dominująca pozostałe pod względem jakości oferowanego produktu (EX2)

# Efekt 1: różnicowanie jakości produktu powoduje, że rynek jest bardziej efektywny - suma nadwyżek producentów i konsumentów jest wyższa
# Efekt 2: w przypadku rynku z różnicowaniem produktu, komunikacja reklamowa pozwala istotnie zwiększyć efektywność, podczas gdy na rynku z jednakową jakością produktu efekt nie jest obserwowany

ex1_total_surplus_eq = []
ex1_producer_surplus_eq = []
ex1_consumer_surplus_eq = []
ex1_price_eq = []
ex1_quantity_eq = []
ex1_advertising_eq = []

ex2_total_surplus_dq = []
ex2_producer_surplus_dq = []
ex2_consumer_surplus_dq = []
ex2_price_dq = []
ex2_quantity_dq = []
ex2_advertising_dq = []
ex2_advertising_highest = []

for i in 1:100

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    sim_with_ads = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.0, 1.0, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], num_links = 1000)
    push!(ex1_total_surplus_eq, calculate_total_surplus(sim_with_ads, "total"))
    push!(ex1_producer_surplus_eq, calculate_total_surplus(sim_with_ads, "producer"))
    push!(ex1_consumer_surplus_eq, calculate_total_surplus(sim_with_ads, "consumer"))
    push!(ex1_price_eq, mean(calculate_price_history.(sim_with_ads.sellers)))
    push!(ex1_quantity_eq, mean(getfield.(sim_with_ads.sellers, :quantity_history)))
    push!(ex1_advertising_eq, mean(getfield.(sim_with_ads.sellers, :advertising_history)))

    sim_with_ads = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.6, 1.3, 1.0, 0.7], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], num_links = 1000, q_init = [1.6, 1.3, 1.0, 0.7])
    push!(ex2_total_surplus_dq, calculate_total_surplus(sim_with_ads, "total"))
    push!(ex2_producer_surplus_dq, calculate_total_surplus(sim_with_ads, "producer"))
    push!(ex2_consumer_surplus_dq, calculate_total_surplus(sim_with_ads, "consumer"))
    push!(ex2_price_dq, mean(calculate_price_history.(sim_with_ads.sellers)))
    push!(ex2_quantity_dq, mean(getfield.(sim_with_ads.sellers, :quantity_history)))
    push!(ex2_advertising_dq, mean(getfield.(sim_with_ads.sellers, :advertising_history)))

end

ex1_p1 = plot(sort(ex1_total_surplus_eq), (1:length(ex1_total_surplus_eq))./length(ex1_total_surplus_eq), xlabel = "Total Market surplus", ylabel = "Probability", title = "Empirical Cumluative Distribution", label = "equal quality", legend=:bottomright)
plot!(sort(ex2_total_surplus_dq), (1:length(ex2_total_surplus_dq))./length(ex2_total_surplus_dq), xlabel = "Total Market surplus", ylabel = "Probability", title = "Empirical Cumluative Distribution", label = "non-equal quality")

ex1_p2 = plot(sort(ex1_producer_surplus_eq), (1:length(ex1_producer_surplus_eq))./length(ex1_producer_surplus_eq), xlabel = "Producer surplus", ylabel = "Probability", title = "Empirical Cumluative Distribution", label = "equal quality", legend=:bottomright)
plot!(sort(ex2_producer_surplus_dq), (1:length(ex2_producer_surplus_dq))./length(ex2_producer_surplus_dq), xlabel = "Producer surplus", ylabel = "Probability", title = "Empirical Cumluative Distribution", label = "non-equal quality")

ex1_p3 = plot(sort(ex1_consumer_surplus_eq), (1:length(ex1_consumer_surplus_eq))./length(ex1_consumer_surplus_eq), xlabel = "Consumer surplus", ylabel = "Probability", title = "Empirical Cumluative Distribution", label = "equal quality", legend=:bottomright)
plot!(sort(ex2_consumer_surplus_dq), (1:length(ex2_consumer_surplus_dq))./length(ex2_consumer_surplus_dq), xlabel = "Consumer surplus", ylabel = "Probability", title = "Empirical Cumluative Distribution", label = "non-equal quality")

plot(sort(mean.(ex1_price_eq)), (1:length(ex1_price_eq))./length(ex1_price_eq), xlabel = "Price", ylabel = "Probability", title = "Empirical Cumluative Distribution", label = "equal quality", legend=:bottomright)
plot!(sort(mean.(ex2_price_dq)), (1:length(ex2_price_dq))./length(ex2_price_dq), xlabel = "Price", ylabel = "Probability", title = "Empirical Cumluative Distribution", label = "non-equal quality", legend=:bottomright)

plot(sort(mean.(ex1_advertising_eq)), (1:length(ex1_advertising_eq))./length(ex1_advertising_eq), xlabel = "Advertising", ylabel = "Probability", title = "Empirical Cumluative Distribution", label = "equal quality", legend=:bottomright)
plot!(sort(mean.(ex2_advertising_dq)), (1:length(ex2_advertising_dq))./length(ex2_advertising_dq), xlabel = "Advertising", ylabel = "Probability", title = "Empirical Cumluative Distribution", label = "non-equal quality", legend=:bottomright)

p = plot(mean(ex1_price_eq), 
    xlabel = "T", ylabel = "Price", 
    title = "Average price", label = "equal quality", legend=:topright)

plot!(mean(ex2_price_dq), 
    xlabel = "T", ylabel = "Price", 
    title = "Average price", label = "not equal quality", legend=:topright)

EqualVarianceTTest(mean(ex1_price_eq), mean(ex2_price_dq))

p = plot(mean(ex1_quantity_eq), 
    xlabel = "T", ylabel = "Quantity", 
    title = "Average quantity", label = "equal quality", legend=:bottomright)
    
plot!(mean(ex2_quantity_dq), 
    xlabel = "T", ylabel = "Quantity", 
    title = "Average quantity", label = "non-equal quality")

EqualVarianceTTest(mean(ex1_quantity_eq), mean(ex2_quantity_dq))

p = plot(mean(ex1_advertising_eq), 
    xlabel = "T", ylabel = "Advertising", 
    title = "Average advertising", label = "equal quality", legend=:bottomright)
    
plot!(mean(ex2_advertising_dq), 
    xlabel = "T", ylabel = "Advertising", 
    title = "Average advertising", label = "non-equal quality")

EqualVarianceTTest(mean(ex1_advertising_eq), mean(ex2_advertising_dq))

# Cel: weryfikacja jak persuasiveness reklamy wpływa na wybory konsumentów

persuasiveness = LinRange(0.05:0.05:0.25)
max_iter = 50
persuasiveness_sensitivity_eq = []
persuasiveness_sensitivity_dq = []

for p in persuasiveness
    for iter in 1:max_iter
        println((p, iter))
        sim_with_obs_13 = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.0, 1.0, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.6, 0.6, 0.6, 0.6], ϵ = [0.33, 0.33, 0.33, 0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], p = [p, p, p, p], q_init = [1.0, 1.0, 1.0, 1.0], num_links = 1000)
        push!(persuasiveness_sensitivity_eq, (price = calculate_price_history.(sim_with_obs_13.sellers), quantity = getfield.(sim_with_obs_13.sellers, :quantity_history), profit = calculate_profit_history.(sim_with_obs_13.sellers, sim_with_obs_13.function_args.γ, sim_with_obs_13.function_args.num_buyers), advertising = getfield.(sim_with_obs_13.sellers, :advertising_history), surplus = calculate_total_surplus(sim_with_obs_13)))
        sim_with_obs_13 = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.6, 1.3, 1.0, 0.7], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], p = [p, p, p, p], q_init = [1.6, 1.3, 1.0, 0.7], num_links = 1000)
        push!(persuasiveness_sensitivity_dq, (price = calculate_price_history.(sim_with_obs_13.sellers), quantity = getfield.(sim_with_obs_13.sellers, :quantity_history), profit = calculate_profit_history.(sim_with_obs_13.sellers, sim_with_obs_13.function_args.γ, sim_with_obs_13.function_args.num_buyers), advertising = getfield.(sim_with_obs_13.sellers, :advertising_history), surplus = calculate_total_surplus(sim_with_obs_13)))
    end
end

persuasiveness_unique = repeat(collect(persuasiveness), inner = max_iter)

plot(collect(persuasiveness), [mean(mean(mean(getfield.(persuasiveness_sensitivity_eq, :price)[persuasiveness_unique .∈ upu]))) for upu in unique(persuasiveness_unique)], label = "eq", ylabel = "price")
plot!(collect(persuasiveness), [mean(mean(mean(getfield.(persuasiveness_sensitivity_dq, :price)[persuasiveness_unique .∈ upu]))) for upu in unique(persuasiveness_unique)])

plot(collect(persuasiveness), [mean(mean(mean(getfield.(persuasiveness_sensitivity_eq, :quantity)[persuasiveness_unique .∈ upu]))) for upu in unique(persuasiveness_unique)], label = "eq", ylabel = "price")
plot!(collect(persuasiveness), [mean(mean(mean(getfield.(persuasiveness_sensitivity_dq, :quantity)[persuasiveness_unique .∈ upu]))) for upu in unique(persuasiveness_unique)])

plot(collect(persuasiveness), [mean(mean(mean(getfield.(persuasiveness_sensitivity_eq, :surplus)[persuasiveness_unique .∈ upu]))) for upu in unique(persuasiveness_unique)], label = "eq", ylabel = "price")
plot!(collect(persuasiveness), [mean(mean(mean(getfield.(persuasiveness_sensitivity_dq, :surplus)[persuasiveness_unique .∈ upu]))) for upu in unique(persuasiveness_unique)])


# Cel: czy punkt startowy m = 0.2 ma znaczenie? Czy wynik symulacji jest niezależny od punktu startowego?

# Przykład dla równych Q

init_margin_sensitivity = []
tested_margins = LinRange(0:0.1:1.0)

for m_init in tested_margins
    println(m_init)
    for iter in 1:10
        sim_with_obs_12 = TO_GO(4, 250, 500, 0.25, 0.25, "deterministic"; q = [1.0, 1.0, 1.0, 1.0], m = [m_init, m_init, m_init, m_init], c = [0.7,0.7,0.7,0.7], ϵ = [0.33,0.33,0.33,0.33], a = [0.01, 0.01, 0.01, 0.01], r = [0.4, 0.4, 0.4, 0.4], num_links = 1000)
        push!(init_margin_sensitivity, sim_with_obs_12)
    end
end

init_margins = getindex.(getfield.(getfield.(init_margin_sensitivity, :function_args), :m), 1)
price_trajectories = calculate_price_history.(getindex.(getfield.(init_margin_sensitivity, :sellers), 1))
quantity_trajectories = getfield.(getindex.(getfield.(init_margin_sensitivity, :sellers), 1), :quantity_history)
profit_trajectories = calculate_profit_history.(getindex.(getfield.(init_margin_sensitivity, :sellers), 1), 0.05, 250)

plot([mean(price_trajectories[init_margins .∈ uim]) for uim in unique(init_margins)], labels = reshape("Margin =  " .* string.(collect(tested_margins)),1,length(tested_margins)), xlabel = "Time", ylabel = "Price")
plot([mean(quantity_trajectories[init_margins .∈ uim]) for uim in unique(init_margins)], labels = reshape("Margin =  " .* string.(collect(tested_margins)),1,length(tested_margins)), xlabel = "Time", ylabel = "Quantity")
plot([mean(profit_trajectories[init_margins .∈ uim]) for uim in unique(init_margins)], labels = reshape("Margin =  " .* string.(collect(tested_margins)),1,length(tested_margins)), xlabel = "Time", ylabel = "Profit")

# Przykład dla różnych Q

init_margin_sensitivity = []
tested_margins = LinRange(0:0.1:1.0)

for m_init in tested_margins
    println(m_init)
    for iter in 1:10
        sim_with_obs_12 = TO_GO(4, 250, 500, 0.25, 0.25, "deterministic"; q = [1.6, 1.2, 1.0, 1.0], m = [m_init, m_init, m_init, m_init], c = [0.7,0.7,0.7,0.7], ϵ = [0.33,0.33,0.33,0.33], a = [0.01, 0.01, 0.01, 0.01], r = [0.4, 0.4, 0.4, 0.4], num_links = 1000)
        push!(init_margin_sensitivity, sim_with_obs_12)
    end
end

init_margins = getindex.(getfield.(getfield.(init_margin_sensitivity, :function_args), :m), 1)
price_trajectories = calculate_price_history.(getindex.(getfield.(init_margin_sensitivity, :sellers), 1))
quantity_trajectories = getfield.(getindex.(getfield.(init_margin_sensitivity, :sellers), 1), :quantity_history)
profit_trajectories = calculate_profit_history.(getindex.(getfield.(init_margin_sensitivity, :sellers), 1), 0.05, 250)

plot([mean(price_trajectories[init_margins .∈ uim]) for uim in unique(init_margins)], labels = reshape("Margin =  " .* string.(collect(tested_margins)),1,length(tested_margins)), xlabel = "Time", ylabel = "Price")
plot([mean(quantity_trajectories[init_margins .∈ uim]) for uim in unique(init_margins)], labels = reshape("Margin =  " .* string.(collect(tested_margins)),1,length(tested_margins)), xlabel = "Time", ylabel = "Quantity")
plot([mean(profit_trajectories[init_margins .∈ uim]) for uim in unique(init_margins)], labels = reshape("Margin =  " .* string.(collect(tested_margins)),1,length(tested_margins)), xlabel = "Time", ylabel = "Profit")

# Cel: pokazać jak wygląda total surplus w zależności od zróżnicowania jakości na rynku

q_diff_sensitivity = []
tested_diff = LinRange(0:0.1:0.5)
maxiter = 25

for q_diff in tested_diff
    println(q_diff)
    for iter in 1:maxiter
        if q_diff == 0
            q_sampled = fill(1.0,4)
        else
            q_sampled = rand(Uniform(1 - q_diff, 1 + q_diff),4)
        end
        sim_with_obs_12 = TO_GO(4, 250, 500, 0.25, 0.25, "deterministic"; q = q_sampled, m = [0.2, 0.2, 0.2, 0.2], c = [0.7,0.7,0.7,0.7], ϵ = [0.33,0.33,0.33,0.33], a = [0.01, 0.01, 0.01, 0.01], r = [0.4, 0.4, 0.4, 0.4], q_init = q_sampled, num_links = 1000)
        push!(q_diff_sensitivity, sim_with_obs_12)
    end
end

init_diff = repeat(tested_diff, inner = maxiter)
price_trajectories = calculate_price_history.(getindex.(getfield.(q_diff_sensitivity, :sellers), 1))
quantity_trajectories = getfield.(getindex.(getfield.(q_diff_sensitivity, :sellers), 1), :quantity_history)
profit_trajectories = calculate_profit_history.(getindex.(getfield.(q_diff_sensitivity, :sellers), 1), 0.05, 250)

plot([mean(price_trajectories[init_diff .∈ udf])[2:end] for udf in unique(init_diff)], labels = reshape("Q diff =  " .* string.(collect(tested_diff)),1,length(tested_diff)), xlabel = "Time", ylabel = "Price")
plot([mean(quantity_trajectories[init_diff .∈ udf]) for udf in unique(init_diff)][2:end], labels = reshape("Q diff =  " .* string.(collect(tested_diff)),1,length(tested_diff)), xlabel = "Time", ylabel = "Quantity")
plot([sum(profit_trajectories[init_diff .∈ udf]) for udf in unique(init_diff)][2:end], labels = reshape("Q diff =  " .* string.(collect(tested_diff)),1,length(tested_diff)), xlabel = "Time", ylabel = "Profit")

####

quality_diff_on_advertising = []
quality_best = LinRange(1.0:0.10:1.6)
maxiter = 100

for qb in quality_best
    for iter in 1:maxiter
        println((qb,iter))
        sim_with_obs_12 = TO_GO(4, 500, 500, 0.25, 0.25, "deterministic"; q = [qb, 1.0, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], q_init = [qb, 1.0, 1.0, 1.0], num_links = 1000)
        push!(quality_diff_on_advertising, (price = calculate_price_history.(sim_with_obs_12.sellers), quantity = getfield.(sim_with_obs_12.sellers, :quantity_history), profit = calculate_profit_history.(sim_with_obs_12.sellers, sim_with_obs_12.function_args.γ, sim_with_obs_12.function_args.num_buyers), advertising = getfield.(sim_with_obs_12.sellers, :advertising_history)))
    end
end



quality_best_unique = repeat(collect(quality_best), inner = maxiter)
plot([getindex(mean(getfield.(quality_diff_on_advertising, :price)[quality_best_unique .∈ uqbu]),1) for uqbu in unique(quality_best_unique)])
plot([getindex(mean(getfield.(quality_diff_on_advertising, :advertising)[quality_best_unique .∈ uqbu]),1) for uqbu in unique(quality_best_unique)])
plot([getindex(mean(getfield.(quality_diff_on_advertising, :quantity)[quality_best_unique .∈ uqbu]),1)[2:end] for uqbu in unique(quality_best_unique)])
plot([getindex(mean(getfield.(quality_diff_on_advertising, :profit)[quality_best_unique .∈ uqbu]),1)[2:end] for uqbu in unique(quality_best_unique)])

# Przykład bez zapominania, zostaje na rynku jeden producent

sim_with_obs_12 = TO_GO(4, 500, 5000, 0.25, 0.25, "deterministic"; q = [1.6, 1.3, 1.0, 0.7], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], q_init = [1.6, 1.3, 1.0, 0.7], num_links = 1000)

plot_margin(sim_with_obs_12)
plot_advertising(sim_with_obs_12)
plot_quantity(sim_with_obs_12)
plot_quality_expectation(sim_with_obs_12)
plot_profit_history(sim_with_obs_12)
calculate_total_surplus(sim_with_obs_12)

sim_with_obs_13 = TO_GO(4, 500, 10000, 1.0, 0.0, "deterministic"; q = [1.6, 1.3, 1.0, 0.7], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0],r = [0.4, 0.4, 0.4, 0.4], q_init = [1.6, 1.3, 1.0, 0.7], num_links = 0)

plot_margin(sim_with_obs_13)
plot_quantity(sim_with_obs_13)
plot_quality_expectation(sim_with_obs_13)

sim_with_obs_14 = TO_GO(4, 500, 10000, 1.0, 0.0, 0.0, 0.01; q = [1.6, 1.3, 1.0, 0.7], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0],r = [0.4, 0.4, 0.4, 0.4], q_init = [1.6, 1.3, 1.0, 0.7], num_links = 0)

plot_margin(sim_with_obs_14)
plot_quantity(sim_with_obs_14)
plot_quality_expectation(sim_with_obs_14)


######### NO MARKET FOR THE HIGHEST QUALITY, no comms, p2p comms

# Cel: pokazać nieefektywność rynku, pozbawionego komunikacji reklamowej. Nawet jeśli konsumenci mogą komunikować się ze sobą, to jeśli na rynku występuje heterogeniczność produktów i wertykalne zróżnicowanie (występują produkty wysokiej i niskiej jakości) to brak komunikacji reklamowej może prowadzić do braku popytu na produkty wysokiej jakości. Nieefektywność wynika z faktu, że konsumenci byliby skłonni kupować produkt wysokiej jakości, ale nie potrafią prawidłowo oszacować jego jakości.

# Przykład bez zapominania, z warunkiem początkowym e(q) = 1.0. Producent o najwyższej jakości nigdy nie wchodzi na rynek.

sim_with_obs_21 = TO_GO(4, 500, 1000, 0.10, 0.10, 0.10, 0.0; q = [1.6, 1.2, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.7,0.7,0.7,0.7], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], num_links = 1000)

plot_phasediagram(sim_with_obs_21.sellers[1], sim_with_obs_21.function_args)
plot_margin(sim_with_obs_21)
plot_advertising(sim_with_obs_21)
plot_quantity(sim_with_obs_21)
plot_quality_expectation(sim_with_obs_21)
plot_profit_history(sim_with_obs_21)




calculate_total_surplus(sim_with_obs_21)

# Przykład bez zapominania, z e(q) = q̂. Producent o najwyższej jakości istnieje

sim_with_obs_22 = TO_GO(4, 500, 1000, 0.25, 0.25, 0.0, 0.0; q = [1.6, 1.2, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.7,0.7,0.7,0.7], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], q_init = [1.6, 1.2, 1.0, 1.0], num_links = 1000)

plot_margin(sim_with_obs_22)
plot_quantity(sim_with_obs_22)
plot_quality_expectation(sim_with_obs_22)

# Przykład z zapominaniem i e(q) = q̂. Producent o najwyższej jakości pojawia się i znika

sim_with_obs_23 = TO_GO(4, 500, 2000, 0.25, 0.25, 0.0, 0.04; q = [1.35, 1.2, 1.0, 0.8], m = [0.2, 0.2, 0.2, 0.2], c = [0.7,0.7,0.7,0.7], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], q_init = [1.35, 1.2, 1.0, 0.8], num_links = 1000)

plot_margin(sim_with_obs_23)
plot_quantity(sim_with_obs_23)
plot_quality_expectation(sim_with_obs_23)

plot(getindex.(mean(getfield.(sim_with_obs_23.buyers, :quality_expectation_history)),1), ylim = (0,1.5))
plot!(twinx(), getindex.(reduce(+,getfield.(sim_with_obs_23.buyers, :unit_bought_history)),1), ylim=(0,300))

######### MARKET FOR THE HIGHEST QUALITY, comms, p2p comms

# Cel: pokazać, że występowanie na rynku komunikacji reklamowej, pozwala producentowi dobra o wysokiej jakości wygenerować popyt na swój produkt.

# Przykład bez zapominania, wstępnie e(q) = 1.0. Producent dobra o najwyższej jakości istnieje na rynku

sim_with_obs_31 = TO_GO(4, 500, 1000, 0.25, 0.25, 0.25, 0.; q = [1.6, 1.2, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.7,0.7,0.7,0.7], ϵ = [0.33,0.33,0.33,0.33], a = [0.01, 0.01, 0.01, 0.01], r = [0.4, 0.4, 0.4, 0.4], num_links = 1000)

plot_margin(sim_with_obs_31)
plot_quantity(sim_with_obs_31)
plot_quality_expectation(sim_with_obs_31)

calculate_total_surplus(sim_with_obs_31)

plot(getindex.(mean(getfield.(sim_with_obs_31.buyers, :quality_expectation_history)),1), ylim = (1,1.5),color = "black", legend = nothing)
plot!(twinx(),getindex.(reduce(+,getfield.(sim_with_obs_31.buyers, :ad_received_history)),1), color = RGBA(1,0,0,0.5), legend = nothing)
plot!(twinx(), getindex.(reduce(+,getfield.(sim_with_obs_31.buyers, :unit_bought_history)),1), ylim=(0,100), legend=nothing, color = :green)

# Przykład z zapominaniem, wstępnie e(q) = 1.0. Producent dobra o najwyższej jakości istnieje na rynku

sim_with_obs_32 = TO_GO(4, 500, 500, 0.25, 0.25, 0.25, 0.025; q = [1.5, 1.2, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.7,0.6,0.5,0.4], ϵ = [0.33,0.33,0.33,0.33], a = [0.10, 0.01, 0.01, 0.01], r = [0.4, 0.4, 0.4, 0.4], num_links = 1000)

Plots.plot([getfield.(sim_with_obs_32.sellers, :margin_history)[x] for x in 1:length(sim_with_obs_32.sellers)])

plot_phasediagram(sim_with_obs_32.sellers[1], sim_with_obs_32.function_args)

sum(getfield(sim_with_obs_32.sellers[1], :quantity_history)) 

groupedbar(reduce(hcat,[getfield.(sim_with_obs_32.sellers, :quantity_history)[x] for x in 1:length(sim_with_obs_32.sellers)]),linecolor=nothing, bar_position = :stack)

# Przykład z zapominaniem

sim_with_obs_33 = TO_GO(4, 500, 500, 0.25, 0.25, 0.25, 0.025; q = [1.5, 1.2, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.7,0.6,0.5,0.4], ϵ = [0.33,0.33,0.33,0.33], a = [0.10, 0.01, 0.01, 0.01], r = [0.4, 0.4, 0.4, 0.4], q_init = [1.5, 1.2, 1.0, 1.0], num_links = 1000)

Plots.plot([getfield.(sim_with_obs_33.sellers, :margin_history)[x] for x in 1:length(sim_with_obs_33.sellers)])

plot_phasediagram(sim_with_obs_33.sellers[1], sim_with_obs_33.function_args)

sum(getfield(sim_with_obs_33.sellers[1], :quantity_history)) 

groupedbar(reduce(hcat,[getfield.(sim_with_obs_33.sellers, :quantity_history)[x] for x in 1:length(sim_with_obs_33.sellers)]),linecolor=nothing, bar_position = :stack)



####################### EXPERIMENT 2 ############################################

num_links_sensitivity = []

for nl in LinRange(0:100:1000)
    println(nl)
    for i in 1:50
        nl_res = TO_GO(4, 100, 250, 0.25, 0.25, 0.25, 0.0; q = [1.60, 1.30, 0.90, 0.70], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a=[0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], num_links = nl) 
        push!(num_links_sensitivity, nl_res)
        nl_res = TO_GO(4, 100, 250, 0.25, 0.25, 0.25, 0.0; q = [1., 1., 1., 1.], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a=[0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], num_links = nl) 
        push!(num_links_sensitivity, nl_res)
    end
end

avg_quantity = [mean(mean(getfield.(x, :quantity_history))) for x in getfield.(num_links_sensitivity, :sellers)]

num_links = getfield.(getfield.(num_links_sensitivity, :function_args), :num_links)

quality_diff = getindex.(getfield.(getfield.(num_links_sensitivity, :function_args), :q), 1) .== 1.6

plot(sort(unique(num_links)), map(x -> mean(avg_quantity[(num_links .== x) .& (quality_diff .== 1)]), sort(unique(num_links))), xlabel = "Number of links", ylabel = "Average quantity", label = "Market w/ vertical diff")
plot!(sort(unique(num_links)), map(x -> mean(avg_quantity[(num_links .== x) .& (quality_diff .== 0)]), sort(unique(num_links))), label = "Market w/o vertical diff")

EqualVarianceTTest(avg_quantity[quality_diff .== 1], avg_quantity[quality_diff .== 0])

####### EXPERIMENT 3 ####################

# Cel: w jaki sposób zróżnicowanie jakości producentów wpływa na elastyczność cenową?

using GLM, DataFrames

sim_with_obs_41 = TO_GO(4, 500, 1000, 0.25, 0.25, 0.25, 0.0; q = [1.0, 1.0, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.01, 0.01, 0.01, 0.01], r = [0.4, 0.4, 0.4, 0.4], num_links = 1000)

sim_with_obs_42 = TO_GO(4, 500, 1000, 0.25, 0.25, 0.25, 0.0; q = [1.4, 1.2, 1.0, 0.8], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.01, 0.01, 0.01, 0.01], r = [0.4, 0.4, 0.4, 0.4], num_links = 1000)

function plot_demand_curve(sellers)
    p = plot()
    colors = palette(:tab10)
    market_average_price = mean(calculate_price_history.(sellers))
    k=1
    for _seller in sellers
        prices = calculate_price_history(_seller)
        relative_prices = prices ./ market_average_price
        unique_relative_prices = sort(unique(relative_prices))
        unique_quantities = map(x -> mean(_seller.quantity_history[relative_prices .== x]), unique_relative_prices)
        is_95 = (unique_relative_prices .>= percentile(unique_relative_prices, 5)) .& (unique_relative_prices .<= percentile(unique_relative_prices, 95))
        df = DataFrame(x=unique_relative_prices, y = unique_quantities)
        df_mod = df[is_95,:]
        demand_curve_model = GLM.lm(@formula(y~x), df_mod)
        demand_curve = predict(demand_curve_model, df)
        p = scatter!(unique_relative_prices, unique_quantities, xlabel = "Price", ylabel = "Quantity", legend = nothing, color=colors[k])
        p = plot!(df.x,demand_curve, color=colors[k])
        k += 1
    end
    return p
end

plot_demand_curve(sim_with_obs_41.sellers)
plot_demand_curve(sim_with_obs_42.sellers)

####### ELASTICITY OF DEMAND ###############

ex5_eqQ = []
ex5_dfQ = []

for iter in 1:400
    println(iter)

    a = rand() * 0.05

    sim_with_obs_51 = TO_GO(4, 500, 250, 0.25, 0.25, 0.25, 0.0; q = [1.0, 1.0, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.7,0.7,0.7,0.7], ϵ = [0.33,0.33,0.33,0.33], a = [a, a, a, a], r = [0.4, 0.4, 0.4, 0.4], num_links = 1000)

    push!(ex5_eqQ, (sim_with_obs_51.function_args.a, calculate_price_elasticity.(sim_with_obs_51.sellers)[1]))

    sim_with_obs_52 = TO_GO(4, 500, 250, 0.25, 0.25, 0.25, 0.0; q = [1.6, 1.2, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.7,0.7,0.7,0.7], ϵ = [0.33,0.33,0.33,0.33], a = [a, a, a, a], r = [0.4, 0.4, 0.4, 0.4], num_links = 1000)

    push!(ex5_dfQ, (sim_with_obs_52.function_args.a, calculate_price_elasticity.(sim_with_obs_52.sellers)[1]))

end

function calculate_price_elasticity(_seller)
    y = _seller.quantity_history
    x = calculate_price_history(_seller)
    df = DataFrame(log_y=log.(1 .+ y), log_x=log.(1 .+ x))
    model_fit = GLM.lm(@formula(log_y~log_x), df)
    return coef(model_fit)[2]
end

scatter(getindex.(getindex.(ex5_eqQ,1),1), getindex.(ex5_eqQ,2), smooth=true)
scatter!(getindex.(getindex.(ex5_dfQ,1),1), getindex.(ex5_dfQ,2), smooth = true)



##############

price_strategy_qd = []
price_strategy_qe = []

for r in LinRange(0:0.1:0.5)
    for rc in LinRange(0:0.1:0.5)
        println((r,rc))
        for iter in 1:50
            ps_res = TO_GO(4, 250, 250, 0.25, 0.25, 0.25, 0.0; q = [1.40, 1.20, 1.0, 0.80], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [r, rc, rc, rc], num_links = 500) 
            my_r = getindex(ps_res.function_args.r, 1)
            comp_r = getindex(ps_res.function_args.r, 2)
            profit = sum.(calculate_profit_history.(ps_res.sellers, 0.05, 250))[1]
            push!(price_strategy_qd, (my_r, comp_r, profit))
            ps_res = TO_GO(4, 250, 250, 0.25, 0.25, 0.25, 0.0; q = [1.0, 1.0, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [r, rc, rc, rc], num_links = 500) 
            my_r = getindex(ps_res.function_args.r, 1)
            comp_r = getindex(ps_res.function_args.r, 2)
            profit = sum.(calculate_profit_history.(ps_res.sellers,0.05,250))[1]
            push!(price_strategy_qe, (my_r, comp_r, profit))
        end
    end
end

my_r_qd = getindex.(price_strategy_qd, 1)
comp_r_qd = getindex.(price_strategy_qd, 2)
profit_qd = getindex.(price_strategy_qd, 3)
profit_matrix_qd = [mean(profit_qd[(my_r_qd .== r) .& (comp_r_qd .== rc)]) for r in sort(unique(my_r_qd)), rc in sort(unique(comp_r_qd))]
p_qd = heatmap(sort(unique(my_r_qd)), sort(unique(comp_r_qd)), profit_matrix_qd', xlabel = "High-quality producer's agressiveness", ylabel = "Competitors' agressiveness", levels=8,clim = (0,4000))

my_r_qe = getindex.(price_strategy_qe, 1)
comp_r_qe = getindex.(price_strategy_qe, 2)
profit_qe = getindex.(price_strategy_qe, 3)
profit_matrix_qe = [mean(profit[(my_r_qe .== r) .& (comp_r_qe .== rc)]) for r in sort(unique(my_r_qe)), rc in sort(unique(comp_r_qe))]
p_qe = heatmap(sort(unique(my_r_qe)), sort(unique(comp_r_qe)), profit_matrix_qe', xlabel = "High-quality producer's agressiveness", ylabel = "Competitors' agressiveness", levels=8,clim = (0,4000))
#### Experiments plans

# popracować nad ekperymentem z elastyznością cenową popytu
# Cel: quality diff vs. advertising -> profit
# Cel: quality diff -> when does product starts to sell?