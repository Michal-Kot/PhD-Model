

# DEAD MARKET, no comms, no p2p comms

# Cel: pokazać zbieżność proponowanego modelu z modelem Izquierdo, który w identycznej sytuacji prowadzi do tej samej konkluzji. Brak komunikacji pomiędzy agentami prowadzi, w długim terminie do zaniku rynku, czyli braku popytu na produkty.

# Przykład bez zapominania, rynek przestaje istnieć.

# Adjustment ceny jest zbyt wolny, żeby przeciwstawić się negatywnemu wpływowi braku komunikacji p2p.

sim_with_obs_11 = TO_GO(4, 500, 10000, 1.0, 0.0, "deterministic"; q = [1.0, 1.0, 1.0, 1.0], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], num_links = 0, variant_advertising = true, allow_negative_margin = true, δ = 0.005)

plot_margin(sim_with_obs_11)
plot_advertising(sim_with_obs_11)
plot_quantity(sim_with_obs_11)
plot_quality_expectation(sim_with_obs_11)



####################################### OPTIMAL PRODUCT DURATION FOR HIGH QUALITY PRODUCER ######

# Cel: weryfikacja, czy wyższa jakość pozwala producentowi oferować niższe duration, bez spadku profitu

duration_test = collect(1:1:25)
max_iter = 10000
duration_sensitivity_eq = []
duration_sensitivity_dq = []

for iter in 1:max_iter
    println(iter)
    d1 = sample(duration_test)
    d2 = sample(duration_test)
    sim_with_obs_13 = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.0, 1.0], m = [0.2, 0.2], c = [0.6, 0.6], ϵ = [0.33, 0.33], a = [0.0, 0.0], r = [0.4, 0.4], q_init = [1.0, 1.0], d = [d1,d2], d_init = [Float64(d1),Float64(d2)], num_links = 1000)
    push!(duration_sensitivity_eq, (d = getfield(getfield(sim_with_obs_13, :function_args), :d), profit = sum.(calculate_profit_history.(sim_with_obs_13.sellers, sim_with_obs_13.function_args.γ, sim_with_obs_13.function_args.num_buyers))))

    sim_with_obs_13 = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.3, 0.7], m = [0.2, 0.2], c = [0.6, 0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], q_init = [1.3, 0.7], d = [d1,d2], d_init = [Float64(d1),Float64(d2)], num_links = 1000)
    push!(duration_sensitivity_dq, (d = getfield(getfield(sim_with_obs_13, :function_args), :d), profit = sum.(calculate_profit_history.(sim_with_obs_13.sellers, sim_with_obs_13.function_args.γ, sim_with_obs_13.function_args.num_buyers))))
end

d_eq = getfield.(duration_sensitivity_eq, :d)
p_eq = getfield.(duration_sensitivity_eq, :profit)
d_dq = getfield.(duration_sensitivity_dq, :d)
p_dq = getfield.(duration_sensitivity_dq, :profit)

d_eq_firm1 = [mean(getindex.(p_eq, 1)[d_eq .== fill([d1,d2], length(d_eq))]) for d1 in duration_test, d2 in duration_test]
d_eq_firm2 = [mean(getindex.(p_eq, 2)[d_eq .== fill([d1,d2], length(d_eq))]) for d1 in duration_test, d2 in duration_test]
d_dq_firm1 = [mean(getindex.(p_dq, 1)[d_dq .== fill([d1,d2], length(d_dq))]) for d1 in duration_test, d2 in duration_test]
d_dq_firm2 = [mean(getindex.(p_dq, 2)[d_dq .== fill([d1,d2], length(d_dq))]) for d1 in duration_test, d2 in duration_test]

Plots.plot([argmax(c) for c in eachcol(d_dq_firm1)], duration_test, label = "BR firm 1, eq. K", xlabel = "Duration of firm 1", ylabel = "Duration of firm 2", xlim = (0,25), xticks = 0:1:25, yticks = 0:1:25, aspect_ratio = :equal)
Plots.plot!(duration_test, [argmax(c) for c in eachrow(d_dq_firm2)], label = "BR firm 2, eq. K")
Plots.plot!([argmax(c) for c in eachcol(d_eq_firm1)], duration_test, label = "BR firm 1, df. K")
Plots.plot!(duration_test, [argmax(c) for c in eachrow(d_eq_firm2)], label = "BR firm 2, df. K")

sing_min = minimum([minimum(mtr) for mtr in (d_eq_firm1, d_eq_firm2, d_dq_firm1, d_dq_firm2)])
sing_max = maximum([maximum(mtr) for mtr in (d_eq_firm1, d_eq_firm2, d_dq_firm1, d_dq_firm2)])

dbl_min = minimum([minimum(mtr) for mtr in (d_eq_firm1 .+ d_eq_firm2, d_dq_firm1 .+ d_dq_firm2)])
dbl_max = maximum([maximum(mtr) for mtr in (d_eq_firm1 .+ d_eq_firm2, d_dq_firm1 .+ d_dq_firm2)])

heatmap(d_eq_firm1', xlabel = "Firm 1 duration", ylabel = "Firm 2 duration", title = "Profit of firm 1", clim=(sing_min,sing_max))
heatmap(d_eq_firm2', xlabel = "Firm 1 duration", ylabel = "Firm 2 duration", title = "Profit of firm 2", clim=(sing_min,sing_max))
heatmap((d_eq_firm1 .+ d_eq_firm2)', xlabel = "Firm 1 duration", ylabel = "Firm 2 duration", title = "Profit of firm 2",clim=(dbl_min, dbl_max))

heatmap(d_dq_firm1', xlabel = "Firm 1 duration", ylabel = "Firm 2 duration", title = "Profit of firm 1", clim=(sing_min,sing_max))
heatmap(d_dq_firm2', xlabel = "Firm 1 duration", ylabel = "Firm 2 duration", title = "Profit of firm 2", clim=(sing_min,sing_max))
heatmap((d_dq_firm1 .+ d_dq_firm2)', xlabel = "Firm 1 duration", ylabel = "Firm 2 duration", title = "Profit of firm 2",clim=(dbl_min, dbl_max))

#################### DURATION TO QUALITY RELATION

duration_test = collect(1:1:25)
quality_test = collect(0.5:0.1:1.5)
max_iter = 5000
duration_quality_sensitivity_eq = []

for iter in 1:max_iter
    println(iter)
    q1 = sample(quality_test)
    q2 = sample(quality_test)
    d1 = sample(duration_test)
    d2 = sample(duration_test)
    sim_with_obs_13 = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [q1, q2], m = [0.2, 0.2], c = [0.6, 0.6], ϵ = [0.33, 0.33], a = [0.0, 0.0], r = [0.4, 0.4], q_init = [q1, q2], d = [d1,d2], d_init = [Float64(d1),Float64(d2)], num_links = 1000)
    push!(duration_quality_sensitivity_eq, (d = getfield(getfield(sim_with_obs_13, :function_args), :d), q = getfield(getfield(sim_with_obs_13, :function_args), :q), profit = sum.(calculate_profit_history.(sim_with_obs_13.sellers, sim_with_obs_13.function_args.γ, sim_with_obs_13.function_args.num_buyers)), csurplus = calculate_total_surplus(sim_with_obs_13, "consumer")))

end

d = getfield.(duration_quality_sensitivity_eq, :d)
q = getfield.(duration_quality_sensitivity_eq, :q)
p = getfield.(duration_quality_sensitivity_eq, :profit)
c = getfield.(duration_quality_sensitivity_eq, :csurplus)

d_eq_firm1 = [mean(getindex.(p, 1)[(getindex.(d,1) .== d1) .& (getindex.(q,1) .== q1)]) for d1 in duration_test, q1 in 1 .+ collect(-0.5:0.1:0.5)]
heatmap(d_eq_firm1', xlabel = "Durability", ylabel = "Quality", yticks = (1:11, 0.5:0.1:1.5), xticks = (1:25))

d_eq_firm2 = [mean(getindex.(p, 2)[(getindex.(d,2) .== d2) .& (getindex.(q,2) .== q2)]) for d2 in duration_test, q2 in 1 .+ collect(-0.5:0.1:0.5)]
heatmap(d_eq_firm2', xlabel = "Durability", ylabel = "Quality", yticks = (1:11, 0.5:0.1:1.5), xticks = (1:25))

d_eq_firm1 = [mean(getindex.(c, 1)[(getindex.(d,1) .== d1) .& (getindex.(q,1) .== q1)]) for d1 in duration_test, q1 in 1 .+ collect(-0.5:0.1:0.5)]
heatmap(d_eq_firm1', xlabel = "Durability", ylabel = "Quality", yticks = (1:11, 0.5:0.1:1.5), xticks = (1:25))

d_eq_firm3 = [mean(getindex.(p, 1)[(getindex.(q,1) .== q1) .& (getindex.(q,2) .== q2)]) for q1 in 1 .+ collect(-0.5:0.1:0.5), q2 in 1 .+ collect(-0.5:0.1:0.5)]
heatmap(d_eq_firm3', xlabel = "Quality", ylabel = "Quality", yticks = (1:11, 0.5:0.1:1.5), xticks = (1:11, 0.5:0.1:1.5))

d_eq_firm3 = [mean(getindex.(c, 1)[(getindex.(q,1) .== q1) .& (getindex.(q,2) .== q2)]) for q1 in 1 .+ collect(-0.5:0.1:0.5), q2 in 1 .+ collect(-0.5:0.1:0.5)]
heatmap(d_eq_firm3', xlabel = "Quality", ylabel = "Quality", yticks = (1:11, 0.5:0.1:1.5), xticks = (1:11, 0.5:0.1:1.5))

d_eq_firm3 = [mean(sum.(p)[(getindex.(q,1) .== q1) .& (getindex.(q,2) .== q2)]) for q1 in 1 .+ collect(-0.5:0.1:0.5), q2 in 1 .+ collect(-0.5:0.1:0.5)]
heatmap(d_eq_firm3', xlabel = "Quality", ylabel = "Quality", yticks = (1:11, 0.5:0.1:1.5), xticks = (1:11, 0.5:0.1:1.5))

d_eq_firm3 = [mean(sum.(p)[(getindex.(d,1) .== d1) .& (getindex.(d,2) .== d2)]) for d1 in duration_test, d2 in duration_test]
heatmap(d_eq_firm3', xlabel = "Quality", ylabel = "Quality")

d_eq_firm3 = [mean(c[(getindex.(d,1) .== d1) .& (getindex.(d,2) .== d2)]) for d1 in duration_test, d2 in duration_test]
heatmap(d_eq_firm3', xlabel = "Quality", ylabel = "Quality")

d_eq_firm3 = [mean((sum.(p).+c)[(getindex.(d,1) .== d1) .& (getindex.(d,2) .== d2)]) for d1 in duration_test, d2 in duration_test]
heatmap(d_eq_firm3', xlabel = "Quality", ylabel = "Quality")

d_eq_firm33 = [mean((sum.(p) .+ c)[(getindex.(q,1) .== q1) .& (getindex.(q,2) .== q2)]) for q1 in quality_test, q2 in quality_test]
heatmap(d_eq_firm33')

d_eq_firm3 = [mean(c[(getindex.(d,1) .== d1) .& (getindex.(d,2) .== d2) .& (getindex.(q,1) .<= 1.0) .& (getindex.(q,2) .<= 1.0)]) for d1 in duration_test, d2 in duration_test]
heatmap(d_eq_firm3', xlabel = "Quality", ylabel = "Quality")

d_eq_firm3 = [mean(c[(getindex.(d,1) .== d1) .& (getindex.(d,2) .== d2) .& (getindex.(q,1) .>= 1.0) .& (getindex.(q,2) .>= 1.0)]) for d1 in duration_test, d2 in duration_test]
heatmap(d_eq_firm3', xlabel = "Quality", ylabel = "Quality")

d_eq_firm3 = [mean(c[(getindex.(d,1) .== d1) .& (getindex.(d,2) .== d2) .& (getindex.(q,1) .>= 1.0) .& (getindex.(q,2) .<= 1.0)]) for d1 in duration_test, d2 in duration_test]
heatmap(d_eq_firm3', xlabel = "Quality", ylabel = "Quality")

d_eq_firm3 = [mean(c[(getindex.(d,1) .== d1) .& (getindex.(d,2) .== d2) .& (getindex.(q,1) .<= 1.0) .& (getindex.(q,2) .>= 1.0)]) for d1 in duration_test, d2 in duration_test]
heatmap(d_eq_firm3', xlabel = "Quality", ylabel = "Quality")

######################################################



#######################################################

riskiness_test = collect(0:0.05:0.5)
quality_test = collect(0.5:0.1:1.5)
max_iter = 1000
riskiness_persuasiveness_sensitivity_eq = []

for iter in 1:max_iter
    println(iter)
    q1 = sample(quality_test)
    q2 = sample(quality_test)
    r1 = sample(riskiness_test)
    r2 = sample(riskiness_test)
    sim_with_obs_13 = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [q1, q2], m = [0.2, 0.2], c = [0.6, 0.6], ϵ = [0.33, 0.33], a = [0.0, 0.0], r = [r1, r2], q_init = [q1, q2], d = [1,1], d_init = [1., 1.], num_links = 1000)
    push!(riskiness_persuasiveness_sensitivity_eq, (d = getfield(getfield(sim_with_obs_13, :function_args), :r), q = getfield(getfield(sim_with_obs_13, :function_args), :q), profit = sum.(calculate_profit_history.(sim_with_obs_13.sellers, sim_with_obs_13.function_args.γ, sim_with_obs_13.function_args.num_buyers)), csurplus = calculate_total_surplus(sim_with_obs_13, "consumer")))

end

r = getfield.(riskiness_persuasiveness_sensitivity_eq, :d)
q = getfield.(riskiness_persuasiveness_sensitivity_eq, :q)
p = getfield.(riskiness_persuasiveness_sensitivity_eq, :profit)
c = getfield.(riskiness_persuasiveness_sensitivity_eq, :csurplus)

d_eq_firm1 = [mean(getindex.(c, 1)[(getindex.(r,1) .== r1) .& (getindex.(q,1) .== q1)]) for r1 in riskiness_test, q1 in quality_test]
heatmap(d_eq_firm1', xlabel = "Risk seeking", ylabel = "Quality")

d_eq_firm2 = [mean((sum.(p) .+ c)[(getindex.(r,1) .== r1) .& (getindex.(r,2) .== r2)]) for r1 in riskiness_test, r2 in riskiness_test]
heatmap(d_eq_firm2', xlabel = "Risk seeking", ylabel = "Quality")

d_eq_firm2 = [mean(sum.(p)[(getindex.(r,1) .== r1) .& (getindex.(r,2) .== r2) .& (getindex.(q,1) .<= 1.0) .& (getindex.(q,2) .>= 1.0)]) for r1 in riskiness_test, r2 in riskiness_test]
heatmap(d_eq_firm2', xlabel = "Risk seeking", ylabel = "Quality")

using CairoMakie

function plot_contour(xs,ys,zs)
    f = Figure()
    Axis(f[1, 1])
    CairoMakie.contour!(xs, ys, zs)
    f
end

plot_contour(quality_test, riskiness_test, [mean(sum.(p)[(getindex.(q,1) .== q1) .& (getindex.(q,2) .== q2)]) for q1 in quality_test, q2 in quality_test])


################# OPTIMAL DURATION FOR 

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

sim_with_obs_12 = TO_GO(4, 500, 500, 0.25, 0.25, "deterministic"; q = [1.6, 1.3, 1.0, 0.7], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], q_init = [1.6, 1.3, 1.0, 0.7], d = [2, 2, 2, 2], d_init = [2.0,2.0,2.0,2.0], num_links = 1000)

sum(getindex.((getfield.(sim_with_obs_12.ut_his, :wtp) .- getfield.(sim_with_obs_12.ut_his, :p)),argmax.(getfield.(sim_with_obs_12.ut_his,:u))))

Plots.plot(calculate_total_surplus(sim_with_obs_12, "consumer", false))
Plots.plot!(calculate_total_surplus(sim_with_obs_12, "producer", false))

plot_margin(sim_with_obs_12)
plot_advertising(sim_with_obs_12)
plot_quantity(sim_with_obs_12)
plot_quality_expectation(sim_with_obs_12)
plot_profit_history(sim_with_obs_12)
plot_price(sim_with)
calculate_total_surplus(sim_with_obs_12, "consumer")

plot(calculate_price_history.(sim_with_obs_12.sellers))

plot([getindex.(mean(getfield.(sim_with_obs_12.buyers, :durability_expectation_history)), x) for x in 1:4])

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