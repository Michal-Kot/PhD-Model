
function calculate_state_profit(K::Float64, D::Float64, M::Float64, Q::Float64, o_K::Float64, o_D::Float64, o_P::Float64, N::Int64, μ_c::Float64, cc::Float64, return_type::String)

    """

    Funkcja licząca oczekiwany zysk z danego stanu. Wykorzystywana przez firmę do szacowania efektu zmiany stanu.

    """

    s = LinRange(0,1,N) # standard reservation price, założone dla N klientów
    o_U = s .* sum_of_geom_series_infinite(o_K, o_D) .- o_P # użyteczność dobra konkurencji, jeśli liczba konkurentów > 1, to o_k, o_D i o_P są średnimi
    U = s .* sum_of_geom_series_infinite(K, D)  .- cost_coefficient(K, D, cc) .* sum_of_geom_series_infinite(K, D) .* M # użyteczność mojego dobra przy parametrach K, D, M
    demand = sum((U .> 0) .& (U .> o_U) .& (rand(N) .< (1-D))) # szacowany popyt. warunek 1: moja użyteczność > 0, warunek 2: moja użyteczność wyższa niż użyteczność dobra konkurencyjnego, warunek 3: oczekiwana liczba klientów poszukujących dobra - skalowanie dla dóbr trwałych > 1 okres

    margin_amount = cost_coefficient(K, D, cc) * sum_of_geom_series_infinite(K, D) * (M - 1) # marża na 1 sprzedanym produkcie

    profit = demand * margin_amount + max(0, Q - demand) * (-μ_c) * cost_coefficient(K, D, cc) * sum_of_geom_series_infinite(K, D) # oczekiwany zysk firmy

    if return_type == "profit"
        return profit
    elseif return_type == "demand"
        return demand
    end

end

k = 0.5
d = 2/3
m = 1.1
q = 10.0
N = 100

δ_k = 0.05 * rand()
durability_in_days = 1 / (1 - d)
d_up = durability_in_days / (durability_in_days + 1)
d_down = (durability_in_days - 2) / (durability_in_days - 1)
δ_m = 0.05 * rand()
δ_q = 1

quality_range = [0.2, 0.8]
durability_range = [0.5, 0.95]
margin_range = [0.8, 2.0]

K_range = [max(quality_range[1], k - δ_k), k, min(quality_range[2], k + δ_k)]

D_range = [max(durability_range[1], d_down), d, min(durability_range[2], d_up)]

M_range = [max(margin_range[1], m - δ_m), m, min(margin_range[2], m + δ_m)]

Q_range = [max(0, q - δ_q), q, min(N, q + δ_q)]

o_k = 0.5
o_d = 2/3
o_p = 0.5 * o_k / (1 - o_d) * 1.1

μ_c = 0.25

cc = 0.5

expected_profit_around = [calculate_state_profit(k, d, m, q, o_k, o_d, o_p, N, μ_c, cc, "profit") for k in K_range, d in D_range, m in M_range, q in Q_range]