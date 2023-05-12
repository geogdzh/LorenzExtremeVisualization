module LorenzExtremeVisualization

export RungeKutta4
struct RungeKutta4{S, T, U}
    k⃗::S
    x̃::T
    xⁿ⁺¹::T
    t::U
end
RungeKutta4(n) = RungeKutta4(zeros(n, 4), zeros(n), zeros(n), [0.0])

function (step::RungeKutta4)(f, x, dt)
    @inbounds let
        @. step.x̃ = x
        step.k⃗[:, 1] .= f(step.x̃, step.t[1])
        @. step.x̃ = x + step.k⃗[:, 1] * dt / 2
        @. step.t += dt / 2
        step.k⃗[:, 2] .= f(step.x̃, step.t[1])
        @. step.x̃ = x + step.k⃗[:, 2] * dt / 2
        step.k⃗[:, 3] .= f(step.x̃, step.t[1])
        @. step.x̃ = x + step.k⃗[:, 3] * dt
        @. step.t += dt / 2
        step.k⃗[:, 4] .= f(step.x̃, step.t[1])
        @. step.xⁿ⁺¹ = x + (step.k⃗[:, 1] + 2 * step.k⃗[:, 2] + 2 * step.k⃗[:, 3] + step.k⃗[:, 4]) * dt / 6
    end
    return nothing
end

end # module LorenzExtremeVisualization
