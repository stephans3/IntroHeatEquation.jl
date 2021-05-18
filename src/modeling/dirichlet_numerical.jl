### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 37b018c4-34b0-11eb-234d-4fbeaa683873
using LinearAlgebra, SparseArrays, Plots, DifferentialEquations

# ╔═╡ baed3ace-34ad-11eb-0cf6-19357238b1de
md"## Heat equation with Dirichlet boundary conditions

The domain is given by the one-dimensional rod $\Omega = (0,L)$ and the time by $t \in (0,T).$ The rod has the physical properties
- thermal conductivity $\lambda$,
- specific heat capacity $c$,
- mass density $\rho$.

### Heat Equation

The diffusivity constant is given by $\alpha = \frac{\lambda}{c \rho}$. The heat equation is noted as

$\dot{\vartheta}(t,x) = \alpha \frac{\partial^2}{\partial x^2} \vartheta(t,x)$

for $(t,x) \in (0,T) \times \Omega$ with initial condition

$\vartheta(0,x) = \vartheta_{0}(x)$ 

for $x \in \overline{\Omega}$ and boundary conditions  $\qquad \vartheta(\cdot,0) = \vartheta(\cdot,L) = 0$.
"

# ╔═╡ 8f5d25dc-34ae-11eb-2aa7-39e601f648a1
md"#### Code: implementing physical constants"

# ╔═╡ 8f493cca-34ae-11eb-3e60-31ab68674568
# Length of rod
L = 0.5

# ╔═╡ 8f149f6a-34ae-11eb-01f2-5983e256a7de
# Thermal conductivity
λ = 45.0 

# ╔═╡ 8efa0286-34ae-11eb-1952-e790d141e19b
# Specific heat capacity
c = 480.0

# ╔═╡ 8ecb6bba-34ae-11eb-2c16-819675412725
# Mass density
ρ = 7800.0

# ╔═╡ e6433430-34af-11eb-318d-2199461e5ec9
# Thermal diffusivity
α = λ/(c*ρ)

# ╔═╡ 37c91fae-34b0-11eb-0f82-639173602ac0
md"### Spatial approximation

The one-dimensional rod is discretized as a one-dimensional grid with N points and the finite discretization is noted by

$\Delta x = \frac{L}{N-1}.$

The resulting grid points have the position $x^{0} = 0$, $x^{n} = n~\Delta x$ and $x^{N-1} = L$. The second order derivative is approximated using the [Taylor series](https://en.wikipedia.org/wiki/Taylor_series) as 

$\left. \frac{\partial^2 f(x)}{\partial x^2} \right\rvert_{x = \tilde{x}} \approx \frac{1}{\Delta x^2} \left[ f(\tilde{x} - \Delta x) - 2 f(\tilde{x}) + f(\tilde{x} + \Delta x)  \right].$

This [Finite Difference](https://en.wikipedia.org/wiki/Finite_difference_method) scheme is written in matrix-vector form as 

$\frac{\partial^2 \vartheta(t,x)}{\partial x^2} \approx \frac{1}{\Delta x^2}
\begin{pmatrix}
-2 & 1 & 0 & \cdots &  & 0 \\
1 & -2 & 1 & 0  & \cdots & \vdots \\
0 & 1 & -2 & 1 & \ddots &  \\
 & & \ddots & \ddots & \ddots &  & \\
& & & 1 & -2 & 1 \\
& & &  & 1 & -2 \\
\end{pmatrix}
~
\begin{pmatrix}
\vartheta(t, x^{0}) \\
\vartheta(t, x^{1}) \\
\vdots \\
\vartheta(t, x^{N-1}) \\
\end{pmatrix}.$
"

# ╔═╡ 65049b1a-3549-11eb-33ee-9dfae609164d
# Number of grid elements
N = 101

# ╔═╡ a36fd440-3549-11eb-0e93-0d36244400ef
# Finite discretization
Δx = L/(N-1)   

# ╔═╡ c18a976a-3549-11eb-2991-799f558e4493
# 1-dimensional grid
xspan = 0 : Δx : L

# ╔═╡ d7d5be46-3549-11eb-3a12-afab5ab4bce5
# Diffusion matrix
M = spdiagm(-1 => ones(N-1), 0 => -2*ones(N), 1 => ones(N-1));

# ╔═╡ 045066e2-354a-11eb-164b-09fc7095de47
Matrix(M)[1:5,1:5]

# ╔═╡ 37822ec8-34b0-11eb-175d-7745382c2955
md"### Boundary conditions

The Dirichlet boundary condition states that the value at both sides is fixed as 

$\qquad \vartheta(\cdot,0) = \vartheta(\cdot, x^{0}) = 0$

on the left boundary and 

$\vartheta(\cdot,L) = \vartheta(\cdot, x^{N-1})  = 0$

on the right boundary. This means, there is no dynamical (e.g. diffusive) behaviour at both boundaries and thus the first and the last row of the diffusion matrix M is set to 0. "

# ╔═╡ 3768b7e0-34b0-11eb-0772-f5f4c95cf8f5
# First row
M[1,1] = M[1,2] = 0;

# ╔═╡ 37529500-34b0-11eb-1ee4-b5837bea23c9
# Last row
M[end,end] = M[end,end-1] = 0;

# ╔═╡ 373652f0-34b0-11eb-2fea-4f14b3ee9140
md"### Heat equation as ODE

After the spatial approximation the heat equation has a form of 

$\dot{\theta}(t) = \frac{\alpha}{\Delta x^2} ~ M ~ \theta(t)$

with M as diffusion matrix and $\quad \theta(t) = \left( \vartheta(t, x^{0}), \cdots, \vartheta(t, x^{N-1}) \right)^{\top}$. This is an ordinary differential equation (ODE) that can be solved with common solvers like forward Euler method or Runge-Kutta scheme."

# ╔═╡ 371c6d52-34b0-11eb-12ab-a51e10fe33c5
# Heat Equation as ODE
function heat_eq(dθ, θ, p, t)
    return dθ .= (1/Δx^2) * α * M * θ
end

# ╔═╡ 37013e12-34b0-11eb-3d23-9d43ca3f06a3
md"### Initial Conditions

The initial data of the original heat equation is assumed as $\quad \vartheta_{0}(x) := m \left[ L~x - x^2 \right]~$ with  $m > 0$. This function is approximated to gain the initial conditions of the ODE 

$\theta(0) = \left( \vartheta_{0}(x^{0}), \cdots , \vartheta_{0}(x^{N-1}) \right)^{\top}.$
"

# ╔═╡ 78b7fa20-354c-11eb-3754-d5d515d0a2d5
#Initial conditions: θ(0) = m * (L * x - x^2)
function initial_data(x, params)
    # α = params[1] # unused here 
    L = params[2]
    m = params[3]

    return m*( L * x - x^2)
end

# ╔═╡ 77906a42-354c-11eb-397e-9dca847ab964
# Initial heat distribution
θ₀ = zeros(N);

# ╔═╡ 93696cf4-354d-11eb-0af3-9b27a8186fbd
# Amplification of initial data
@bind amp html"<input type='range' min='0' max='100' step='1'>"

# ╔═╡ a91c8b80-354d-11eb-03ec-8f4386a69708
m = amp

# ╔═╡ 86f0f21c-354d-11eb-2a32-299183f9118b
# Parameters
param = [α, L, amp]

# ╔═╡ 666a8878-354d-11eb-22dc-e94a72ede22d
let

	for i = 1 : length(xspan)
		θ₀[i] = initial_data(xspan[i], param)
	end
	
	plot(xspan, θ₀, xlabel="Position x", ylabel="Temperature", legend=false)
end

# ╔═╡ cfc0e592-354d-11eb-2b06-458d198c38e5
md"### Numerical stability

Before the approximated heat equation can be simulated, the [numerical stability](https://en.wikipedia.org/wiki/Numerical_stability) has to be proven with the [von Neumann stability analysis](https://en.wikipedia.org/wiki/Von_Neumann_stability_analysis). 

In case of the forward Euler method one finds the approach

$\vartheta(t^{k+1}, x^{n}) = \vartheta(t^{k}, x^{n}) + \alpha \frac{\Delta t}{\Delta x^2} \left( \vartheta(t^{k}, x^{n-1}) - 2 ~ \vartheta(t^{k}, x^{n}) + \vartheta(t^{k}, x^{n+1}) \right)$

or equivalent

$\vartheta(t^{k+1}, x^{n}) = \left(1 - 2 \alpha \frac{\Delta t}{\Delta x^2} \right) \vartheta(t^{k}, x^{n}) + \alpha \frac{\Delta t}{\Delta x^2}  \left( \vartheta(t^{k}, x^{n-1}) + \vartheta(t^{k}, x^{n+1}) \right).$

The sampling time $\Delta t$ has to be chosen to guarantee 

$\left(1 - 2 \alpha \frac{\Delta t}{\Delta x^2} \right) > 0.$

Therefore, the upper limit of the sampling time is given by 

$\Delta t < \frac{1}{2} \frac{\Delta x^2}{\alpha}.$
"

# ╔═╡ 25789968-3552-11eb-1e25-35c0d1b1c461
# Upper limit of sampling time
ul = 0.5 * Δx^2 / α

# ╔═╡ 56daee6e-3554-11eb-37cd-49ad4ac56805
# Sampling time
Δt = 0.8

# ╔═╡ 56bb225a-3554-11eb-0062-6d03b8ea772f
# time discretization
tspan = (0.0, 1000.0)

# ╔═╡ 569a20f0-3554-11eb-0479-f1e583620682
md"### Code: Solving the ODE"

# ╔═╡ 56790ed8-3554-11eb-3c51-b9b31a4c340e
# Solving the ODE
let
	prob = ODEProblem( heat_eq, θ₀, tspan ) # ODE Problem
	
	sol = solve(prob,Euler(),dt=Δt,progress=true, save_everystep=false, 
				save_start=true) # Solving the ODE
	
	plot(xspan, sol.u[2], xlabel = "Position x", ylabel="Temperature", legend=false)
end

# ╔═╡ Cell order:
# ╟─baed3ace-34ad-11eb-0cf6-19357238b1de
# ╟─8f5d25dc-34ae-11eb-2aa7-39e601f648a1
# ╠═8f493cca-34ae-11eb-3e60-31ab68674568
# ╠═8f149f6a-34ae-11eb-01f2-5983e256a7de
# ╠═8efa0286-34ae-11eb-1952-e790d141e19b
# ╠═8ecb6bba-34ae-11eb-2c16-819675412725
# ╠═e6433430-34af-11eb-318d-2199461e5ec9
# ╟─37c91fae-34b0-11eb-0f82-639173602ac0
# ╠═37b018c4-34b0-11eb-234d-4fbeaa683873
# ╠═65049b1a-3549-11eb-33ee-9dfae609164d
# ╠═a36fd440-3549-11eb-0e93-0d36244400ef
# ╠═c18a976a-3549-11eb-2991-799f558e4493
# ╠═d7d5be46-3549-11eb-3a12-afab5ab4bce5
# ╠═045066e2-354a-11eb-164b-09fc7095de47
# ╟─37822ec8-34b0-11eb-175d-7745382c2955
# ╠═3768b7e0-34b0-11eb-0772-f5f4c95cf8f5
# ╠═37529500-34b0-11eb-1ee4-b5837bea23c9
# ╟─373652f0-34b0-11eb-2fea-4f14b3ee9140
# ╠═371c6d52-34b0-11eb-12ab-a51e10fe33c5
# ╟─37013e12-34b0-11eb-3d23-9d43ca3f06a3
# ╠═78b7fa20-354c-11eb-3754-d5d515d0a2d5
# ╠═77906a42-354c-11eb-397e-9dca847ab964
# ╠═93696cf4-354d-11eb-0af3-9b27a8186fbd
# ╠═a91c8b80-354d-11eb-03ec-8f4386a69708
# ╠═86f0f21c-354d-11eb-2a32-299183f9118b
# ╠═666a8878-354d-11eb-22dc-e94a72ede22d
# ╟─cfc0e592-354d-11eb-2b06-458d198c38e5
# ╠═25789968-3552-11eb-1e25-35c0d1b1c461
# ╠═56daee6e-3554-11eb-37cd-49ad4ac56805
# ╠═56bb225a-3554-11eb-0062-6d03b8ea772f
# ╟─569a20f0-3554-11eb-0479-f1e583620682
# ╠═56790ed8-3554-11eb-3c51-b9b31a4c340e
