### A Pluto.jl notebook ###
# v0.12.15

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

# ╔═╡ 271a0ba8-3557-11eb-339e-d754bca4b807
using LinearAlgebra, SparseArrays, Plots, DifferentialEquations

# ╔═╡ a25cf69c-3555-11eb-0cc8-a5a605bd4581
md"## Heat equation with Neumann boundary conditions

The domain is given by the one-dimensional rod $\Omega = (0,L)$ and the time by $t \in (0,T).$ The rod has the physical properties
- thermal conductivity $\lambda$,
- specific heat capacity $c$,
- mass density $\rho$.

The diffusivity constant is given by $\alpha = \frac{\lambda}{c \rho}$. The heat equation is noted as

$\dot{\vartheta}(t,x) = \alpha \frac{\partial^2}{\partial x^2} \vartheta(t,x)$

for $(t,x) \in (0,T) \times \Omega$ with initial condition

$\vartheta(0,x) = \vartheta_{0}(x)$ 

for $x \in \overline{\Omega}$ and boundary conditions  

$\lambda ~ \left.\frac{\partial}{\partial x} \vartheta(\cdot, x) \right\rvert_{x = 0} \cdot \vec{n} ~=~ \lambda ~  \left. \frac{\partial}{\partial x} \vartheta(\cdot, x) \right\rvert_{x = L} \cdot \vec{n} ~=~ 0$

with outer [normal vector](https://en.wikipedia.org/wiki/Normal_(geometry)) $\vec{n}$ on the left or right boundary. Here, the normal vector $\vec{n} = -1$ on the left boundary at $x=0$, and  $\vec{n} = 1$ on the right boundary at $x=L$."

# ╔═╡ c1c415c2-3555-11eb-1957-4f6b5ce551e0
md"#### Code: implementing physical constants"

# ╔═╡ 909d656c-3556-11eb-0507-317403a6d74a
# Length of rod
L = 0.5

# ╔═╡ 125d530a-3557-11eb-15fe-51c2df3e812d
# Thermal conductivity
λ = 45.0 

# ╔═╡ 1581566c-3557-11eb-21a7-b936097d94f7
# Specific heat capacity
c = 480.0

# ╔═╡ 197db882-3557-11eb-183d-1d8ae847f443
# Mass density
ρ = 7800.0

# ╔═╡ 1e4c2b50-3557-11eb-1e2a-7f5d8dc28eb2
# Thermal diffusivity
α = λ/(c*ρ)

# ╔═╡ 24899264-3557-11eb-2c72-07ab0d214f4a
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

# ╔═╡ 2700c51c-3557-11eb-28dd-dd5e5c6073c4
# Number of grid elements
N = 101

# ╔═╡ 26e4c8f8-3557-11eb-0964-236aac988cfe
# Finite discretization
Δx = L/(N-1)   

# ╔═╡ 26cd7ffe-3557-11eb-261b-87ba7df942c5
# 1-dimensional grid
xspan = 0 : Δx : L

# ╔═╡ 26b11922-3557-11eb-1340-694c36380b5e
# Diffusion matrix
M = spdiagm(-1 => ones(N-1), 0 => -2*ones(N), 1 => ones(N-1));

# ╔═╡ 269cf7d2-3557-11eb-310b-2591cd6a22c9
Matrix(M)[1:5,1:5]

# ╔═╡ 2681a336-3557-11eb-1aaf-87de2d8f4711
md"### Boundary conditions

The spatial approximation of the heat equation 

$\frac{\partial^2 \vartheta(t,x)}{\partial x^2} \approx \frac{1}{\Delta x^2} \left( \vartheta(\cdot, x^{n-1}) - 2 ~ \vartheta(\cdot, x^{n}) + \vartheta(\cdot, x^{n+1}) \right)$

can not be evaluated directly at the grid points $x^{0}$ and $x^{N-1}$ because they depend on values from the not-existing grid points $x^{-1}$ and $x^{N}$. Therefore, the Neumann boundary condition is used to find $\vartheta(\cdot, x^{-1})$ and $\vartheta(\cdot, x^{N})$.

The Neumann boundary condition states that the heat flux at both sides is fixed as 

$\lambda \left. \frac{\partial \vartheta(\cdot,x)}{\partial x}  \right\rvert_{x = 0} \cdot \vec{n} = -1 \cdot \lambda \frac{\vartheta(\cdot, x^{1}) - \vartheta(\cdot, x^{-1})}{2 \Delta x}   = 0$

on the left boundary and 

$\lambda \left. \frac{\partial \vartheta(\cdot,x)}{\partial x}  \right\rvert_{x = L} \cdot \vec{n} = \lambda \frac{\vartheta(\cdot, x^{N}) - \vartheta(\cdot, x^{N-2})}{2 \Delta x}  = 0$

on the right boundary. Reformulating both equations lead to 

$\vartheta(\cdot, x^{-1}) = \vartheta(\cdot, x^{1}) \quad \text{and} \quad \vartheta(\cdot, x^{N}) = \vartheta(\cdot, x^{N-2}).$

On the left boundary one yields 

$\left. \frac{\partial^2 \vartheta(t,x)}{\partial x^2} \right\rvert_{x = 0} \approx \frac{1}{\Delta x^2} \left( - 2 ~ \vartheta(\cdot, x^{0}) + 2 ~ \vartheta(\cdot, x^{1}) \right)$

and analog on the right boundary

$\left. \frac{\partial^2 \vartheta(t,x)}{\partial x^2} \right\rvert_{x = 0} \approx \frac{1}{\Delta x^2} \left(2 ~ \vartheta(\cdot, x^{N-2}) - 2 ~ \vartheta(\cdot, x^{N-1}) \right).$

This means, there is no (induced/emitted) heat flux at both boundaries and thus, the diffusion matrix has the form 

$M = 
\begin{pmatrix}
-2 & 2 & 0 & \cdots &  & 0 \\
1 & -2 & 1 & 0  & \cdots & \vdots \\
0 & 1 & -2 & 1 & \ddots &  \\
 & & \ddots & \ddots & \ddots &  & \\
& & & 1 & -2 & 1 \\
& & &  & 2 & -2 \\
\end{pmatrix}.$"

# ╔═╡ 26672cfe-3557-11eb-0c59-a380e7746021
# First row
M[1,2] = 2;

# ╔═╡ 264b0bdc-3557-11eb-208a-4fd7b5aebc78
# Last row
M[end,end-1] = 2;

# ╔═╡ 262d874a-3557-11eb-0b14-7bc4207b4702
md"### Heat equation as ODE

After the spatial approximation the heat equation has a form of 

$\dot{\theta}(t) = \frac{\alpha}{\Delta x^2} ~ M ~ \theta(t)$

with M as diffusion matrix and $\quad \theta(t) = \left( \vartheta(t, x^{0}), \cdots, \vartheta(t, x^{N-1}) \right)^{\top}$. This is an ordinary differential equation (ODE) that can be solved with common solvers like forward Euler method or Runge-Kutta scheme."

# ╔═╡ 260b8f48-3557-11eb-10f9-1b3ce6fb5b3d
# Heat Equation as ODE
function heat_eq(dθ, θ, p, t)
    return dθ .= (1/Δx^2) * α * M * θ
end

# ╔═╡ 25e6bd08-3557-11eb-1e5a-bfc0dafb8967
md"### Initial Conditions

The initial data of the original heat equation is assumed as $\quad \vartheta_{0}(x) := m \left[ L~x - x^2 \right]~$ with  $m > 0$. This function is approximated to gain the initial conditions of the ODE 

$\theta(0) = \left( \vartheta_{0}(x^{0}), \cdots , \vartheta_{0}(x^{N-1}) \right)^{\top}.$
"

# ╔═╡ 25bd5f6c-3557-11eb-343e-035104baca78
#Initial conditions: θ(0) = m * (L * x - x^2)
function initial_data(x, params)
    # α = params[1] # unused here 
    L = params[2]
    m = params[3]

    return m*( L * x - x^2)
end

# ╔═╡ 259e1b5c-3557-11eb-1acb-490bb548bd7e
# Initial heat distribution
θ₀ = zeros(N);

# ╔═╡ 183acba8-3567-11eb-0c79-cb464e53ca96
# Amplification of initial data
@bind amp html"<input type='range' min='0' max='100' step='1'>"

# ╔═╡ 181893c6-3567-11eb-04c4-976b37dae40c
m = amp

# ╔═╡ 17fda6ce-3567-11eb-13dc-95a7282429d3
# Parameters
param = [α, L, amp]

# ╔═╡ 17e196b4-3567-11eb-36ee-6567a654e031
let

	for i = 1 : length(xspan)
		θ₀[i] = initial_data(xspan[i], param)
	end
	
	plot(xspan, θ₀, xlabel="Position x", ylabel="Temperature", legend=false)
end

# ╔═╡ 17c86858-3567-11eb-0b76-53c802d1168e
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

# ╔═╡ 17aacdaa-3567-11eb-3d11-f3f1c6d90751
# Upper limit of sampling time
ul = 0.5 * Δx^2 / α

# ╔═╡ 178d5ffe-3567-11eb-0190-f3c2331e3e86
# Sampling time
Δt = 0.8

# ╔═╡ 0058fe28-3568-11eb-05b0-bbbb624657a8
# Final time
@bind Tf html"<input type='range' min='0' max='1000' step='1'>"

# ╔═╡ 1747cdfe-3567-11eb-319d-171b7c1b69b7
# time discretization
tspan = (0.0, Tf)

# ╔═╡ 172485c4-3567-11eb-2b49-a9e0a4731acc
# Solving the ODE
let
	prob = ODEProblem( heat_eq, θ₀, tspan ) # ODE Problem
	sol = solve(prob,Euler(),dt=Δt,progress=true, save_everystep=false, save_start=true) # Solving the ODE
	plot(xspan, sol.u[2], xlabel = "Position x", ylabel="Temperature", legend=false)
end

# ╔═╡ Cell order:
# ╠═a25cf69c-3555-11eb-0cc8-a5a605bd4581
# ╠═c1c415c2-3555-11eb-1957-4f6b5ce551e0
# ╠═909d656c-3556-11eb-0507-317403a6d74a
# ╠═125d530a-3557-11eb-15fe-51c2df3e812d
# ╠═1581566c-3557-11eb-21a7-b936097d94f7
# ╠═197db882-3557-11eb-183d-1d8ae847f443
# ╠═1e4c2b50-3557-11eb-1e2a-7f5d8dc28eb2
# ╠═24899264-3557-11eb-2c72-07ab0d214f4a
# ╠═271a0ba8-3557-11eb-339e-d754bca4b807
# ╠═2700c51c-3557-11eb-28dd-dd5e5c6073c4
# ╠═26e4c8f8-3557-11eb-0964-236aac988cfe
# ╠═26cd7ffe-3557-11eb-261b-87ba7df942c5
# ╠═26b11922-3557-11eb-1340-694c36380b5e
# ╠═269cf7d2-3557-11eb-310b-2591cd6a22c9
# ╠═2681a336-3557-11eb-1aaf-87de2d8f4711
# ╠═26672cfe-3557-11eb-0c59-a380e7746021
# ╠═264b0bdc-3557-11eb-208a-4fd7b5aebc78
# ╠═262d874a-3557-11eb-0b14-7bc4207b4702
# ╠═260b8f48-3557-11eb-10f9-1b3ce6fb5b3d
# ╠═25e6bd08-3557-11eb-1e5a-bfc0dafb8967
# ╠═25bd5f6c-3557-11eb-343e-035104baca78
# ╠═259e1b5c-3557-11eb-1acb-490bb548bd7e
# ╠═183acba8-3567-11eb-0c79-cb464e53ca96
# ╠═181893c6-3567-11eb-04c4-976b37dae40c
# ╠═17fda6ce-3567-11eb-13dc-95a7282429d3
# ╠═17e196b4-3567-11eb-36ee-6567a654e031
# ╠═17c86858-3567-11eb-0b76-53c802d1168e
# ╠═17aacdaa-3567-11eb-3d11-f3f1c6d90751
# ╠═178d5ffe-3567-11eb-0190-f3c2331e3e86
# ╠═0058fe28-3568-11eb-05b0-bbbb624657a8
# ╠═1747cdfe-3567-11eb-319d-171b7c1b69b7
# ╠═172485c4-3567-11eb-2b49-a9e0a4731acc
