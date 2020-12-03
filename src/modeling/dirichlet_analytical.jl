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

# ╔═╡ e982a7dc-3486-11eb-3d00-c581b6204ba4
using Plots

# ╔═╡ 3a6b8abe-33b6-11eb-1b17-e5499677cbbb
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

# ╔═╡ 067b4ae4-34a6-11eb-0ab0-8baaaafebe91
md"#### Code: implementing physical constants"

# ╔═╡ de9f801a-33bc-11eb-1f48-759320188ecb
# Length of rod
L = 0.5

# ╔═╡ da31f9f4-33bc-11eb-3deb-c752080f0364
# Final time
Tf = 10000.0  

# ╔═╡ 31b447dc-3484-11eb-09a4-fb4fdfae50ea
# Thermal conductivity
λ = 45.0 

# ╔═╡ 2c493a82-3484-11eb-38a0-3f9186ef4651
# Specific heat capacity
c = 480.0

# ╔═╡ fc7a2b6e-33b5-11eb-0394-2f6610423fc4
# Mass density
ρ = 7800.0

# ╔═╡ d8ddb41c-3484-11eb-1dbd-5923e303ce44
# Thermal diffusivity
α = λ/(c*ρ)

# ╔═╡ 1702ea48-33b6-11eb-225f-690266bff18b
md"### Separation of Variables

It is assumed that $\quad \vartheta(t,x) = f(t) ~ g(x) \quad$ holds. Applying [separation of variables](https://en.wikipedia.org/wiki/Separation_of_variables) approach into the heat equation, one yields

$\dot{\vartheta}(t,x) = \dot{f}(t) ~ g(x) =  \alpha ~ f(t) ~ \frac{\partial^2}{\partial x^2} g(x).$

Separating the functions results in

$\frac{\dot{f}(t)}{f(t)} = \alpha ~ \frac{g''(x)}{g(x)} = p.$

The ordinary differential equation $\quad \dot{f}(t) = p ~ f(t) \quad$ is solved by $\quad f(t) = C_{f} ~ \exp\left(p~t\right)$. And the second order ODE $\quad g''(x) = q ~ g(x) \quad$ with $q = \frac{p}{\alpha}$ is solved by 

$g(x) ~=~ A_{g} ~ \sin(\sqrt{-q} x) + B_{g} \cos(\sqrt{-q} x).$"

# ╔═╡ 44439e16-3489-11eb-044c-8fd1eb5d9d96
md"#### Code: Temporal Approximation"

# ╔═╡ cd1c81cc-34a1-11eb-236e-034aa56d0630
# Number of time steps
Nt = 1001

# ╔═╡ 71b83c18-3488-11eb-345f-a7f32779e4be
# Temporal discretization
Δt = Tf / (Nt - 1)

# ╔═╡ 077e8d24-3489-11eb-163d-191a00eead4c
tspan = 0 : Δt : Tf

# ╔═╡ 71358fac-3488-11eb-0534-a7d68c35547e
# Assumed constant
p = -0.001

# ╔═╡ 70d2ee30-3488-11eb-1ace-8b8a2d141df2
# Assumed gain
Cf = 1.0

# ╔═╡ d117be54-3488-11eb-1b04-fd594e88d629
let 
	f(t) = Cf * exp(p*t)
	plot( tspan, f.(tspan), xlabel="Time t", ylabel="Temporal approx f", legend=false) 
end

# ╔═╡ 625eb216-3489-11eb-3d04-6ff714e33ad0
md"#### Code: Spatial Approximation"

# ╔═╡ ce7b84a8-3485-11eb-3d54-7d8382f09525
# Number of spatial elements
Nx = 2000

# ╔═╡ ce57afa4-3485-11eb-01c2-af85042c6bb8
# Finite discretization
Δx = L/(Nx-1)

# ╔═╡ a44877ea-3485-11eb-3bf9-2fb6f8e904e0
# Discretized x positions
xspan = 0 : Δx : L

# ╔═╡ a428eb3c-3485-11eb-2123-af0ace14b3d2
# Coefficients
Ag = 1.0

# ╔═╡ 8530d152-3485-11eb-148d-ab2016be12de
Bg = 1.0

# ╔═╡ 43f049c6-3486-11eb-12cd-6f9a5c9851f0
# Assumed Eigenvalue
q = p/α

# ╔═╡ 5deef368-3486-11eb-3913-f39efa98c3ea
# Spatial approximation
let 
	g(x) = Ag * sin(sqrt(-q)*x) + Bg *cos(sqrt(-q)*x)
	plot( xspan, g.(xspan),  xlabel="Position x", ylabel="Spatial approx g", legend=false)
end

# ╔═╡ e91452e4-3486-11eb-0331-e33181d11d27


# ╔═╡ e8c3a556-3486-11eb-1f98-f5f885b83b9f


# ╔═╡ 84e6c5fa-3485-11eb-01d2-05b27470788d
md"### Eigenvalues and Eigenvectors

Applying the initial condition at the boundaries leads to 

$\vartheta(0,0) = f(0) ~ g(0) ~=~ C_{f} ~ B_{g} ~=~ 0$ 

implying $B_{g} = 0$, and

$\vartheta(0,L) = f(0) ~ g(L) ~=~ C_{f} ~ B_{g} \sin\left( \sqrt{-q} L \right) ~=~ 0.$ 

Thus, the Eigenvalues are calculated by $\quad \sin\left( \sqrt{-q}_{i} L \right) = 0 \quad$ or equivalent $\quad q_{i} =  -\left(\frac{i ~ \pi}{L}\right)^{2}$.

Variable $A_{i} = A_{g,i}$ as coefficient of spatial solution

$g_{i}(x) = A_{g,i} \sin\left( \sqrt{-q}_{i} x \right)$ 

is introduced and used to find the Eigenvectors or Eigenfunctions

$\Phi_{i}(x) ~=~ A_{i} ~ \sin\left( \sqrt{-q_{i}} x \right) ~=~ A_{i} ~ \sin\left( i ~ \pi ~ \frac{x}{L} \right)$

"

# ╔═╡ 04bb2616-34a1-11eb-33a8-bba7eac56bc9
# Maximum index i
max_order = 10

# ╔═╡ 3ea9ace2-348a-11eb-2eaa-31beeaec252a
# Eigenvalues
let
	i = 0 : 1 : max_order
	q = -1 * (i*pi/L).^2
	scatter(q, label="q", xaxis="Index i", yaxis="Eigenvalue")
end

# ╔═╡ f4035e12-348a-11eb-0d97-85087cc013da
@bind istep html"<input type='range' min='1' max='10' step='1'>"

# ╔═╡ 3e85b56c-348a-11eb-34e9-f1138c4cfc04
let 
	qi = - (istep*pi/L)^2
	Ai = 1.0
	Φi(x) = Ai * sin( sqrt(-qi)*x)
	plot(xspan, Φi.(xspan), xlabel="Position x", ylabel="Eigenvector Phi")
end

# ╔═╡ 3684e970-348a-11eb-0e5a-abf04a015492


# ╔═╡ 849684be-3485-11eb-16dd-d9ba541937de
md"#### Orthonormal basis
The Eigenvectors form an [orthonormal basis](https://en.wikipedia.org/wiki/Orthonormal_basis) and thus 

$\langle \Phi_{i}, \Phi_{j} \rangle ~=~ \int\limits_{0}^{L} A_{i}  \sin\left( i ~ \pi ~ \frac{x}{L} \right) ~ A_{j}  \sin\left( j ~ \pi ~ \frac{x}{L} \right) \operatorname{d}x ~=~ \delta_{i,j}$

with $\delta_{i,j} ~=~ 
	\begin{cases}
	1 \quad \text{for} \quad i = j,\\
	0 \quad \text{otherwise}.
	\end{cases}$

In the case $i = j$ one finds with 

$\langle \Phi_{i}, \Phi_{i} \rangle = A_{i}^2 \int\limits_{0}^{L}  \sin\left( i ~ \pi ~ \frac{x}{L} \right)^2 \operatorname{d}x = 1$ 

the coefficients $A_{i} ~=~ \sqrt{\frac{2}{L}}$.

In the case $i \neq j$ the equality 

$\langle \Phi_{i}, \Phi_{j} \rangle ~=~ A_{i} A_{j}  \int\limits_{0}^{L} \sin\left( i ~ \pi ~ \frac{x}{L} \right) ~  \sin\left( j ~ \pi ~ \frac{x}{L} \right) \operatorname{d}x = 0$ 

holds because $\int\limits_{0}^{L} \sin\left( i ~ \pi ~ \frac{x}{L} \right) ~  \sin\left( j ~ \pi ~ \frac{x}{L} \right) \operatorname{d}x = 0.$ 

Consequently, the Eigenvectors (or Eigenfunctions) are given by

$\Phi_{i}(x) = \sqrt{\frac{2}{L}} ~ \sin\left( i \pi \frac{x}{L}\right).$"

# ╔═╡ f65eac2e-348b-11eb-3672-59d4d7031578
# Found Eigenvectors
let 
	Ai = sqrt(2/L)
	Φi(x) = Ai * sin( istep*pi*x/L)
	plot(xspan, Φi.(xspan),  xlabel="Position x", ylabel="Eigenvector Phi", legend=false)
end

# ╔═╡ f6402fe2-348b-11eb-2340-7fa6145176b4


# ╔═╡ 1bbfd9c4-33b6-11eb-1f64-5ff08ab43cd4
md"### Initial Conditions

The preliminary solution of the heat equation is noted by

$\vartheta(t,x) = \sum\limits_{i = 1}^{\infty} T_{i}(t) ~ \Phi_{i}(x) ~=~ \sum\limits_{i = 1}^{\infty} C_{f,i} ~ \exp\left( \alpha ~ q_{i} ~ t \right) \sqrt{\frac{2}{L}} \sin\left( i \pi \frac{x}{L} \right)$

The coefficients $C_{f,i}$ are found with the initial conditions

$\vartheta(0,x) = \sum\limits_{i = 1}^{\infty} C_{i} ~ \Phi_{i}(x) ~=~ \vartheta_{0}(x).$

Both sides of the last equation are multiplied by the Eigenvector $\Phi_{i}$ and the inner product is computed elementwise to yield

$C_{i} \langle \Phi_{i}, \Phi_{i} \rangle ~=~ \langle \Phi_{i}, \vartheta_{0} \rangle$

or equivalent 

$C_{i} = \frac{\langle \Phi_{i}, \vartheta_{0} \rangle}{\langle \Phi_{i}, \Phi_{i} \rangle} = \langle \Phi_{i}, \vartheta_{0} \rangle$ 

because of the orthonormal property $\langle \Phi_{i}, \Phi_{i} \rangle = 1$.

Now, an example initial data is assumed with $\quad \vartheta_{0}(x) := m \left[ L~x - x^2 \right] \quad$ with  $m > 0$.

This initial data leads to the coefficients

$C_{i} = \sqrt{\frac{2}{L}}  ~ 2 ~ m ~ \left(\frac{L}{i \pi}\right)^3 \left[1 - \left(-1\right)^{i} \right] ~=~
	\begin{cases}
	0 \quad &\text{if} \quad i \text{ is even,} \\
	4 ~ m ~ \sqrt{\frac{2}{L}} ~ \left(\frac{L}{i \pi}\right)^3 \quad &\text{if} \quad i \text{ is odd.} 
\end{cases}$"

# ╔═╡ 2bf32e60-349a-11eb-3c46-2dd16252ad4a
md"#### Code: Initial data"

# ╔═╡ 2b13b0d2-349a-11eb-143a-d324c45eb6a4
@bind amp html"<input type='range' min='0' max='100' step='0.1'>"

# ╔═╡ 60e47128-349b-11eb-2c62-cd740543936b
amp # Amplification/gain m

# ╔═╡ 89501474-349a-11eb-1c05-8914b1459b1a
let 
	ϑ₀(x) = amp * (L*x - x^2)
	plot(xspan, ϑ₀.(xspan), xlabel="Position x", ylabel="Temperature", legend=false)
end

# ╔═╡ d1b9a26e-349c-11eb-1717-210e2fa10582
md"#### Code: Coefficients"

# ╔═╡ db7f3228-349c-11eb-2a1f-ff88bdf7766a
let
	i = 1 : 1 : max_order
	Cf(i) = sqrt(2/L)*amp*(L/(i*pi))^3 * (1 - (-1)^i)
	scatter(i, Cf.(i), label="Cf", xaxis="i")
end

# ╔═╡ 2af10a1e-349a-11eb-258f-db2f85abb895
md"###  Final solution

Finally, the solution of the heat equation with Dirichlet boundary conditions is given by

$\vartheta(t,x) = 4 m ~ \frac{L^2}{\pi^3} \sum\limits_{i = 1}^{\infty} \frac{1}{i^3} \left[1 - \left(-1\right)^{i} \right] \exp\left(- \alpha \left(\frac{i \pi }{L}\right)^2 ~ t \right) \sin\left( i \pi \frac{x}{L}\right)$ 

and with odd index $i = (2 k - 1)$ one yields 

$\vartheta(t,x) = 8 m ~ \frac{L^2}{\pi^3}  ~ \sum\limits_{k = 1}^{\infty} \frac{1}{(2~k - 1)^3} ~ \exp\left(- \alpha ~ \left(2 ~ k - 1\right)^2 ~ \left(\frac{\pi }{L}\right)^2 ~ t \right) ~  \sin\left( \left(2 ~ k - 1\right) ~ \pi \frac{x}{L}\right)$"

# ╔═╡ 00952038-34a0-11eb-3a2e-b3e796c8d076
md"### Simulation

#### Temporal approximation"

# ╔═╡ ffc0ab5a-349f-11eb-2fe9-c1bc1ccc3e90
function temporal_approx(elem, t, params)
    k = elem
    α = params[1]
    L = params[2]
    return exp( -α * (2*k - 1)^2*(π/L)^2 *t)
end

# ╔═╡ 298cf7fe-34a0-11eb-0df2-413b1defdb13
md"#### Spatial approximation"

# ╔═╡ 296eaf88-34a0-11eb-359e-29f76a278305
function spatial_approx(elem, x, params)
    k = elem
    L = params[2]
    return sin((2*k - 1) * π*(x/L))
end

# ╔═╡ 2954b574-34a0-11eb-1efb-ff60574966f1
md"#### Heat equation"

# ╔═╡ 293a1610-34a0-11eb-0cd3-e1fc499e3b00
function heat_eqn(t,x, order, params)
    α = params[1]
    L = params[2]
    m = params[3]

    sum = 0;
    for i = 1 : order
        sum = sum + (1/(2*i - 1)^3) * temporal_approx(i, t, params) * spatial_approx(i, x, params)
    end

    return 8*m*L^2/(π^3) * sum
end

# ╔═╡ 52883704-34a0-11eb-1dca-47a7e02d61cd
# Parameter
params = [α, L, amp]

# ╔═╡ 51c22e2e-34a0-11eb-06ad-f9a163aad9e2
# Solution array
θsol = zeros(length(xspan), length(tspan));

# ╔═╡ 69dd8602-34a0-11eb-3b27-c142552aa78b
md"Computation of solution"

# ╔═╡ 696f527a-34a0-11eb-180b-939677013514
for i = 1 : length(tspan)
    for k = 1 : length(xspan)
        θsol[k,i] = heat_eqn(tspan[i], xspan[k], max_order, params)
    end
end

# ╔═╡ 7918ea80-34a0-11eb-06f4-3bc07e22c81c
md"#### Visualisation: spatial distribution"

# ╔═╡ 7868e7fc-34a0-11eb-1506-47ab55e7cf6c
@bind tstep html"<input type='range' min='1' max='1001' step='1'>"

# ╔═╡ d1a2460c-34a2-11eb-24f1-694d916e6fbb
Time = tspan[tstep]

# ╔═╡ 784c423c-34a0-11eb-327e-e93f5492dd1f
plot(xspan, θsol[:,tstep], xlabel="Position x", ylabel="Temperature", legend=false)

# ╔═╡ 782d483c-34a0-11eb-39bd-c1c98e011761
md"#### Visualisation: evolution in time"

# ╔═╡ 1821a288-34a2-11eb-33a0-9bef5aa30455
@bind xstep html"<input type='range' min='1' max='2000' step='1'>"

# ╔═╡ 47a780fe-34a2-11eb-156d-ab4c511095f6
Position = xspan[xstep]

# ╔═╡ 180812aa-34a2-11eb-0efa-d1b0911b9722
plot(tspan, θsol[xstep,:], xlabel="Time t", ylabel="Temperature", legend=false)

# ╔═╡ Cell order:
# ╠═3a6b8abe-33b6-11eb-1b17-e5499677cbbb
# ╠═067b4ae4-34a6-11eb-0ab0-8baaaafebe91
# ╠═de9f801a-33bc-11eb-1f48-759320188ecb
# ╠═da31f9f4-33bc-11eb-3deb-c752080f0364
# ╠═31b447dc-3484-11eb-09a4-fb4fdfae50ea
# ╠═2c493a82-3484-11eb-38a0-3f9186ef4651
# ╠═fc7a2b6e-33b5-11eb-0394-2f6610423fc4
# ╠═d8ddb41c-3484-11eb-1dbd-5923e303ce44
# ╟─1702ea48-33b6-11eb-225f-690266bff18b
# ╟─44439e16-3489-11eb-044c-8fd1eb5d9d96
# ╠═cd1c81cc-34a1-11eb-236e-034aa56d0630
# ╠═71b83c18-3488-11eb-345f-a7f32779e4be
# ╠═077e8d24-3489-11eb-163d-191a00eead4c
# ╠═71358fac-3488-11eb-0534-a7d68c35547e
# ╟─70d2ee30-3488-11eb-1ace-8b8a2d141df2
# ╠═d117be54-3488-11eb-1b04-fd594e88d629
# ╠═625eb216-3489-11eb-3d04-6ff714e33ad0
# ╠═ce7b84a8-3485-11eb-3d54-7d8382f09525
# ╠═ce57afa4-3485-11eb-01c2-af85042c6bb8
# ╠═a44877ea-3485-11eb-3bf9-2fb6f8e904e0
# ╠═a428eb3c-3485-11eb-2123-af0ace14b3d2
# ╠═8530d152-3485-11eb-148d-ab2016be12de
# ╠═43f049c6-3486-11eb-12cd-6f9a5c9851f0
# ╠═e982a7dc-3486-11eb-3d00-c581b6204ba4
# ╠═5deef368-3486-11eb-3913-f39efa98c3ea
# ╠═e91452e4-3486-11eb-0331-e33181d11d27
# ╠═e8c3a556-3486-11eb-1f98-f5f885b83b9f
# ╠═84e6c5fa-3485-11eb-01d2-05b27470788d
# ╠═04bb2616-34a1-11eb-33a8-bba7eac56bc9
# ╠═3ea9ace2-348a-11eb-2eaa-31beeaec252a
# ╠═f4035e12-348a-11eb-0d97-85087cc013da
# ╠═3e85b56c-348a-11eb-34e9-f1138c4cfc04
# ╠═3684e970-348a-11eb-0e5a-abf04a015492
# ╠═849684be-3485-11eb-16dd-d9ba541937de
# ╠═f65eac2e-348b-11eb-3672-59d4d7031578
# ╠═f6402fe2-348b-11eb-2340-7fa6145176b4
# ╠═1bbfd9c4-33b6-11eb-1f64-5ff08ab43cd4
# ╠═2bf32e60-349a-11eb-3c46-2dd16252ad4a
# ╠═2b13b0d2-349a-11eb-143a-d324c45eb6a4
# ╠═60e47128-349b-11eb-2c62-cd740543936b
# ╠═89501474-349a-11eb-1c05-8914b1459b1a
# ╠═d1b9a26e-349c-11eb-1717-210e2fa10582
# ╠═db7f3228-349c-11eb-2a1f-ff88bdf7766a
# ╠═2af10a1e-349a-11eb-258f-db2f85abb895
# ╠═00952038-34a0-11eb-3a2e-b3e796c8d076
# ╠═ffc0ab5a-349f-11eb-2fe9-c1bc1ccc3e90
# ╠═298cf7fe-34a0-11eb-0df2-413b1defdb13
# ╠═296eaf88-34a0-11eb-359e-29f76a278305
# ╠═2954b574-34a0-11eb-1efb-ff60574966f1
# ╠═293a1610-34a0-11eb-0cd3-e1fc499e3b00
# ╠═52883704-34a0-11eb-1dca-47a7e02d61cd
# ╠═51c22e2e-34a0-11eb-06ad-f9a163aad9e2
# ╠═69dd8602-34a0-11eb-3b27-c142552aa78b
# ╠═696f527a-34a0-11eb-180b-939677013514
# ╠═7918ea80-34a0-11eb-06f4-3bc07e22c81c
# ╠═7868e7fc-34a0-11eb-1506-47ab55e7cf6c
# ╠═d1a2460c-34a2-11eb-24f1-694d916e6fbb
# ╠═784c423c-34a0-11eb-327e-e93f5492dd1f
# ╠═782d483c-34a0-11eb-39bd-c1c98e011761
# ╠═1821a288-34a2-11eb-33a0-9bef5aa30455
# ╠═47a780fe-34a2-11eb-156d-ab4c511095f6
# ╠═180812aa-34a2-11eb-0efa-d1b0911b9722
