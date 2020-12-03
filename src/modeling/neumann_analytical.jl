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

# ╔═╡ d5dc9294-34a7-11eb-09b9-e36684398036
using Plots

# ╔═╡ f75649f2-33e9-11eb-0a65-0f48aa7209a6
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

$\left.\frac{\partial}{\partial x} \vartheta(\cdot, x) \right\rvert_{x = 0} \cdot \vec{n} ~=~ \left. \frac{\partial}{\partial x} \vartheta(\cdot, x) \right\rvert_{x = L} \cdot \vec{n} ~=~ 0$

with outer [normal vector](https://en.wikipedia.org/wiki/Normal_(geometry)) $\vec{n}$ on the left or right boundary. Here, the normal vector $\vec{n} = -1$ on the left boundary at $x=0$, and  $\vec{n} = 1$ on the right boundary at $x=L$."

# ╔═╡ 038013f8-34a6-11eb-37b1-69150fede3db
md"#### Code: implementing physical constants"

# ╔═╡ 0366d38c-34a6-11eb-2cfe-05e6f8edf17b
# Length of rod
L = 0.5

# ╔═╡ 03494950-34a6-11eb-32f5-7b485d899c42
# Final time
Tf = 10000.0  

# ╔═╡ 032aaeb6-34a6-11eb-2742-b12ab7682f17
# Thermal conductivity
λ = 45.0 

# ╔═╡ 03002010-34a6-11eb-2b2f-9dd02826187b
# Specific heat capacity
c = 480.0

# ╔═╡ 02e3e238-34a6-11eb-2d89-310870da991b
# Mass density
ρ = 7800.0

# ╔═╡ 44a6fa0c-34a6-11eb-2917-779a9aa42543
# Thermal diffusivity
α = λ/(c*ρ)

# ╔═╡ b1a58804-34a5-11eb-3d49-25b877f17b75
md"### Separation of variables

It is assumed that $\quad \vartheta(t,x) = f(t) ~ g(x) \quad$ holds. Applying [separation of variables](https://en.wikipedia.org/wiki/Separation_of_variables) approach into the heat equation, one yields

$\dot{\vartheta}(t,x) = \dot{f}(t) ~ g(x) =  \alpha ~ f(t) ~ \frac{\partial^2}{\partial x^2} g(x).$

Separating the functions results in

$\frac{\dot{f}(t)}{f(t)} = \alpha ~ \frac{g''(x)}{g(x)} = p.$


The ordinary differential equation $\quad \dot{f}(t) = p ~ f(t) \quad$ is solved by $\quad f(t) = C_{f} ~ \exp\left(p~t\right)$. And the second order ODE $\quad g''(x) = q ~ g(x) \quad$ with $q = \frac{p}{\alpha}$ is solved by 

$g(x) ~=~ A_{g} ~ \sin(\sqrt{-q} x) + B_{g} \cos(\sqrt{-q} x).$

### Eigenvalues and Eigenvectors

Applying the initial condition at the boundaries leads to 

$\left.\frac{\partial}{\partial x} \vartheta(0, x) \right\rvert_{x = 0} ~=~ f(0) ~ g'(0) = C_{f} ~ \sqrt{-q} ~ A_{g} ~=~ 0$

implying $A_{g}=0$, and 

$\left. \frac{\partial}{\partial x} \vartheta(\cdot, x) \right\rvert_{x = L} ~=~ f(0) g'(L) = C_{f} ~ \sqrt{-q} ~ (-B_{g}) ~  \sin(\sqrt{-q} x) ~=~ 0.$

Thus, the Eigenvalues are calculated by $\quad \sqrt{-q}_{i}~\sin\left( \sqrt{-q}_{i} L \right) = 0 \quad$ or equivalent $\quad q_{i} = 0$ for all $i \in \mathbb{N}$ or  $\quad q_{i} =  \left(\frac{i ~ \pi}{L}\right)^{2}$.

Variable $B_{i} = B_{g,i}$ as coefficient of spatial solution

$g_{i}(x) = B_{g,i} \cos\left( \sqrt{-q}_{i} x \right)$ 

is introduced and used to find the Eigenvectors or Eigenfunctions

$\Phi_{i}(x) ~=~  B_{i} ~ \cos\left( \sqrt{-q_{i}} x \right) ~=~  B_{i}  \cos\left( i ~ \pi ~ \frac{x}{L} \right)$

The Eigenvectors form an [orthonormal basis](https://en.wikipedia.org/wiki/Orthonormal_basis) and thus 

$\langle \Phi_{i}, \Phi_{j} \rangle ~=~ \int\limits_{0}^{L} B_{i}  \cos\left( i ~ \pi ~ \frac{x}{L} \right) ~ B_{j}  \cos\left( j ~ \pi ~ \frac{x}{L} \right) \operatorname{d}x ~=~ \delta_{i,j}$

with $\delta_{i,j} ~=~ 
	\begin{cases}
	1 \quad \text{for} \quad i = j,\\
	0 \quad \text{otherwise}.
	\end{cases}$

In the case $i = j$ one finds with 

$\langle \Phi_{i}, \Phi_{i} \rangle = B_{i}^2 \int\limits_{0}^{L}  \cos\left( i ~ \pi ~ \frac{x}{L} \right)^2 \operatorname{d}x = 1$ 

the coefficients $B_{i} ~=~ \sqrt{\frac{2}{L}}$.

In the case $i \neq j$ the equality 

$\langle \Phi_{i}, \Phi_{j} \rangle ~=~ B_{i} B_{j}  \int\limits_{0}^{L} \cos\left( i ~ \pi ~ \frac{x}{L} \right) ~  \cos\left( j ~ \pi ~ \frac{x}{L} \right) \operatorname{d}x = 0$ 

holds because $\int\limits_{0}^{L} \cos\left( i ~ \pi ~ \frac{x}{L} \right) ~  \cos\left( j ~ \pi ~ \frac{x}{L} \right) \operatorname{d}x = 0.$ 

Consequently, the Eigenvectors (or Eigenfunctions) are given by

$\Phi_{i}(x) = \sqrt{\frac{2}{L}} ~ \cos\left( i \pi \frac{x}{L}\right).$"

# ╔═╡ 137ff0a4-34a7-11eb-14cf-d5f7c5ad29f7
md"#### Code: Eigenvalues and Eigenvectors"

# ╔═╡ a1008ff4-34a9-11eb-2b19-270ae9e857ef
# Maximum index i
max_order = 10

# ╔═╡ 1365b150-34a7-11eb-3c6e-39d36506722f
# Eigenvalues
let
	i = 0 : 1 : max_order
	q = -1 * (i*pi/L).^2
	scatter(q, label="q", xaxis="Index i", yaxis="Eigenvalue")
end

# ╔═╡ 13426eaa-34a7-11eb-0ee6-f3282b7dcd20
@bind istep html"<input type='range' min='1' max='10' step='1'>"

# ╔═╡ 14490a74-34a8-11eb-0129-ed4370a3dc51
# Number of spatial elements
Nx = 2000

# ╔═╡ 142cef62-34a8-11eb-2df7-898cfde39332
# Finite discretization
Δx = L/(Nx-1)

# ╔═╡ 13ed2e20-34a8-11eb-233e-fd36ac4ceeca
# Discretized x positions
xspan = 0 : Δx : L

# ╔═╡ b1528564-34a5-11eb-20ec-3d3c2933cb36
let 
	qi = - (istep*pi/L)^2
	Ai = sqrt(2/L)
	Φi(x) = Ai * cos( sqrt(-qi)*x)
	plot(xspan, Φi.(xspan), xlabel="Position x", ylabel="Eigenvector Phi", legend=false)
end

# ╔═╡ 13cfca10-34a8-11eb-140e-47b3a1b55a56


# ╔═╡ 0412e44a-34a8-11eb-3c66-e18f3e9f8789


# ╔═╡ 8d3b4652-34a5-11eb-1f85-cb954b678ea0
md"### Initial Conditions

The preliminary solution of the heat equation is noted by

$\vartheta(t,x) = \eta_{0} + \sum\limits_{i = 1}^{\infty} T_{i}(t) ~ \Phi_{i}(x) ~=~ \eta_{0} + \sum\limits_{i = 1}^{\infty} C_{f,i} ~ \exp\left( \alpha ~ q_{i} ~ t \right) \sqrt{\frac{2}{L}} \cos\left( i \pi \frac{x}{L} \right)$

**Remark:** the addional coefficient $\eta_{0}$ describes the behaviour of $T_{i}(t) ~ \Phi_{i}(x)$ at index $i=0$. 

The coefficients $C_{f,i}$ are found with the initial conditions

$\vartheta(0,x) = \eta_{0} + \sum\limits_{i = 1}^{\infty} C_{i} ~ \Phi_{i}(x) ~=~ \vartheta_{0}(x)$

Both sides of the last equation are multiplied by the Eigenvector $\Phi_{i}$ and the inner product is computed elementwise to yield

$C_{i} \langle \Phi_{i}, \Phi_{i} \rangle ~=~ \langle \Phi_{i}, \vartheta_{0} -  \eta_{0}\rangle$

or equivalent 

$C_{i} = \frac{\langle \Phi_{i}, \vartheta_{0}  -  \eta_{0}\rangle}{\langle \Phi_{i}, \Phi_{i} \rangle} = \langle \Phi_{i}, \vartheta_{0}  -  \eta_{0} \rangle$ 

because of the orthonormal property $\langle \Phi_{i}, \Phi_{i} \rangle = 1$.

Now, an example initial data is assumed with $\quad \vartheta_{0}(x) := m \left[ L~x - x^2 \right] \quad$ with  $m > 0$.

This iniial data leads to the coefficients

$C_{i} ~=~ \sqrt{\frac{2}{L}}  ~ (-1) ~ m ~ L ~ \left(\frac{L}{i \pi}\right)^2 \left[\left(-1\right)^{i} + 1 \right] ~=~
\begin{cases}
-2 ~ m ~ L ~ \sqrt{\frac{2}{L}} ~ \left(\frac{L}{i \pi}\right)^2 \quad &\text{if} \quad i \text{ is even} \\
 0 \quad &\text{if} \quad i \text{ is odd} 
\end{cases}$"

# ╔═╡ 3ba830e4-34a9-11eb-002b-4952d88fbaad
md"#### Code: Coefficients"

# ╔═╡ 482641dc-34a9-11eb-1794-792fd102aeab
# Amplification (initial data)
@bind amp html"<input type='range' min='0' max='100' step='0.1'>"

# ╔═╡ 48426022-34a9-11eb-2baa-8dfb566d9cce
let
	max_order = 10
	i = 1 : 1 : max_order
	Cf(i) = sqrt(2/L)*amp*(L/(i*pi))^3 * (1 - (-1)^i)
	scatter(i, Cf.(i), label="Cf", xaxis="i")
end

# ╔═╡ 3b902486-34a9-11eb-1529-45521466a6c8


# ╔═╡ 8d573ab0-34a5-11eb-1f0f-59bd9cb87092
md"###  Final solution

Finally, the solution of the heat equation with Dirichlet boundary conditions is given by

$\vartheta(t,x) = \eta_{0} - 2 m ~ \frac{L^2}{\pi^2} \sum\limits_{i = 1}^{\infty} \frac{1}{i^2} \left[\left(-1\right)^{i} + 1 \right] ~ \exp\left(- \alpha \left(\frac{i \pi }{L}\right)^2 ~ t \right) ~  \cos\left( i \pi \frac{x}{L}\right)$ 

and with odd index $i = 2 k$ one yields 

$\vartheta(t,x) = \eta_{0} - m ~ \frac{L^2}{\pi^2}  ~ \sum\limits_{k = 1}^{\infty} \frac{1}{k^2} ~ \exp\left(- \alpha ~ 4 ~ \left(\frac{k \pi }{L}\right)^2 ~ t \right) ~  \cos\left( 2 ~k ~ \pi \frac{x}{L}\right)$"


# ╔═╡ 8ce0b9a8-34a5-11eb-16e2-b1858f9d89c4
md"#### Code: spatial and temporal approximation"

# ╔═╡ 8cc10270-34a5-11eb-0958-b566b8f488f9
function temporal_approx(elem, t, params)
    k = elem
    α = params[1]
    L = params[2]
    return exp( -4*α *(k*π/L)^2 *t)
end

# ╔═╡ 8ca3d774-34a5-11eb-1915-b31f19aad0d5
function spatial_approx(elem, x, params)
    k = elem
    L = params[2]
    return cos( 2*k* π*(x/L))
end

# ╔═╡ 8d1cb34a-34a5-11eb-0a91-cf59727f69b7
md"#### Coefficient $\eta_{0}$

The value of offset $\eta_{0}$ is found by computing the temperature at $(t,x)=(0,0)$ with

$\vartheta(t,x) = \eta_{0} - m ~ \frac{L^2}{\pi^2}  ~ \sum\limits_{k = 1}^{\infty} \frac{1}{k^2} = \vartheta_{0}(0).$

The sum of $\frac{1}{k^2}$ corresponds with the [Riemannian Zeta function](https://en.wikipedia.org/wiki/Riemann_zeta_function) at position 2  and can be used in

$\eta_{0} = \vartheta_{0}(0) + m ~ \frac{L^2}{\pi^2}  ~ \sum\limits_{k = 1}^{\infty} \frac{1}{k^2} = \vartheta_{0}(0) + m ~ \frac{L^2}{\pi^2}  ~ \zeta(2).$

According to the [Basel problem](https://en.wikipedia.org/wiki/Basel_problem), the Zeta function $\zeta(2) = \frac{\pi^2}{6}$ and consequently the coeffient is found with

$\eta_{0} \approx  \vartheta_{0}(0) + m ~ \frac{L^2}{6}.$"

# ╔═╡ 8d00d274-34a5-11eb-23e4-9b9a6999abe9
# Cofficient η₀
η₀ = amp * L^2/6

# ╔═╡ 8c81adb4-34a5-11eb-0221-2f9093992726
md"#### Code: Heat equation and solution"

# ╔═╡ 8c62f89c-34a5-11eb-3136-35c0b89ca4ad
function heat_eq(t,x, order, params)
    α = params[1]
    L = params[2]
    m = params[3]

    sum = 0;
    for i = 1 : order
        sum = sum + (1/i^2) * temporal_approx(i, t, params) * spatial_approx(i, x, params)
    end

    η = m * L^2/6;

    return η  -  m*L^2/(π^2) * sum
end

# ╔═╡ 8c46f208-34a5-11eb-3c16-37149eb6a2c7
# Parameter
params = [α, L, amp]

# ╔═╡ 8beb67c8-34a5-11eb-2d1c-9ba0b4dc71d9
# Number of time steps
Nt = 1001

# ╔═╡ 8bc809ac-34a5-11eb-3fe5-2756444c0e52
# Temporal discretization
Δt = Tf / (Nt - 1)

# ╔═╡ 8bad1f72-34a5-11eb-11b9-fddc7e984050
tspan = 0 : Δt : Tf

# ╔═╡ 8c2ae3c8-34a5-11eb-0a68-950ae952a525
# Solution array
θsol = zeros(length(xspan), length(tspan));

# ╔═╡ 8c0c6054-34a5-11eb-06e1-8335ef9544fb
for i = 1 : length(tspan)
    for k = 1 : length(xspan)
        θsol[k,i] = heat_eq(tspan[i], xspan[k], max_order, params)
    end
end

# ╔═╡ 8b92737c-34a5-11eb-28b1-3dc2e0aa2a64
md"#### Visualisation: spatial distribution"

# ╔═╡ 8b73e0b8-34a5-11eb-048c-892c233c0b05
@bind tstep html"<input type='range' min='1' max='1001' step='1'>"

# ╔═╡ 8b58e428-34a5-11eb-14a8-ed421bee29ba
Time = tspan[tstep]

# ╔═╡ 8b3e3880-34a5-11eb-3dbe-09e2003a6803
plot(xspan, θsol[:,tstep], xlabel="Position x", ylabel="Temperature", legend=false)

# ╔═╡ 5b4d530a-34ab-11eb-3bfb-7985e1617652
md"#### Visualisation: evolution in time"

# ╔═╡ 5ae93d48-34ab-11eb-0444-fd8ebb0108a3
@bind xstep html"<input type='range' min='1' max='2000' step='1'>"

# ╔═╡ ea54dbdc-33fb-11eb-3955-eb0e29d4e3a2
Position = xspan[xstep]

# ╔═╡ 69b29af4-34ab-11eb-26d3-a58d6f015fa1
plot(tspan, θsol[xstep,:], xlabel="Time t", ylabel="Temperature", legend=false)

# ╔═╡ Cell order:
# ╠═f75649f2-33e9-11eb-0a65-0f48aa7209a6
# ╠═038013f8-34a6-11eb-37b1-69150fede3db
# ╠═0366d38c-34a6-11eb-2cfe-05e6f8edf17b
# ╠═03494950-34a6-11eb-32f5-7b485d899c42
# ╠═032aaeb6-34a6-11eb-2742-b12ab7682f17
# ╠═03002010-34a6-11eb-2b2f-9dd02826187b
# ╠═02e3e238-34a6-11eb-2d89-310870da991b
# ╠═44a6fa0c-34a6-11eb-2917-779a9aa42543
# ╟─b1a58804-34a5-11eb-3d49-25b877f17b75
# ╟─137ff0a4-34a7-11eb-14cf-d5f7c5ad29f7
# ╠═d5dc9294-34a7-11eb-09b9-e36684398036
# ╠═a1008ff4-34a9-11eb-2b19-270ae9e857ef
# ╠═1365b150-34a7-11eb-3c6e-39d36506722f
# ╠═13426eaa-34a7-11eb-0ee6-f3282b7dcd20
# ╠═14490a74-34a8-11eb-0129-ed4370a3dc51
# ╠═142cef62-34a8-11eb-2df7-898cfde39332
# ╠═13ed2e20-34a8-11eb-233e-fd36ac4ceeca
# ╠═b1528564-34a5-11eb-20ec-3d3c2933cb36
# ╠═13cfca10-34a8-11eb-140e-47b3a1b55a56
# ╠═0412e44a-34a8-11eb-3c66-e18f3e9f8789
# ╟─8d3b4652-34a5-11eb-1f85-cb954b678ea0
# ╠═3ba830e4-34a9-11eb-002b-4952d88fbaad
# ╠═482641dc-34a9-11eb-1794-792fd102aeab
# ╠═48426022-34a9-11eb-2baa-8dfb566d9cce
# ╠═3b902486-34a9-11eb-1529-45521466a6c8
# ╠═8d573ab0-34a5-11eb-1f0f-59bd9cb87092
# ╠═8ce0b9a8-34a5-11eb-16e2-b1858f9d89c4
# ╠═8cc10270-34a5-11eb-0958-b566b8f488f9
# ╠═8ca3d774-34a5-11eb-1915-b31f19aad0d5
# ╠═8d1cb34a-34a5-11eb-0a91-cf59727f69b7
# ╠═8d00d274-34a5-11eb-23e4-9b9a6999abe9
# ╠═8c81adb4-34a5-11eb-0221-2f9093992726
# ╠═8c62f89c-34a5-11eb-3136-35c0b89ca4ad
# ╠═8c46f208-34a5-11eb-3c16-37149eb6a2c7
# ╠═8beb67c8-34a5-11eb-2d1c-9ba0b4dc71d9
# ╠═8bc809ac-34a5-11eb-3fe5-2756444c0e52
# ╠═8bad1f72-34a5-11eb-11b9-fddc7e984050
# ╠═8c2ae3c8-34a5-11eb-0a68-950ae952a525
# ╠═8c0c6054-34a5-11eb-06e1-8335ef9544fb
# ╠═8b92737c-34a5-11eb-28b1-3dc2e0aa2a64
# ╠═8b73e0b8-34a5-11eb-048c-892c233c0b05
# ╠═8b58e428-34a5-11eb-14a8-ed421bee29ba
# ╠═8b3e3880-34a5-11eb-3dbe-09e2003a6803
# ╠═5b4d530a-34ab-11eb-3bfb-7985e1617652
# ╠═5ae93d48-34ab-11eb-0444-fd8ebb0108a3
# ╠═ea54dbdc-33fb-11eb-3955-eb0e29d4e3a2
# ╠═69b29af4-34ab-11eb-26d3-a58d6f015fa1
