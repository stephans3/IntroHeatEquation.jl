### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 7c62a2c4-9159-11eb-1b54-452017ea90d9
using LinearAlgebra, SparseArrays

# ╔═╡ 7ad16586-9162-11eb-1d94-5188eb8d4c5d
using DifferentialEquations, Plots

# ╔═╡ 44df6004-9158-11eb-1a63-651d7dca957b
md"# Heat equation with Robin boundary conditions 

Consider a one-dimensional rod $\Omega = (0,L)$ and time $t \in (0,T).$ 
"

# ╔═╡ aff71d6e-9158-11eb-0559-cb7334f4a624
L = 0.1 # Length of rod

# ╔═╡ 3771dab8-9159-11eb-3a69-cd9b94a10399
Tf = 60.0 # Final simulation time

# ╔═╡ 44dbf81c-9159-11eb-2920-27ee60a7db9a
md"The rod has the physical properties
- thermal conductivity $\lambda$,
- specific heat capacity $c$,
- mass density $\rho$.

The diffusivity constant is given by $\alpha = \frac{\lambda}{c \rho}$. The heat equation is noted as

$\dot{\vartheta}(t,x) = \alpha \frac{\partial^2}{\partial x^2} \vartheta(t,x)$

for $(t,x) \in (0,T) \times \Omega$ with initial condition

$\vartheta(0,x) = \vartheta_{0}(x)$ 

for $x \in \overline{\Omega}$." 

# ╔═╡ 86a40dc8-9158-11eb-16a3-1f776d6bde4e
begin
	λ = 45.0;    # Thermal conductivity
	ρ = 7800.0;  # Mass density
	cap = 480.0; # Specific heat capacitivity
end

# ╔═╡ 9e8ee0a2-9158-11eb-0210-4513221b7823
α = λ/(cap * ρ)   # Diffusivity

# ╔═╡ a4f1e1e2-9158-11eb-0754-a70399a341ba
md"## Boundary conditions  

On both sides of the rod, heat transfer

$\phi_{t}(t,x) = - h ~ (\vartheta(t,x) - \vartheta_{amb})$

and heat radiation

$\phi_{r}(t,x) = - \epsilon \varrho ~ (\vartheta(t,x)^4 - \vartheta_{amb}^4)$

is assumed to operate. Parameter $\vartheta_{amb}$ describes the ambient temperature, $h>0$ is called the heat transfer coefficient, $\epsilon \in (0,1)$ is called emissivity and $\varrho\approx 5.67 \cdot 10^{-8}$ is known as Stefan-Boltzmann constant. For simplicity, $k=\epsilon \cdot \varrho$. 

The sum of $\phi_{t}$ and $\phi_{r}$ 

$\phi_{out}(t,x) = \phi_{t}(t,x) + \phi_{r}(t,x) = - h ~ (\vartheta(t,x) - \vartheta_{amb}) - k ~ (\vartheta(t,x)^4 - \vartheta_{amb}^4)$

is the flux from the rod to the environment. Thus, the natural (Robin) boundary conditions are noted as 

$\lambda ~ \left.\frac{\partial}{\partial x} \vartheta(\cdot, x) \right\rvert_{x = 0} \cdot \vec{n} ~=~ \lambda ~  \left. \frac{\partial}{\partial x} \vartheta(\cdot, x) \right\rvert_{x = L} \cdot \vec{n} ~=~ \phi_{out}(t,x)$

with outer [normal vector](https://en.wikipedia.org/wiki/Normal_(geometry)) $\vec{n}$ on the left or right boundary. Here, the normal vector $\vec{n} = -1$ on the left boundary at $x=0$, and  $\vec{n} = 1$ on the right boundary at $x=L$."

# ╔═╡ 7c9e9a4a-9159-11eb-3dda-59d2aa15f32a
begin
	h = 10.0; # Heat transfer coefficient
	ϵ = 0.6;  # Emissivity
	sb = 5.67*10^(-8) # Stefan-Boltzmann constant
	k = ϵ * sb; 	  # Radiation coefficient
end

# ╔═╡ c5da798e-9160-11eb-1d16-75ce792ef2fe
θamb = 298.0 # Ambient temperature in Kelvin

# ╔═╡ 7c80218c-9159-11eb-301a-95a428fb76d9
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

# ╔═╡ 3487893a-915c-11eb-3a04-1fdd31009413
N = 101 # Number of grid elements

# ╔═╡ 39663578-915c-11eb-108c-2718b706729e
Δx = L/(N-1) # Finite discretization   

# ╔═╡ 3468092a-915c-11eb-1719-d1b9b04b68e0
# Diffusion matrix
M = spdiagm(-1 => ones(N-1), 0 => -2*ones(N), 1 => ones(N-1));

# ╔═╡ 344c1bb6-915c-11eb-2be7-35f340f4caa8
Matrix(M)[1:5,1:5]

# ╔═╡ 342fe676-915c-11eb-3f22-1d20d5b502ed
md"### Approximated boundary conditions

The spatial approximation of the heat equation 

$\frac{\partial^2 \vartheta(t,x)}{\partial x^2} \approx \frac{1}{\Delta x^2} \left( \vartheta(\cdot, x^{n-1}) - 2 ~ \vartheta(\cdot, x^{n}) + \vartheta(\cdot, x^{n+1}) \right)$

can not be evaluated directly at the grid points $x^{0}$ and $x^{N-1}$ because they depend on values from the not-existing grid points $x^{-1}$ and $x^{N}$. Therefore, the Neumann boundary condition is used to find $\vartheta(\cdot, x^{-1})$ and $\vartheta(\cdot, x^{N})$.

The Robin boundary condition states that the heat flux at both sides is described by

$\lambda \left. \frac{\partial \vartheta(\cdot,x)}{\partial x}  \right\rvert_{x = 0} \cdot \vec{n} = -1 \cdot \lambda \frac{\vartheta(\cdot, x^{1}) - \vartheta(\cdot, x^{-1})}{2 \Delta x}   = \phi_{out}(t,x^{0})$

on the left boundary and 

$\lambda \left. \frac{\partial \vartheta(\cdot,x)}{\partial x}  \right\rvert_{x = L} \cdot \vec{n} = \lambda \frac{\vartheta(\cdot, x^{N}) - \vartheta(\cdot, x^{N-2})}{2 \Delta x}  = \phi_{out}(t,x^{N-1})$

on the right boundary. Reformulating both equations lead to 

$\vartheta(\cdot, x^{-1}) = \vartheta(\cdot, x^{1}) + 2 \Delta x ~ \phi_{out}(t,x^{0})$

and

$\vartheta(\cdot, x^{N}) = \vartheta(\cdot, x^{N-2})  + 2 \Delta x ~ \phi_{out}(t,x^{N-1}).$

On the left boundary one yields 

$\left. \frac{\partial^2 \vartheta(t,x)}{\partial x^2} \right\rvert_{x = 0} \approx \frac{1}{\Delta x^2} \left( - 2 ~ \vartheta(\cdot, x^{0}) + 2 ~ \vartheta(\cdot, x^{1}) \right) + \frac{2}{\Delta x} \phi_{out}(t,x^{0})$

and analog on the right boundary

$\left. \frac{\partial^2 \vartheta(t,x)}{\partial x^2} \right\rvert_{x = 0} \approx \frac{1}{\Delta x^2} \left(2 ~ \vartheta(\cdot, x^{N-2}) - 2 ~ \vartheta(\cdot, x^{N-1}) \right) + \frac{2}{\Delta x} \phi_{out}(t,x^{N-1}).$

The diffusion matrix and the approximated flux are noted as

$M = 
\begin{pmatrix}
-2 & 2 & 0 & \cdots &  & 0 \\
1 & -2 & 1 & 0  & \cdots & \vdots \\
0 & 1 & -2 & 1 & \ddots &  \\
 & & \ddots & \ddots & \ddots &  & \\
& & & 1 & -2 & 1 \\
& & &  & 2 & -2 \\
\end{pmatrix}, \quad

\Phi_{out}(t) = 
\begin{pmatrix}
\phi_{out}(t, x^{0}) \\
0 \\
\vdots \\
\\
0 \\
\phi_{out}(t, x^{N-1}) \\
\end{pmatrix}.$
"


# ╔═╡ 3414b70c-915c-11eb-2606-4d1ad425d3ab
# First row
M[1,2] = 2;

# ╔═╡ 33cd51a0-915c-11eb-3c62-efc367e08126
# Last row
M[end,end-1] = 2;

# ╔═╡ 33af0164-915c-11eb-30f0-313d926c2dae
md"### Heat equation as ODE

After the spatial approximation the heat equation has a form of 

$\dot{\theta}(t) =  \frac{\alpha}{\Delta x^2} ~ M ~ \theta(t) + \frac{2 \alpha}{\Delta x} ~ \Phi_{out}(t)$

with M as diffusion matrix, temperature vector $\quad \theta(t) = \left( \vartheta(t, x^{0}), \cdots, \vartheta(t, x^{N-1}) \right)^{\top}$. This is an ordinary differential equation (ODE) that can be solved with common solvers like forward Euler method or Runge-Kutta scheme."

# ╔═╡ 328a1ed0-9160-11eb-0892-7d0a9d05a54d
# Heat Equation as ODE
function heat_eq!(dθ, θ, p, t)
	
	N = size(θ)
	Φout = zeros(N);
	
	Φout[1] = -h * (θ[1] - θamb) - k*(θ[1]^4 - θamb^4)

	Φout[end] = -h * (θ[end] - θamb) - k*(θ[end]^4 - θamb^4)

    dθ .= (1/Δx^2) * α * M * θ + (2/Δx) * α * Φout
end

# ╔═╡ 66d5aee4-9161-11eb-1ef1-2193d9c2f6a2
md"### Initial Conditions

For simplicity, the initial data of the heat equation is assumed as $\vartheta_{0}(x)= 1000$ Kelvin for all $x \in \left[0,L\right]$. This function is approximated to gain the initial conditions of the ODE 

$\theta(0) = \left( 10^3, \cdots , 10^3 \right)^{\top}.$
"

# ╔═╡ 91a272be-9162-11eb-1046-8dd9e33653e5
θ₀ = 10^3 * ones(N)

# ╔═╡ 1571f3f2-9162-11eb-0b63-ff4e25a07572
md"The upper limit of the sampling time has to 

$\Delta t < \frac{1}{2} \frac{\Delta x^2}{\alpha}.$

to guarantee numerical stability."

# ╔═╡ f85a2d4a-9161-11eb-0e0f-851c0750832a
Δt_ul = (0.5*Δx^2)/α # upper limit of sampling time

# ╔═╡ 5aa78f74-9162-11eb-1bd4-c96f2b90e586
Δt = 10^(-2) # Sampling time

# ╔═╡ 66ba8b22-9162-11eb-13b0-b14dfcead641
md"## Simulation"

# ╔═╡ 861018de-9162-11eb-2a01-2354acd868e2
tspan = (0.0, Tf)

# ╔═╡ e0e033e8-9162-11eb-3fe6-05565991ce77
# 1-dimensional grid
xspan = 0 : Δx : L

# ╔═╡ 01b373f0-9163-11eb-3f35-3f645825720f
prob = ODEProblem( heat_eq!, θ₀, tspan ) # ODE Problem

# ╔═╡ 0103888c-9163-11eb-16ae-371cd6d04c45
sol = solve(prob,Euler(),dt=Δt,progress=true, save_everystep=false, save_start=true) # Solving the ODE

# ╔═╡ 08dcb90e-9163-11eb-2ff0-7129add62934
plot(xspan, sol.u[2], xlabel = "Position x", ylabel="Temperature", legend=false)

# ╔═╡ Cell order:
# ╟─44df6004-9158-11eb-1a63-651d7dca957b
# ╠═aff71d6e-9158-11eb-0559-cb7334f4a624
# ╠═3771dab8-9159-11eb-3a69-cd9b94a10399
# ╟─44dbf81c-9159-11eb-2920-27ee60a7db9a
# ╠═86a40dc8-9158-11eb-16a3-1f776d6bde4e
# ╠═9e8ee0a2-9158-11eb-0210-4513221b7823
# ╟─a4f1e1e2-9158-11eb-0754-a70399a341ba
# ╠═7c9e9a4a-9159-11eb-3dda-59d2aa15f32a
# ╠═c5da798e-9160-11eb-1d16-75ce792ef2fe
# ╟─7c80218c-9159-11eb-301a-95a428fb76d9
# ╠═7c62a2c4-9159-11eb-1b54-452017ea90d9
# ╠═3487893a-915c-11eb-3a04-1fdd31009413
# ╠═39663578-915c-11eb-108c-2718b706729e
# ╠═3468092a-915c-11eb-1719-d1b9b04b68e0
# ╠═344c1bb6-915c-11eb-2be7-35f340f4caa8
# ╟─342fe676-915c-11eb-3f22-1d20d5b502ed
# ╠═3414b70c-915c-11eb-2606-4d1ad425d3ab
# ╠═33cd51a0-915c-11eb-3c62-efc367e08126
# ╟─33af0164-915c-11eb-30f0-313d926c2dae
# ╠═328a1ed0-9160-11eb-0892-7d0a9d05a54d
# ╟─66d5aee4-9161-11eb-1ef1-2193d9c2f6a2
# ╠═91a272be-9162-11eb-1046-8dd9e33653e5
# ╟─1571f3f2-9162-11eb-0b63-ff4e25a07572
# ╠═f85a2d4a-9161-11eb-0e0f-851c0750832a
# ╠═5aa78f74-9162-11eb-1bd4-c96f2b90e586
# ╠═66ba8b22-9162-11eb-13b0-b14dfcead641
# ╠═7ad16586-9162-11eb-1d94-5188eb8d4c5d
# ╠═861018de-9162-11eb-2a01-2354acd868e2
# ╠═e0e033e8-9162-11eb-3fe6-05565991ce77
# ╠═01b373f0-9163-11eb-3f35-3f645825720f
# ╠═0103888c-9163-11eb-16ae-371cd6d04c45
# ╠═08dcb90e-9163-11eb-2ff0-7129add62934
