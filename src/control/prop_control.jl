### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 892fb737-6fc0-4438-b586-039a440b01f9
using LinearAlgebra, SparseArrays

# ╔═╡ 069cf190-77cc-4757-9db3-1945ddd07438
using DifferentialEquations, Plots

# ╔═╡ 7bc3e874-b7be-11eb-0792-cf4e21e8c8b2
md"# Heat equation with Proportional control

Consider a one-dimensional rod $\Omega = (0,L)$ and time $t \in (0,T).$ 
"

# ╔═╡ 1158de2e-b7a6-43f3-a4e6-2ef9598a4ca5
L = 0.1 # Length of rod

# ╔═╡ 30cebec3-9c70-4c6a-94a8-fdcbe6cda783
Tf = 2000.0 # Final simulation time

# ╔═╡ 93fb0f90-c454-410f-aa5e-cd8e6c055dc0
md"The rod has the physical properties
- thermal conductivity $\lambda$,
- specific heat capacity $c$,
- mass density $\rho$.

The diffusivity constant is given by $\alpha = \frac{\lambda}{c \rho}$. The heat equation is noted as

$\dot{\vartheta}(t,x) = \alpha \frac{\partial^2}{\partial x^2} \vartheta(t,x)$

for $(t,x) \in (0,T) \times \Omega$ with initial condition

$\vartheta(0,x) = \vartheta_{0}(x)$ 

for $x \in \overline{\Omega}$." 

# ╔═╡ b294bbfc-45e7-4a15-a839-5ceb56746be7
begin
	λ = 45.0;    # Thermal conductivity
	ρ = 7800.0;  # Mass density
	cap = 480.0; # Specific heat capacitivity
end

# ╔═╡ 543b3f13-4cef-46f1-9141-172b81e7d48b
α = λ/(cap * ρ)   # Diffusivity

# ╔═╡ 1deb117d-21f7-4329-a7db-2be26e43129a
md"## Boundary conditions  

Only on the right side of the rod, heat transfer

$\phi_{t}(t,x=L) = - h (\vartheta(t,x=L) - \vartheta_{amb})$

and heat radiation

$\phi_{r}(t,x=L) = - \epsilon \varrho (\vartheta(t,x=L)^4 - \vartheta_{amb}^4)$

is assumed to operate. Parameter $\vartheta_{amb}$ describes the ambient temperature, $h>0$ is called the heat transfer coefficient, $\epsilon \in (0,1)$ is called emissivity and $\varrho\approx 5.67 \cdot 10^{-8}$ is known as Stefan-Boltzmann constant. For simplicity, $k=\epsilon \cdot \varrho$. 

The sum of $\phi_{t}$ and $\phi_{r}$ 

$\phi_{out}(t,L) = \phi_{t}(t,L) + \phi_{r}(t,L) = - h (\vartheta(t,L) - \vartheta_{amb}) - k (\vartheta(t,L)^4 - \vartheta_{amb}^4)$

is the flux from the rod to the environment."

# ╔═╡ 78292b42-a203-496c-a143-def3ff128f69
begin
	h = 10.0; # Heat transfer coefficient
	ϵ = 0.6;  # Emissivity
	sb = 5.67*10^(-8) # Stefan-Boltzmann constant
	k = ϵ * sb; 	  # Radiation coefficient
end

# ╔═╡ 0ecea33b-bf24-48eb-b4b4-bbf288e62e6b
θamb = 298.0 # Ambient temperature in Kelvin

# ╔═╡ 3a059f15-d850-468a-bb1d-f42ef5b24609
ϕout(θ) = -h * (θ - θamb)  - k*(θ^4 - θamb^4) 

# ╔═╡ 7202aa44-b393-4e0b-9c34-2a9b9376194c
md"### Induced heat

On the left side of the rod, a heat source is assumed as 

$\phi_{in}(t,x=0) = b ~ u(t)$

with constant $b>0$, input signal 

$u(t) = K_{p} ~ e(t)$

and error

$e(t) = y_{ref}(t) - y(t).$

The output is considered as the temperature on the right side, e.g.

$y(t) = \vartheta(t, x=L)$

and the reference is as an arbitrary fixed value 

$y_{ref}(t) = y_{ref} > 0.$

Thus, the natural (Robin) boundary conditions are noted as 

$\lambda ~ \left.\frac{\partial}{\partial x} \vartheta(\cdot, x) \right\rvert_{x = 0} \cdot \vec{n} ~=~ \phi_{in}(t,0) = b ~ u(t)$

and

$\lambda ~  \left. \frac{\partial}{\partial x} \vartheta(\cdot, x) \right\rvert_{x = L} \cdot \vec{n} ~=~ \phi_{out}(t,x) = - h (\vartheta(t,L) - \vartheta_{amb}) - k (\vartheta(t,L)^4 - \vartheta_{amb}^4)$

with outer [normal vector](https://en.wikipedia.org/wiki/Normal_(geometry)) $\vec{n}$ on the left or right boundary. Here, the normal vector $\vec{n} = -1$ on the left boundary at $x=0$, and  $\vec{n} = 1$ on the right boundary at $x=L$."

# ╔═╡ 2775433b-d8b1-48be-9b60-08641db15803
b = 1;

# ╔═╡ 4c59909b-ffb6-4dfc-94ab-142d41deeeab
Kp = 10^3  # Proportional gain

# ╔═╡ 4f0c3671-d3bf-4e36-a1f3-8e8d0ac0e454
yref = 400.0 # Reference temperature

# ╔═╡ 23343fa0-6033-498d-b089-d43e4382e332
u_in(err) = Kp * err # input signals

# ╔═╡ daa17037-fd28-460c-97f5-826742b9a419
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

# ╔═╡ 0309bf5e-4033-4832-86ed-146aea0863c0
N = 101 # Number of grid elements

# ╔═╡ 5c9ec071-6140-4af6-a418-df79c10daf5e
Δx = L/(N-1) # Finite discretization   

# ╔═╡ 905337c5-9eee-4b67-b215-bfd7ca905e4c
# Diffusion matrix
M = spdiagm(-1 => ones(N-1), 0 => -2*ones(N), 1 => ones(N-1));

# ╔═╡ f0e9f197-303b-48a6-a540-9b11513050d4
Matrix(M)[1:5,1:5]

# ╔═╡ 0e155080-2cfe-45c0-8166-ba91f95d6f8b
md"### Approximated boundary conditions

The spatial approximation of the heat equation 

$\frac{\partial^2 \vartheta(t,x)}{\partial x^2} \approx \frac{1}{\Delta x^2} \left( \vartheta(\cdot, x^{n-1}) - 2 ~ \vartheta(\cdot, x^{n}) + \vartheta(\cdot, x^{n+1}) \right)$

can not be evaluated directly at the grid points $x^{0}$ and $x^{N-1}$ because they depend on values from the not-existing grid points $x^{-1}$ and $x^{N}$. Therefore, the Neumann boundary condition is used to find $\vartheta(\cdot, x^{-1})$ and $\vartheta(\cdot, x^{N})$.

The Robin boundary condition states that the heat flux at both sides is described by

$\lambda \left. \frac{\partial \vartheta(\cdot,x)}{\partial x}  \right\rvert_{x = 0} \cdot \vec{n} = -1 \cdot \lambda \frac{\vartheta(\cdot, x^{1}) - \vartheta(\cdot, x^{-1})}{2 \Delta x}   = \phi_{in}(t,x^{0})$

on the left boundary and 

$\lambda \left. \frac{\partial \vartheta(\cdot,x)}{\partial x}  \right\rvert_{x = L} \cdot \vec{n} = \lambda \frac{\vartheta(\cdot, x^{N}) - \vartheta(\cdot, x^{N-2})}{2 \Delta x}  = \phi_{out}(t,x^{N-1})$

on the right boundary. Reformulating both equations lead to 

$\vartheta(\cdot, x^{-1}) = \vartheta(\cdot, x^{1}) + \frac{2 \Delta x}{\lambda} ~ \phi_{in}(t,x^{0})$

and

$\vartheta(\cdot, x^{N}) = \vartheta(\cdot, x^{N-2})  + \frac{2 \Delta x}{\lambda} ~ \phi_{out}(t,x^{N-1}).$

On the left boundary one yields 

$\left. \frac{\partial^2 \vartheta(t,x)}{\partial x^2} \right\rvert_{x = 0} \approx \frac{1}{\Delta x^2} \left( - 2 ~ \vartheta(\cdot, x^{0}) + 2 ~ \vartheta(\cdot, x^{1}) \right) + \frac{2}{\lambda ~ \Delta x} \phi_{in}(t,x^{0})$

and analog on the right boundary

$\left. \frac{\partial^2 \vartheta(t,x)}{\partial x^2} \right\rvert_{x = 0} \approx \frac{1}{\Delta x^2} \left(2 ~ \vartheta(\cdot, x^{N-2}) - 2 ~ \vartheta(\cdot, x^{N-1}) \right) + \frac{2}{\lambda ~ \Delta x} \phi_{out}(t,x^{N-1}).$

The induced and emitted flux are noted as 

$\Phi_{in}(t) = 
\begin{pmatrix}
b  \\
0 \\
\vdots \\
0 \\
\end{pmatrix}
~ u(t)

\quad \text{and} \quad

\Phi_{out}(t) = 
\begin{pmatrix}
0 \\
\vdots \\
0 \\
1 \\
\end{pmatrix}
~\phi_{out}(t, x^{N-1})$

and the diffusion matrix is described by

$M = 
\begin{pmatrix}
-2 & 2 & 0 & \cdots &  & 0 \\
1 & -2 & 1 & 0  & \cdots & \vdots \\
0 & 1 & -2 & 1 & \ddots &  \\
 & & \ddots & \ddots & \ddots &  & \\
& & & 1 & -2 & 1 \\
& & &  & 2 & -2 \\
\end{pmatrix}.$
"


# ╔═╡ c61962e8-86bd-45ed-a62d-ffc812d8e35e
# First row
M[1,2] = 2;

# ╔═╡ 33baab03-fbc1-4944-bc16-e4022b051c69
# Last row
M[end,end-1] = 2;

# ╔═╡ facdbabf-6e2b-418a-a07e-36188244e585
md"### Heat equation in State-space representation

After the spatial approximation the heat equation has a form of 

$\dot{\theta}(t) =  \frac{\alpha}{\Delta x^2} ~ M ~ \theta(t) + \frac{2 \alpha}{\Delta x} ~ \Phi_{in}(t) + \frac{2 \alpha}{\Delta x} ~ \Phi_{out}(t) \tag{1}$

with M as diffusion matrix, temperature vector $\quad \theta(t) = \left( \vartheta(t, x^{0}), \cdots, \vartheta(t, x^{N-1}) \right)^{\top}$. 

The original ODE $(1)$ is recasted in [state-space representation](https://en.wikipedia.org/wiki/State-space_representation) as

$\dot{\theta}(t) = A ~ \theta(t) + B ~ u(t) + E ~ w(t)$
$y(t) = C \theta(t)$

with 
- system matrix $A =  \frac{\alpha}{\Delta x^2} ~ M$, 
- input vector 
$B = \frac{2 ~ \alpha}{\lambda ~ \Delta x} 
\begin{pmatrix}
b \\
0 \\
\vdots \\
0
\end{pmatrix},$

- disturbance vector 
$E = \frac{2 ~ \alpha}{\lambda ~ \Delta x} 
\begin{pmatrix}
0 \\
\vdots \\
0 \\
1
\end{pmatrix},$

- emitted flux $w(t) =  \phi_{out}(t,x^{N-1})$ and
- output vector $C = (0, \cdots, 0, 1)$

to describe the input-state-output behaviour including the integral controller. The emitted heat is considered as a disturbance. The PI controller is introduced as 

$u(t) = K_{p} ~ e(t)$ 

with error 
$e(t) = y_{ref}(t) - y(t) = y_{ref}(t) - C ~ \theta(t)$


The thermal dynamics with PI control is described by

$\dot{\theta}(t) = A ~ \theta(t) + B ~ K_{p} ~ e(t) + E ~ w(t).$

This is an ordinary differential equation (ODE) that can be solved with common solvers like forward Euler method or Runge-Kutta scheme.

"

# ╔═╡ 5ef85c8f-cb69-42df-a098-ab98903a5f32
A = α/(Δx^2) * M

# ╔═╡ 25bb45ab-55af-4b8b-a11b-05be68631b78
begin
	B = spzeros(N); # Input vector
	B[1] = 2 * (α/(λ*Δx)) * b
	E = spzeros(N); # Disturbance vector
	E[end] = 2 * (α/(λ*Δx)) * b
	C = spzeros(1,N); # Output vector
	C[end] = 1;
end

# ╔═╡ 81a2d204-7386-4ed7-af89-7c076940817b
w(θ) = ϕout(θ)

# ╔═╡ 113d8c87-f4e6-4f1b-b188-879e7ad1e146
function heat_eq!(dx,x,p,t)
    global u_hist
  
    err = yref - x[end]
    u = u_in(err)
   

    dx .= A * x + B * u + E * ϕout( x[end] )    # Integration of temperatures
   
  end

# ╔═╡ d9138271-79fc-4b7c-8501-e8d8f3304cd5
md"### Initial Conditions

For simplicity, the initial data of the heat equation is assumed as $\vartheta_{0}(x)= 1000$ Kelvin for all $x \in \left[0,L\right]$. This function is approximated to gain the initial conditions of the ODE 

$\theta(0) = \left( 300, \cdots , 300 \right)^{\top}.$
"

# ╔═╡ db6248c7-24b8-4f44-82d3-4892203297b5
θ₀ = 300.0 * ones(N)

# ╔═╡ 837632eb-b4db-47ec-8e4a-32a0fe81060b
md"## Simulation"

# ╔═╡ 3f382fe5-2aa6-4d82-b939-977985ae8bbb
Δt = 10^(-2) # Sampling time

# ╔═╡ 02ba374a-15b4-4481-baf6-c1a3ce45ed5e
tspan = (0.0, Tf)

# ╔═╡ 0db82493-4dda-4992-aa5b-63734b556a90
prob = ODEProblem( heat_eq!, θ₀, tspan ) # ODE Problem

# ╔═╡ f7dcb8a6-0eb4-4762-942e-53362d4fe1ba
sol = solve(prob,Euler(),dt=Δt,progress=true, saveat=1.0, save_start=true); # Solving the ODE

# ╔═╡ 9c38b2b3-1622-41f3-a2c2-ff79b552e454
# 1-dimensional grid
xspan = 0 : Δx : L

# ╔═╡ ba66c778-6c1c-4e32-934d-e8aca0eb1c66
heatmap(sol.t, xspan, sol[1:end,1:end], xaxis="Time [s]", yaxis="Position [m]", title="Evolution of temperature")

# ╔═╡ 3dcfed19-8235-4309-b1c7-8d4e30ba6701
plot(sol.t,[sol[1,:], sol[end-1,:]], label=["Left" "Right"], title="Temperature at the left/right end", legend=:bottomright)

# ╔═╡ Cell order:
# ╟─7bc3e874-b7be-11eb-0792-cf4e21e8c8b2
# ╠═1158de2e-b7a6-43f3-a4e6-2ef9598a4ca5
# ╠═30cebec3-9c70-4c6a-94a8-fdcbe6cda783
# ╟─93fb0f90-c454-410f-aa5e-cd8e6c055dc0
# ╠═b294bbfc-45e7-4a15-a839-5ceb56746be7
# ╠═543b3f13-4cef-46f1-9141-172b81e7d48b
# ╟─1deb117d-21f7-4329-a7db-2be26e43129a
# ╠═78292b42-a203-496c-a143-def3ff128f69
# ╠═0ecea33b-bf24-48eb-b4b4-bbf288e62e6b
# ╠═3a059f15-d850-468a-bb1d-f42ef5b24609
# ╟─7202aa44-b393-4e0b-9c34-2a9b9376194c
# ╠═2775433b-d8b1-48be-9b60-08641db15803
# ╠═4c59909b-ffb6-4dfc-94ab-142d41deeeab
# ╠═4f0c3671-d3bf-4e36-a1f3-8e8d0ac0e454
# ╠═23343fa0-6033-498d-b089-d43e4382e332
# ╟─daa17037-fd28-460c-97f5-826742b9a419
# ╠═892fb737-6fc0-4438-b586-039a440b01f9
# ╠═0309bf5e-4033-4832-86ed-146aea0863c0
# ╠═5c9ec071-6140-4af6-a418-df79c10daf5e
# ╠═905337c5-9eee-4b67-b215-bfd7ca905e4c
# ╠═f0e9f197-303b-48a6-a540-9b11513050d4
# ╟─0e155080-2cfe-45c0-8166-ba91f95d6f8b
# ╠═c61962e8-86bd-45ed-a62d-ffc812d8e35e
# ╠═33baab03-fbc1-4944-bc16-e4022b051c69
# ╟─facdbabf-6e2b-418a-a07e-36188244e585
# ╠═5ef85c8f-cb69-42df-a098-ab98903a5f32
# ╠═25bb45ab-55af-4b8b-a11b-05be68631b78
# ╠═81a2d204-7386-4ed7-af89-7c076940817b
# ╠═113d8c87-f4e6-4f1b-b188-879e7ad1e146
# ╟─d9138271-79fc-4b7c-8501-e8d8f3304cd5
# ╠═db6248c7-24b8-4f44-82d3-4892203297b5
# ╠═837632eb-b4db-47ec-8e4a-32a0fe81060b
# ╠═069cf190-77cc-4757-9db3-1945ddd07438
# ╠═3f382fe5-2aa6-4d82-b939-977985ae8bbb
# ╠═02ba374a-15b4-4481-baf6-c1a3ce45ed5e
# ╠═0db82493-4dda-4992-aa5b-63734b556a90
# ╠═f7dcb8a6-0eb4-4762-942e-53362d4fe1ba
# ╠═9c38b2b3-1622-41f3-a2c2-ff79b552e454
# ╠═ba66c778-6c1c-4e32-934d-e8aca0eb1c66
# ╠═3dcfed19-8235-4309-b1c7-8d4e30ba6701
