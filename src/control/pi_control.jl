### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ a986c796-9166-11eb-2a58-cb9aeb65d440
using LinearAlgebra, SparseArrays

# ╔═╡ 501b67a8-921e-11eb-15fe-33cebed6a427
using DifferentialEquations, Plots

# ╔═╡ af2dd4bc-9163-11eb-1991-ef8b7897519a
md"# Heat equation with PI control

Consider a one-dimensional rod $\Omega = (0,L)$ and time $t \in (0,T).$ 
"

# ╔═╡ c0607ff8-9163-11eb-1610-25eab3f919c4
L = 0.1 # Length of rod

# ╔═╡ c040a860-9163-11eb-308a-332e7b19a186
Tf = 2000.0 # Final simulation time

# ╔═╡ c021f60e-9163-11eb-34e3-6d8faa141b4a
md"The rod has the physical properties
- thermal conductivity $\lambda$,
- specific heat capacity $c$,
- mass density $\rho$.

The diffusivity constant is given by $\alpha = \frac{\lambda}{c \rho}$. The heat equation is noted as

$\dot{\vartheta}(t,x) = \alpha \frac{\partial^2}{\partial x^2} \vartheta(t,x)$

for $(t,x) \in (0,T) \times \Omega$ with initial condition

$\vartheta(0,x) = \vartheta_{0}(x)$ 

for $x \in \overline{\Omega}$." 

# ╔═╡ c006067e-9163-11eb-141d-23e5206e7e41
begin
	λ = 45.0;    # Thermal conductivity
	ρ = 7800.0;  # Mass density
	cap = 480.0; # Specific heat capacitivity
end

# ╔═╡ bfec0b16-9163-11eb-0296-9db3fc9aaa79
α = λ/(cap * ρ)   # Diffusivity

# ╔═╡ bfd345ea-9163-11eb-25e4-a7ed1ac077aa
md"## Boundary conditions  

Only on the right side of the rod, heat transfer

$\phi_{t}(t,x=L) = - h (\vartheta(t,x=L) - \vartheta_{amb})$

and heat radiation

$\phi_{r}(t,x=L) = - \epsilon \varrho (\vartheta(t,x=L)^4 - \vartheta_{amb}^4)$

is assumed to operate. Parameter $\vartheta_{amb}$ describes the ambient temperature, $h>0$ is called the heat transfer coefficient, $\epsilon \in (0,1)$ is called emissivity and $\varrho\approx 5.67 \cdot 10^{-8}$ is known as Stefan-Boltzmann constant. For simplicity, $k=\epsilon \cdot \varrho$. 

The sum of $\phi_{t}$ and $\phi_{r}$ 

$\phi_{out}(t,L) = \phi_{t}(t,L) + \phi_{r}(t,L) = - h (\vartheta(t,L) - \vartheta_{amb}) - k (\vartheta(t,L)^4 - \vartheta_{amb}^4)$

is the flux from the rod to the environment."

# ╔═╡ 303f96d6-9166-11eb-3814-ef21a86a4972
begin
	h = 10.0; # Heat transfer coefficient
	ϵ = 0.6;  # Emissivity
	sb = 5.67*10^(-8) # Stefan-Boltzmann constant
	k = ϵ * sb; 	  # Radiation coefficient
end

# ╔═╡ 33f52982-9166-11eb-3de5-c1323f49271d
θamb = 298.0 # Ambient temperature in Kelvin

# ╔═╡ 6b4763f8-921c-11eb-1b6f-23505e52ed31
ϕout(θ) = -h * (θ - θamb)  - k*(θ^4 - θamb^4) 

# ╔═╡ bfb9cef8-9163-11eb-18d5-efa318c120d2
md"### Induced heat

On the left side of the rod, a heat source is assumed as 

$\phi_{in}(t,x=0) = b ~ u(t)$

with constant $b>0$, input signal 

$u(t) = K_{p} ~ e(t) + K_{i} \int_{0}^{t} e(\tau) d\tau$

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

# ╔═╡ 41fd4f98-9166-11eb-1f3c-776b51fec5fb
b = 1;

# ╔═╡ bf9dc0b4-9163-11eb-16fb-69102edbe01d
Kp = 10^3  # Proportional gain

# ╔═╡ bf84c80c-9163-11eb-16a9-25e697ae16fb
Ki = 0.1  # Integral gain

# ╔═╡ bf6ac3f8-9163-11eb-1717-e90a6dee7ea8
yref = 400.0 # Reference temperature

# ╔═╡ b110df24-9166-11eb-14c6-370b5230d543
u_in(err, int_err) = Kp * err + Ki * int_err # Input signal

# ╔═╡ bf4ed6b6-9163-11eb-0c54-995ad34d448e
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

# ╔═╡ dae6964a-9166-11eb-2b78-558153a8831a
N = 101 # Number of grid elements

# ╔═╡ de5782a6-9166-11eb-0aff-0f650914cc94
Δx = L/(N-1) # Finite discretization   

# ╔═╡ e1edca82-9166-11eb-03fe-27d2cf8194c9
# Diffusion matrix
M = spdiagm(-1 => ones(N-1), 0 => -2*ones(N), 1 => ones(N-1));

# ╔═╡ e5d0be3c-9166-11eb-21d4-6fc3fbd1c632
Matrix(M)[1:5,1:5]

# ╔═╡ f06ea6ec-9166-11eb-0185-3136fb3f8a5e
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


# ╔═╡ e4a4f5c2-9167-11eb-1630-f548b4278efa
# First row
M[1,2] = 2;

# ╔═╡ e48b5810-9167-11eb-2677-cf591bfc81af
# Last row
M[end,end-1] = 2;

# ╔═╡ e473a436-9167-11eb-308d-e15fa237bebf
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

$u(t) = K_{p} ~ e(t) + K_{i} \int_{0}^{t} e(\tau) d\tau = K_{p} ~ e(t) + K_{i} \varepsilon(t)$

with error 
$e(t) = y_{ref}(t) - y(t) = y_{ref}(t) - C ~ \theta(t)$

and integral error 

$\varepsilon(t) = \int_{0}^{t} e(\tau) d\tau.$

The thermal dynamics with PI control is described by

$\dot{\theta}(t) = A ~ \theta(t) + B ~ (K_{p} ~ e(t) + K_{i} \varepsilon(t)) + E ~ w(t)$

and

$\dot{\varepsilon}(t) = e(t) = y_{ref}(t) - C ~ \theta(t)$

which is further formulated as 



$\begin{matrix}
\dot{\theta}(t)= \\
\dot{\varepsilon}(t) = 
\end{matrix}
~
\begin{matrix}
(A - B~C~K_{p}) ~ \theta(t) + K_{p}~B~y_{ref} + K_{i}~B~\varepsilon(t)  + E ~ w(t) \\
-C ~ \theta(t) +  y_{ref}(t) 
\end{matrix}$

and finally noted in matrix-vector notation as 

$\begin{pmatrix}
\dot{\theta}(t)\\
\dot{\varepsilon}(t)
\end{pmatrix}
~=~
\begin{pmatrix}
A - K_{p}~B~C &  K_{i}~B \\
-C & 0
\end{pmatrix}
\begin{pmatrix}
\theta(t) \\
\varepsilon(t)
\end{pmatrix}
+
\begin{pmatrix}
K_{p}~B \\
1
\end{pmatrix}
~y_{ref}(t)
+
\begin{pmatrix}
E \\
0
\end{pmatrix}
w(t)$


This is an ordinary differential equation (ODE) that can be solved with common solvers like forward Euler method or Runge-Kutta scheme.

"

# ╔═╡ e455cb5a-9167-11eb-2c2f-d9c4f82d030a
A = α/(Δx^2) * M

# ╔═╡ e439e750-9167-11eb-0570-8593b291f149
begin
	B = spzeros(N); # Input vector
	B[1] = 2 * (α/(λ*Δx)) * b
	E = spzeros(N); # Disturbance vector
	E[end] = 2 * (α/(λ*Δx)) * b
	C = spzeros(1,N); # Output vector
	C[end] = 1;
end

# ╔═╡ e4203f8a-9167-11eb-36b4-f949479484c9
w(θ) = ϕout(θ) # 2 * (α/Δx) * vcat(zeros(N-1), [1]) * ϕout(θ)

# ╔═╡ 80ee2446-921f-11eb-23aa-7918b942b017
begin
	Acl = spzeros(N+1,N+1) # Closed-loop system matrix
	Acl[1:N,1:N] = A - Kp*B*C;
	Acl[1:N,end] = Ki * B
	Acl[end,1:N] = -C
end

# ╔═╡ d8be69e2-9222-11eb-33a8-6b21fe87aff0
eigvals(Matrix(Acl))

# ╔═╡ cd7823ce-921d-11eb-3917-934c7c0fe92d
begin
	Bcl = spzeros(N+1)  # Closed-loop input matrix
	Bcl[1:N] = Kp * B
	Bcl[end] = 1
end

# ╔═╡ 940f40f8-9221-11eb-1595-e32aa8b902d4
begin
	Ecl = spzeros(N+1)  # Closed-loop input matrix
	Ecl[1:N] = E
	Ecl[end] = 0
end

# ╔═╡ 14d10396-921f-11eb-3596-5bdd57220329
# Heat Equation as ODE
function heat_eq!(dz, z, p, t)

	θ = @views z[1:N]
	
	dz .= Acl * z + Bcl*yref + Ecl*ϕout(θ[end])
end

# ╔═╡ 6ffbe0e2-9229-11eb-25eb-9365d813a295
md"### Initial Conditions

For simplicity, the initial data of the heat equation is assumed as $\vartheta_{0}(x)= 300$ Kelvin for all $x \in \left[0,L\right]$. This function is approximated to gain the initial conditions of the ODE 

$\theta(0) = \left( 300, \cdots , 300 \right)^{\top}.$
"

# ╔═╡ 62a00876-9222-11eb-36cf-650b78a94456
θ₀ = 300.0 * ones(N)

# ╔═╡ fb4bdffa-921d-11eb-3c00-e3ffee61390b
md"## Simulation"

# ╔═╡ 261f295e-9222-11eb-01ac-61f68623a7f9
Δt = 10^(-2) # Sampling time

# ╔═╡ 5411bb50-921e-11eb-30e7-53a7b597db62
tspan = (0.0, Tf)

# ╔═╡ 75dc7136-9222-11eb-1634-cf5a3f3beb32
z0 = vcat(θ₀, [0])

# ╔═╡ 356616d2-921f-11eb-30d5-cb503c68d611
prob = ODEProblem( heat_eq!, z0, tspan ) # ODE Problem

# ╔═╡ 68de8dc8-921f-11eb-0a36-5dcc4f1cf721
sol = solve(prob,Euler(),dt=Δt,progress=true, saveat=1.0, save_start=true); # Solving the ODE

# ╔═╡ 09ba7556-9223-11eb-3776-71e66d4b66a3
# 1-dimensional grid
xspan = 0 : Δx : L

# ╔═╡ 94f2ce1f-6959-420f-a9e9-e17455ac70c1
size(sol)

# ╔═╡ 240643c2-9223-11eb-0ab3-b73a7cf9a826
sol[1:end,end]

# ╔═╡ c42bf0fe-9223-11eb-0eee-d3c931d98df8
heatmap(sol.t, xspan, sol[1:end-1,1:end], xaxis="Time [s]", yaxis="Position [m]", title="Evolution of temperature")

# ╔═╡ d8850b74-922a-11eb-3c09-4d93ace2e9a7
plot(sol.t,[sol[1,:], sol[end-1,:]], label=["Left" "Right"], title="Temperature at the left/right boundary", legend=:bottomright)

# ╔═╡ d8529aec-922a-11eb-2e06-17076ac56ebd


# ╔═╡ d83bc27a-922a-11eb-05fc-65acae3baea5


# ╔═╡ d821a52a-922a-11eb-2aa4-49b15b77e3ec


# ╔═╡ Cell order:
# ╟─af2dd4bc-9163-11eb-1991-ef8b7897519a
# ╠═c0607ff8-9163-11eb-1610-25eab3f919c4
# ╠═c040a860-9163-11eb-308a-332e7b19a186
# ╟─c021f60e-9163-11eb-34e3-6d8faa141b4a
# ╠═c006067e-9163-11eb-141d-23e5206e7e41
# ╠═bfec0b16-9163-11eb-0296-9db3fc9aaa79
# ╟─bfd345ea-9163-11eb-25e4-a7ed1ac077aa
# ╠═303f96d6-9166-11eb-3814-ef21a86a4972
# ╠═33f52982-9166-11eb-3de5-c1323f49271d
# ╠═6b4763f8-921c-11eb-1b6f-23505e52ed31
# ╟─bfb9cef8-9163-11eb-18d5-efa318c120d2
# ╠═41fd4f98-9166-11eb-1f3c-776b51fec5fb
# ╠═bf9dc0b4-9163-11eb-16fb-69102edbe01d
# ╠═bf84c80c-9163-11eb-16a9-25e697ae16fb
# ╠═bf6ac3f8-9163-11eb-1717-e90a6dee7ea8
# ╠═b110df24-9166-11eb-14c6-370b5230d543
# ╠═d8be69e2-9222-11eb-33a8-6b21fe87aff0
# ╟─bf4ed6b6-9163-11eb-0c54-995ad34d448e
# ╠═a986c796-9166-11eb-2a58-cb9aeb65d440
# ╠═dae6964a-9166-11eb-2b78-558153a8831a
# ╠═de5782a6-9166-11eb-0aff-0f650914cc94
# ╠═e1edca82-9166-11eb-03fe-27d2cf8194c9
# ╠═e5d0be3c-9166-11eb-21d4-6fc3fbd1c632
# ╟─f06ea6ec-9166-11eb-0185-3136fb3f8a5e
# ╠═e4a4f5c2-9167-11eb-1630-f548b4278efa
# ╠═e48b5810-9167-11eb-2677-cf591bfc81af
# ╟─e473a436-9167-11eb-308d-e15fa237bebf
# ╠═e455cb5a-9167-11eb-2c2f-d9c4f82d030a
# ╠═e439e750-9167-11eb-0570-8593b291f149
# ╠═e4203f8a-9167-11eb-36b4-f949479484c9
# ╠═80ee2446-921f-11eb-23aa-7918b942b017
# ╠═cd7823ce-921d-11eb-3917-934c7c0fe92d
# ╠═940f40f8-9221-11eb-1595-e32aa8b902d4
# ╠═14d10396-921f-11eb-3596-5bdd57220329
# ╟─6ffbe0e2-9229-11eb-25eb-9365d813a295
# ╠═62a00876-9222-11eb-36cf-650b78a94456
# ╠═fb4bdffa-921d-11eb-3c00-e3ffee61390b
# ╠═501b67a8-921e-11eb-15fe-33cebed6a427
# ╠═261f295e-9222-11eb-01ac-61f68623a7f9
# ╠═5411bb50-921e-11eb-30e7-53a7b597db62
# ╠═75dc7136-9222-11eb-1634-cf5a3f3beb32
# ╠═356616d2-921f-11eb-30d5-cb503c68d611
# ╠═68de8dc8-921f-11eb-0a36-5dcc4f1cf721
# ╠═09ba7556-9223-11eb-3776-71e66d4b66a3
# ╠═94f2ce1f-6959-420f-a9e9-e17455ac70c1
# ╠═240643c2-9223-11eb-0ab3-b73a7cf9a826
# ╠═c42bf0fe-9223-11eb-0eee-d3c931d98df8
# ╠═d8850b74-922a-11eb-3c09-4d93ace2e9a7
# ╠═d8529aec-922a-11eb-2e06-17076ac56ebd
# ╠═d83bc27a-922a-11eb-05fc-65acae3baea5
# ╠═d821a52a-922a-11eb-2aa4-49b15b77e3ec
