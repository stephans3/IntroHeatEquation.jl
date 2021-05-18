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

# ╔═╡ 71903690-5f0e-11eb-2435-a56a85d76fe0
using Plots

# ╔═╡ 7eef98bc-5f0e-11eb-0dc2-43a0e5715aa4
using PlutoUI

# ╔═╡ d5cced56-5f0c-11eb-36ce-37a624f02490
md"# Heat kernel

The heat kernel is defined by

$K(t,x,y) = \frac{1}{(4 ~ \pi ~ t)^{d/2}} ~ \exp\left(\frac{-\lvert x - y \rvert^2}{4 t}\right)$

with $d = 1,2,\cdots$ as the number of dimensions, time $t \geq 0$ and space $x \in \mathbb{R}^{d}$.
"

# ╔═╡ 4131c77a-5f0e-11eb-1aae-41015705d66e
function kernel(t, x, y; dim = 1)
    d = dim

    arg = -1 * (abs(x - y)^2)/(4*t)
    den = (4*pi*t)^(d/2)
    
    sol = (1/den) * exp(arg)

    return sol
end

# ╔═╡ 4118304e-5f0e-11eb-3844-5159c87d48b7
tspan = 0.001 : 0.001 : 1.0 # time

# ╔═╡ 40fa23ec-5f0e-11eb-2b02-73a75b7a6049
xspan = -3.0 : 0.1 : 3.0 # space

# ╔═╡ 401f2062-5f0e-11eb-2517-5900dd4b84ba
y = 0.0 # Offset

# ╔═╡ 679b9134-5f0e-11eb-159c-258752c5fcef
sol = zeros(length(xspan), length(tspan));

# ╔═╡ 6f5614c6-5f0e-11eb-3412-5381fbf426a0
for i = 1 : length(tspan)
    for k = 1 : length(xspan)
        sol[k,i] = kernel(tspan[i], xspan[k], y)
    end
end

# ╔═╡ ab00a9dc-5f0e-11eb-36e1-fd6e0d1acb68
@bind pos Slider(1:10:1000, default=1)

# ╔═╡ be7a3cf8-5f0e-11eb-3bb6-0d132e1b8e42
scatter(xspan, sol[:,pos], legend=false)

# ╔═╡ Cell order:
# ╟─d5cced56-5f0c-11eb-36ce-37a624f02490
# ╠═4131c77a-5f0e-11eb-1aae-41015705d66e
# ╠═4118304e-5f0e-11eb-3844-5159c87d48b7
# ╠═40fa23ec-5f0e-11eb-2b02-73a75b7a6049
# ╠═401f2062-5f0e-11eb-2517-5900dd4b84ba
# ╠═679b9134-5f0e-11eb-159c-258752c5fcef
# ╠═6f5614c6-5f0e-11eb-3412-5381fbf426a0
# ╠═71903690-5f0e-11eb-2435-a56a85d76fe0
# ╠═7eef98bc-5f0e-11eb-0dc2-43a0e5715aa4
# ╠═ab00a9dc-5f0e-11eb-36e1-fd6e0d1acb68
# ╠═be7a3cf8-5f0e-11eb-3bb6-0d132e1b8e42
