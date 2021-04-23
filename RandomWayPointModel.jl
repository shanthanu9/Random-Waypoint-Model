### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 412bc176-a3f5-11eb-0b3c-8bdd80245218
# Install required packages
begin
	import Pkg
	Pkg.add("Plots")
	Pkg.add("PlutoUI")
end

# ╔═╡ 4fea2195-e3e4-4281-a81a-04c12e4a7383
begin
	# Import required packages
	
	using Plots
	using PlutoUI
end

# ╔═╡ 89310cac-3f5a-41fe-b304-0d495cd7d714
md"""
# Random Waypoint Model
"""

# ╔═╡ 16754c7f-70e3-4176-8d4a-f6456ff98eeb
md"""
Let's first import the required packages.
"""

# ╔═╡ 90c90ff9-ff96-4de0-8565-c4b636e44a28
md"""
## Random Waypoint Model Implementation

Below is the implementation code for RWM.

Some of the implementations are hidden for brevity. They can be seen by clicking on the "eye" icon beside the cell in the notebook.
"""

# ╔═╡ f1ba21af-1a26-441a-baea-2940ad127012
# Global constansts

begin
	# Boundary of plot
	xlim = (-5, 5)
	ylim = (-5, 5)
end

# ╔═╡ eeffb367-1cf8-4465-b40e-6286e55604b1
# Function for generating required random numbers

begin
	function randx()
		print((xlim[2] - xlim[1])*rand() + xlim[1])
		
		return (xlim[2] - xlim[1])*rand() + xlim[1]
	end
	
	function randy()
		return (ylim[2] - ylim[1])*rand() + ylim[1]
	end
	
	function randv(V)
		return V*rand()
	end
end

# ╔═╡ 1604105c-0985-4fcb-9a77-3685f2623c2c
# Define "Node" parameters
begin
	Base.@kwdef mutable struct Node
		x::Float64
		y::Float64
		fx::Float64
		fy::Float64
		v::Float64
		dt::Float64
		p::Float64
	end

	function createNode(V)
		n = Node(
				randx(), 		# x position
				randy(),  		# y position
				randx(),     	# Final x position
				randy(), 		# Final y position
				randv(V),  		# Speed
				0.02, 			# Time granularity
				0 				# Pause time
			)
	end
end

# ╔═╡ c2a22336-f63d-4da2-9ca7-62316e7807c8
# "Step" function
begin
	function dist(x1, y1, x2, y2)
		return sqrt((y2-y1)^2 + (x2-x1)^2)
	end
	
	function step!(n::Node, V, P)
		# If node is paused
		if n.p > 0
			n.p -= 1
			return
		end
		
		D = dist(n.x, n.y, n.fx, n.fy)

		# If reached destination
		if n.v*n.dt >= D
			n.fx = randx()
			n.fy = randy()
			n.v = randv(V)
			n.p = P
			return
		end
		
		m = (n.fy - n.y)/(n.fx - n.x)	# Slope
		
		# Two choices for x and y
		x1 = n.x + n.v*n.dt/sqrt(1+m^2)
		x2 = n.x - n.v*n.dt/sqrt(1+m^2)
		
		y1 = m*(x1-n.x) + n.y
		y2 = m*(x2-n.x) + n.y
		
		# Choose the one at least distance from destination
		if dist(x1, y1, n.fx, n.fy) < dist(x2, y2, n.fx, n.fy)
			n.x = x1
			n.y = y1
		else
			n.x = x2
			n.y = y2
		end
	end
end

# ╔═╡ 80a7a887-afb6-4e0e-95ff-43a8eb1c7d96
"""
    RandomWayPointModel(n_nodes, V, P[, trace])

Output the simulation of Random Waypoint model as gif.

`n_nodes` is the number of nodes in RWM simulation.
`V` is the maximum velocity of each node.
`P` is the pause time of each node.
If `trace` is true, then the trail of the node's movement is drawn in the plot.

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
function RandomWayPointModel(n_nodes, V, P, epochs, trace=false)
	
	nodes = []
	for i=1:n_nodes
		push!(nodes, createNode(V))
	end
	
	if trace == true
		plt = plot(size(nodes)[1], xlim=xlim, ylim=ylim)
		push!(plt, [n.x for n in nodes], [n.y for n in nodes])
		
		@gif for epoch=1:epochs
			for n in nodes
				step!(n, V, P)
			end
			push!(plt, [n.x for n in nodes], [n.y for n in nodes])
		end
	else
		scatter([n.x for n in nodes], [n.y for n in nodes], xlim=xlim, ylim=ylim)
	
		@gif for epoch=1:epochs
			for n in nodes
				step!(n, V, P)
			end
			scatter([n.x for n in nodes], [n.y for n in nodes], xlim=xlim, ylim=ylim)
		end
	end
end

# ╔═╡ ab3a5884-cf45-48e1-9dd9-47061ce0c0a3
md"""
The above implementation outputs gifs of RWM simulation based on the given input parameters
"""

# ╔═╡ 6b95209e-2a1f-4423-b296-97ea0ee87543
md"""
#### Single node RWM with trace
"""

# ╔═╡ 35f3e884-fb42-4b4d-866e-ae2f3401aeb5
RandomWayPointModel(1, 50, 10, 200, true)

# ╔═╡ 4cc55d97-8ab1-4ab8-8831-f5d8f2e9a3a9
md"""
#### Single node RWM without trail
"""

# ╔═╡ fd0179a7-6a3b-4731-8f2f-42b68fff136b
RandomWayPointModel(1, 50, 10, 200)

# ╔═╡ d2231069-36e1-4a41-996c-f1143ec413a4
md"""
#### Stable RWM with small V and large P

V = 1, P = 10

V is maximum velocity of each node and P is pause time at each transition.
"""

# ╔═╡ 82dbbcec-3903-411f-9de9-61996a454317
RandomWayPointModel(50, 1, 10, 200)

# ╔═╡ 504bfce2-58c2-4da2-8d0f-ce65ffeafd79
md"""
#### Dynmaic RWM with large V and small P

V = 100, P = 1
"""

# ╔═╡ 9b6ffb89-7c26-4de0-b48b-54775412ead5
RandomWayPointModel(50, 100, 1, 200)

# ╔═╡ c3ad8413-9b65-4edf-8ea6-4166f0a88196
md"""
#### Stable RWM with large V and large P

V = 100, P = 75
"""

# ╔═╡ 47d5372b-0a98-4eb7-977b-bb36505c20b9
RandomWayPointModel(50, 100, 75, 1000)

# ╔═╡ be370624-43b1-4dfa-968c-d99d1755d4ee
md"""
#### RWM with 100 nodes

V = 50, P = 10

Observe that over time, there are lot more nodes in the center than in the boundaries.
"""

# ╔═╡ ede983a0-fb84-4501-9274-b4e7c088058a
RandomWayPointModel(100, 50, 10, 1000)

# ╔═╡ e2596a00-e6a4-47f0-bc91-ce8f815d7f42
RandomWayPointModel(0, 1, 0, 1)

# ╔═╡ Cell order:
# ╟─89310cac-3f5a-41fe-b304-0d495cd7d714
# ╟─16754c7f-70e3-4176-8d4a-f6456ff98eeb
# ╠═412bc176-a3f5-11eb-0b3c-8bdd80245218
# ╠═4fea2195-e3e4-4281-a81a-04c12e4a7383
# ╟─90c90ff9-ff96-4de0-8565-c4b636e44a28
# ╠═f1ba21af-1a26-441a-baea-2940ad127012
# ╟─eeffb367-1cf8-4465-b40e-6286e55604b1
# ╟─1604105c-0985-4fcb-9a77-3685f2623c2c
# ╟─c2a22336-f63d-4da2-9ca7-62316e7807c8
# ╠═80a7a887-afb6-4e0e-95ff-43a8eb1c7d96
# ╟─ab3a5884-cf45-48e1-9dd9-47061ce0c0a3
# ╟─6b95209e-2a1f-4423-b296-97ea0ee87543
# ╠═35f3e884-fb42-4b4d-866e-ae2f3401aeb5
# ╟─4cc55d97-8ab1-4ab8-8831-f5d8f2e9a3a9
# ╠═fd0179a7-6a3b-4731-8f2f-42b68fff136b
# ╟─d2231069-36e1-4a41-996c-f1143ec413a4
# ╠═82dbbcec-3903-411f-9de9-61996a454317
# ╟─504bfce2-58c2-4da2-8d0f-ce65ffeafd79
# ╠═9b6ffb89-7c26-4de0-b48b-54775412ead5
# ╟─c3ad8413-9b65-4edf-8ea6-4166f0a88196
# ╠═47d5372b-0a98-4eb7-977b-bb36505c20b9
# ╟─be370624-43b1-4dfa-968c-d99d1755d4ee
# ╠═ede983a0-fb84-4501-9274-b4e7c088058a
# ╠═e2596a00-e6a4-47f0-bc91-ce8f815d7f42
