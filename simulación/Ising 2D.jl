struct MicroEstado
    σ::Array{Int,2}
	# Vamos a suponer que todas las configuraciones son cuadradas
    L::Int
    frio::Array{Int,2}
end

function edo_aleatorio(L::Int)
    σ = ones(Int, (L,L))
    for i in 1:L^2
        if rand() <= 0.5
			σ[i] = -1
        end
	end
    frio = ones(Int,L,L)
    MicroEstado(σ,L,frio)
end

function voltea_espin!(m::MicroEstado, i::Int, j::Int,conf::String)
    if conf == "caliente"
        m.σ[i,j] *= -1
    elseif conf == "fria"
        m.frio[i,j] *= -1
    end
end

function energia_ij(m::MicroEstado, i::Int, j::Int,conf::String)
    L = m.L
	J = 0.5
    if conf == "caliente"
        -2J*m.σ[i,j]*(m.σ[mod1(i-1,L),j] + m.σ[mod1(i+1,L),j] + m.σ[i,mod1(j-1,L)] + m.σ[i,mod1(j+1,L)])
    elseif conf == "fria"
        -2J*m.frio[i,j]*(m.frio[mod1(i-1,L),j]+m.frio[mod1(i+1,L),j]+m.frio[i,mod1(j-1,L)]+m.frio[i,mod1(j+1,L)])
        #-m.frio[i,j]*(m.frio[mod1(i-1,L),j]+m.frio[mod1(i+1,L),j]+m.frio[i,mod1(j-1,L)]+m.frio[i,mod1(j+1,L)])
    end
end
        

function energia_total(m::MicroEstado,conf::String)
	out = 0.
    for i in 1:m.L, j in 1:m.L
		out += energia_ij(m,i,j,conf)
    end
    out
end

function propone_cambio(m::MicroEstado, β::Float64,conf::String)
    i, j = rand(1:m.L), rand(1:m.L)  # Es más rápido que rand(1:m.L, 2)
	ΔE = -energia_ij(m, i, j, conf)

	ΔE, i, j
end

function paso_montecarlo!(m::MicroEstado, β::Float64,conf::String)
	aceptado = false

	while aceptado == false
		ΔE, i, j = propone_cambio(m, β,conf)

		# El parámetro de aceptación
		α = min(1., exp(-β*ΔE))

		if rand() < α
			aceptado = true
			ΔM = -2*m.σ[i,j]
			voltea_espin!(m, i, j,conf)
			return ΔE, ΔM
		end
    end
end

function magnetizacion_total(m::MicroEstado,conf::String)
    if conf == "caliente"
        sum(m.σ)
    elseif conf == "fria"
        sum(m.frio)
    end
end

function simulacion_montecarlo(L::Int, T, num_pasos::Int,conf)
	β = 1/T
	m = edo_aleatorio(L)

	ener = zeros(num_pasos)
	ener[1] = energia_total(m,conf)
	mag = zeros(num_pasos)
	mag[1] = magnetizacion_total(m,conf)

	for i in 1:num_pasos-1
		ΔE, ΔM = paso_montecarlo!(m, β,conf)
		ener[i+1] = ener[i] + ΔE
		mag[i+1] = mag[i] + ΔM
	end

	return ener, mag, m
end

function microEstados_montecarlo(L::Int, T, num_pasos::Int,conf)
	β = 1/T
	m = edo_aleatorio(L)

	out = Array{Int,2}[copy(m.σ)]
	sizehint!(out, num_pasos)

	for i in 1:num_pasos-1
		paso_montecarlo!(m, β,conf)
		push!(out, copy(m.σ))
	end

	out
end
