"""
Definimos un objeto que tenga en sus atributos las configuraciones fría y caliente
para poder trabajar con ellas durante el proceso. Las configuraciones son de 1D y 2D
"""

mutable struct MicroEstados
    N::Int          #Arreglo de N×N de spines
    σ::Array{Int,1} #Arreglo de spines aleatorios
    frio::Array{Int,1} #Arreglo de spines alineados
    redCaliente::Array{Int,2}
    redFria::Array{Int,2}
end

# RED DE SPINES CALIENTE Y FRÍA

"""
Red de spines caliente, donde se puede medir la proporción de spines "abajo" por
medio de f, siendo 1 todos abajo y 0 ningún abajo. Ej f=0.5 50% spines abajo

Parámetros:

N := tamaño de la red
f := flotante en [0,1]
"""

function redSpinesCaliente(N,f::Float64)
    red = ones(Int,N,N)
    for i in 1:N
        for j in 1:N
            if rand() <= f
                red[i,j] = -1
            end
        end
    end
    return red
end

"""Esta función genera una red de spines arriba, es decir, la configuración fría"""
redSpinesFrio(N) = ones(Int,N,N)

# DEFINICIÓN DE CONFIGURACIONES

"""
Configuraciones.

La función le da valores a las entradas del objeto MicroEstados. Esta función genera redes
de spines en su conf fría y caliente. Rellena los atributos del objeto MicroEstados con 
dos arreglos 1D con la configuración caliente y fría, dos arreglos 2D con la configuración
caliente y fría y N el tamaño de la malla.

N := longitud del arreglo cuadrado de N×N
f := porcentaje de spines abajo de la red de spines 
"""
function Configuraciones(N::Int,f::Float64)
    redCaliente = redSpinesCaliente(N::Int,f::Float64)
    lista = []
    for i in 1:N
        for j in 1:N
            push!(lista,redCaliente[i,j])
        end
    end
    redFria = redSpinesFrio(N)
    frio = ones(N^2)
    
    MicroEstados(N,lista,frio,redCaliente,redFria)
end

# ENERGÍA

"""
vecinos

Generamos la lista de los primeros vecinos para un sitio i del vector generado
por la función Configuraciones. La lista vendrá en la siguiente configuración 
[arriba,derecha,abajo,izquierda]

Parámetros:

m := Microestado que guarda la información de la red 
     de spines de N×N y la longitud N.
i := sitio al que nos queremos enfocar
conf := configuracion caliente o fria
"""
function vecinos(m::MicroEstados,i::Int,conf::String)
    N² = (m.N)^2
    N = m.N
    
    #Para el vecino de arriba
    if i <= N²
        arriba = i + N
        if arriba > N²
            arriba = i - N*(N-1)
        end
    end
    
    #Para el vecino de la derecha
    if mod(i,N) != 0
        derecha = i + 1
    else
        derecha = i - N + 1
    end
    
    #Para el vecino de abajo
    if i - N >= 1
        abajo = i - N
    else i - N < 1
        abajo = i + N*(N-1)
    end  
    
    #Para el vecino de la izquierda
    if mod(i-1,N) != 0
        izquierda = i - 1
    else
        izquierda = i + N - 1
    end

    #Regresa los índices, recuerda
    if conf == "caliente"
        return [m.σ[arriba],m.σ[derecha],m.σ[abajo],m.σ[izquierda]]
    elseif conf == "fria"
        [m.frio[arriba],m.frio[derecha],m.frio[abajo],m.frio[izquierda]]
    end
end


"""
enegia_i

Ahora generamos la energía de interacción entre átomos dada por
E = -J⋅∑s_is_j

Parámetros:

m := Micro Estado con configuracion caliente y fria
i := sitio
conf := configuracion fría o caliente

"""
function energia_i(m::MicroEstados,i::Int,conf::String)
    J = 0.5 #se escogió este valor por consenso en clase
    σi = m.σ[i]
    frio = m.frio[i]
    
    primerosVecinos = vecinos(m,i,conf)
    suma = sum(primerosVecinos)
    
    if conf == "caliente"
        return -2J*σi*suma
    elseif conf == "fria"
        return -2J*frio*suma
    end
end  

"""
energiaSitio

Calculamos la energía por cada sitio (i,j) de las configuraciones fría y caliente.

Parámetros:

m := MicroEstados
i := i-ésimo renglón
j := j-ésimo rengón
conf := configuración caliente o fría
"""

function energiaSitio(m::MicroEstados, i::Int, j::Int,conf::String)
    N = m.N
    J = 0.5
    if conf == "caliente"
        -2J*m.σ[i,j]*(m.σ[mod1(i-1,N),j] + m.σ[mod1(i+1,N),j] + m.σ[i,mod1(j-1,N)] + m.σ[i,mod1(j+1,N)])
    elseif conf == "fria"
        -2J*m.frio[i,j]*(m.frio[mod1(i-1,N),j]+m.frio[mod1(i+1,N),j]+m.frio[i,mod1(j-1,N)]+m.frio[i,mod1(j+1,N)])
    end
end

"""
energiaTotal

Energia total de toda la configuración, es decir la suma de todas las energías por sitio
del arreglo 1D

Parámetros: 

m := MicroEstados
conf := configuración caliente o fría
"""

function energiaTotal(m::MicroEstados,conf::String)
    sum = 0
    for i in 1:m.N^2
        sum += energia_i(m,i,conf)
    end
    return sum
end

"""
energiaTotalMalla

Energía total de toda la configuración, en este caso la suma de todas las energías por sitio
del arreglo matricial.

Parámetros:

m := MicroEstados
conf := configuración caliente o fría
"""

function energiaTotalMalla(m::MicroEstados,conf::String)
	out = 0.
    for i in 1:m.N, j in 1:m.N
		out += energia_ij(m,i,j,conf)
    end
    out
end

#TEMPERATURA CRÍTICA

#Igualando la ecuación trascendental a cero y quedándonos con el lado izquierdo
tempCritica(T) = 2*tanh(1/T)^2 - 1 

#Definimos diferencias finitas centradas para la derivada
function Dcentrada(f,a,h)
    derivada = (f(a+h) - f(a-h))/2h
    return derivada
end

#Definimos Newton Rhapson para resolver la ecuación trascendental
function newtonRhapson(f,x_inicial,epsilon,n)
    h=0.0001
    if n==1
        return x_inicial
    elseif n > 1
        iteracion = newtonRhapson(f,x_inicial,epsilon,n-1)-f(newtonRhapson(f,x_inicial,epsilon,n-1))/Dcentrada(f,newtonRhapson(f,x_inicial,epsilon,n-1),h)
    end
    return iteracion
end

# ALGORÍTMO METROPOLIS 1D

"""
Generamos una función para voltear el spín; se va a utilizar para
pasar de una configuración C_j a una C_k después de cierto número de pasos,
si es que la condición que impone Metropolis se cumple.

Parámetros

m := MicroEstado
i := Sitio de la configuración de la lista de spines generada.
j := j-ésimo renglón de la red de spines caliente/fría
k := k-ésima columna de la red de spines caliente/fría
conf := configuración fria o caliente.
"""
function volteaSpin(m::MicroEstados,i::Int,j::Int,k::Int,conf::String) 
    
    if conf == "caliente"
        m.redCaliente[j,k] *= -1
        m.σ[i] *= -1 #De esta manera se voltean los ±1
    elseif conf == "fria"
        m.redFria[j,k] *= -1
        m.frio[i] *= -1
    end
end

"""
Esta función es la misma solo que la ocupo en la función termalización.
"""

function volteaSpin1(m::MicroEstados,i::Int,conf::String) 
    
    if conf == "caliente"
        return m.σ[i] *= -1 #De esta manera se voltean los ±1
    elseif conf == "fria"
        return m.frio[i] *= -1
    end
end

"""
Generamos una función para cambiar de Micro Estado de manera aleatoria. Tomamos un 
sitio i, (j,k) (de ambos arreglos) aleatorio para actualizar el micro estado en el que
se hará el siguiente paso del algorítmo

Parámetros

m := MicroEstado
conf := configuracion caliente o fría
"""

function eleccionDeMicroEstado(m::MicroEstados,conf::String)
    i = rand(1:m.N^2)
    j , k = rand(1:m.N), rand(1:m.N)
    ΔE = -energia_i(m,i,conf)
    
    return ΔE, i, j, k
end

"""
Magnetización total es la suma de todos los spines de la configuración

Parámetros:

m    := MicroEstados
conf := configuración caliente o fría
"""

function magnetizacion(m::MicroEstados,conf::String)
    if conf == "caliente"
        return sum(m.σ)
    elseif conf == "fria"
        return sum(m.frio)
    end
end


"""
Generamos el paso Metropolis-MonteCarlo que consiste en ubicarse en un Micro Estado,
calcular su energía de interacción con sus vecinos y si cumple un par de condiciones
se le voltea el spín a ese sitio.

El micro estado se acepta si cumple que el mínimo entre e^(-β*ΔE) y 1 es mayor a un 
número aleatorio entre 0 y 1. Las pruebas en la simulación se harán para una temperatura
fría y otra caliente.

Parámteros

m := MicroEstado
T := Temperatura
conf := configuración caliente o fria

Comentarios: Esta función puede mejorarse, la idea de mejora es hace que Metropolis
funcione para un solo paso y ya que en simulaciónMontecarlo se haga el resto. Sospecho
que el while esta alentando la función. Y si lo vemos en perspectiva, ese while esta
ahí de adorno nomás. Porque el parámetro de aceptación está muy equis.
"""

function Metropolis(m::MicroEstados,β::Float64,conf::String)
    aceptado = false
    
    while aceptado == false
    #for i in 1:m.N^2
    
        ΔE, i, j, k = eleccionDeMicroEstado(m,conf)
    
        parametro = exp(-β*ΔE) #min(1.,exp(-β*ΔE))
        
        if ΔE < 0 || rand() < parametro #con < se viola la segunda ley de la termodinámica. No digas mamadas
            aceptado = true             #la configuración caliente se ordena a un ferromagneto.
            volteaSpin(m,i,j,k,conf)
        end
        E = energiaTotal(m,conf)
        M = magnetizacion(m,conf)
        return E , M
    end
end


"""
Simulación montecarlo, es utilizar metropolis para un número n de pasos.

Parámetros

m := MicroEstados
T := Temperatura
n := número de pasos
conf := configuración caliente o fria


"""
function simulacionMontecarlo(m::MicroEstados,T,n::Int,conf::String)
    β = 1/T

    energiaFinal = zeros(n)
    magnFinal = zeros(n)
        
    for i in 1:n
        ΔE , ΔM = Metropolis(m,β,conf) 
        energiaFinal[i] = ΔE
        magnFinal[i] = ΔM
    end
    return energiaFinal, magnFinal
end

"""
Es una función alternativa al cálculo de la simulación, voy a ver si mejor la optimizo en dos
partes como las anteriores a ver si funciona porque tarda bastante en compilarse para una 
red de 50: 900 segundos aproximadamente. 

Actualización 13abril: ya no pude optimizar la función. Me voy a quedar por el momento con las
anteriores porque para 100mil pasos me tarda dos minutos y eso es considerablemente menos que
los 900 segundos de termalización.

(Mismo día) Siempre si pude, esta función tarda lo mismo que las anteriores. Quedó solucionado el asunto y 
de hecho tarda aproximadamente lo mismo que la otra alternativa
"""


function termalizacion(m::MicroEstados,T,n::Int,conf::String)
    eFinal = zeros(n)
    mFinal = zeros(n)
    for k in 1:n 
        i = rand(1:m.N^2)
        dE = -energia_i(m,i,conf)
        parametro = exp(-dE/T)
        if dE < 0 || rand() < parametro
            volteaSpin1(m,i,conf)
        end
        ET = energiaTotal(m,conf)
        MT = magnetizacion(m,conf)
        eFinal[k] = ET
        mFinal[k] = MT
    end
    return eFinal, mFinal
end




# CANTIDADES OBSERVABLES

# En esta sección se van a determinar la energia, magnetización, calor específico y susceptibilidad magnética
# dado un conjunto de temperaturas. La idea es utilizar los algorítmos de metropolis para poder obtener 
# los resultados.

"""
Comenzamos con las primeras pruebas una vez que ya obtuvimos buenos resultados de ambas simulaciones montecarlo.
Vamos a comnezar por la energía y la magnetización. La idea es que dado un arreglo de temperaturas, se debe
obtener la energía y magnetización final de cada simulación...

cantidadesFinales de momento es menos eficiente que cantidadesObservables. Tarda considerablemente 
más en compilar sobre todo para cuando se trata de 100 000 iteraciones, se tarda varios minutos más.
Nos quedamos de momento con cantidadesObservables
"""


function cantidadesFinales(m::MicroEstados,T,n::Int,conf::String)
    eFinal = zeros(n)
    mFinal = zeros(n)
    for k in 1:n 
        i = rand(1:m.N^2)
        dE = -energia_i(m,i,conf)
        parametro = exp(-dE/T)
        if dE < 0 || rand() < parametro
            volteaSpin1(m,i,conf)
        end
    end
    return energiaTotal(m,conf), magnetizacion(m,conf)
end

function cantidadesObservables(m::MicroEstados,T,n::Int,conf::String)
    eFinal = zeros(length(T))
    mFinal = zeros(length(T))
    calorEsp = zeros(length(T))
    suscepMagn = zeros(length(T))
    medidas = 400
    for t in 1:length(T)
        E = zeros(medidas)
        M = zeros(medidas)
        for i in 1:medidas
            ener, magn = cantidadesFinales(m,T[t],n,conf)
            E[i] = ener
            M[i] = magn
        end
        eFinal[t] = mean(E)/m.N^2
        mFinal[t] = mean(M)/m.N^2
        calorEsp[t] = c_H(E,T[t])
        suscepMagn[t] = χ(M,T[t])
    end
    return eFinal, mFinal, calorEsp, suscepMagn
end       

"""
Nos quedamos con cantidades observables para trabajar. Hasta este punto, nos hace falta definir la
capacidad calorífica y la susceptibilidad magnética.
"""

c_H(energias,T::Float64) = (1/T^2)*var(energias)

χ(magn,T::Float64) = var(magn)/T


# ALGORÍTMO METROPOLIS 2D

#= function propone_cambio(m::MicroEstados, β::Float64,conf::String)
    i, j = rand(1:m.L), rand(1:m.L)  # Es más rápido que rand(1:m.L, 2)
    J = 0.5
	ΔE = -2J*energia_ij(m, i, j, conf)

	ΔE, i, j
end =#



function paso_montecarlo!(m::MicroEstados, β::Float64,conf::String)
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

function magnetizacion_total(m::MicroEstados,conf::String)
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


function microEstados_montecarlo(L::Int,f, T, num_pasos::Int,conf)
	β = 1/T
	m = Configuraciones(L,f)

	out = Array{Int,2}[copy(m.σ)]
	sizehint!(out, num_pasos)

	for i in 1:num_pasos-1
		paso_montecarlo!(m, β,conf)
		push!(out, copy(m.σ))
	end

	out
end
