
"""
Definimos un objeto que contenga la información de las configuraciones caliente, fría
y el número N del tamaño de la red de spines
"""

mutable struct MicroEstados
    σ::Array{Int,1} #Arreglo de spines aleatorios
    N::Int          #Arreglo de N×N de spines
    frio::Array{Int,1} #Arreglo de spines alineados
end

"""
Red de spines caliente, donde se puede medir la proporción de spines "abajo" por
medio de f, siendo 1 todos abajo y 0 ningún abajo
"""

#Definimos una función que nos genere la red de spines con valores ∈{1,-1}
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
    
redSpinesFrio(N) = -ones(Int,N,N)





"""
La función devuelve un vector con las entradas de la matriz de espines

N := longitud del arreglo cuadrado de N×N
f := porcentaje de spines abajo de la red de spines 
"""
function Lista(N::Int,f::Float64)
    matriz = redSpinesCaliente(N::Int,f::Float64)
    lista = []
    for i in 1:N
        for j in 1:N
            push!(lista,matriz[i,j])
        end
    end
    
    frio = ones(N^2)
    
    MicroEstados(lista,N,frio)
end


"""
Generamos la lista de los primeros vecinos para un sitio (i,j) en la matriz y para
un sitio i del vector generado por la función Lista. La lista vendrá en la siguiente
configuración [arriba,derecha,abajo,izquierda]

Parámetros:

m := Microestado que guarda la información de la red 
     de spines de N×N y la longitud N.
i := sitio al que nos queremos enfocar

"""
function vecinos(m::MicroEstados,i::Int)
    N² = (m.N)^2
    N = m.N
    
    #=matriz = [1 2 3;4 5 6;7 8 9]
    N² = length(matriz)
    N = Int(sqrt(N²))
    lista = []
    for i in 1:N
        for j in 1:N
            push!(lista,matriz[i,j])
        end
    end
    lista=#                            #Para hacer pruebas de que esta bien
    
    
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
    return [m.σ[arriba],m.σ[derecha],m.σ[abajo],m.σ[izquierda]],
    [m.frio[arriba],m.frio[derecha],m.frio[abajo],m.frio[izquierda]]
    #[lista[arriba],lista[abajo],lista[derecha],lista[izquierda]]
end


"""
Ahora generamos la energía de interacción entre átomos dada por
E = -J⋅∑s_is_j

Parámetros:

m := Micro Estado con configuracion caliente y fria
i := sitio

"""
function energia_i(m::MicroEstados,i::Int)
    J = 0.5 #se escogió este valor por consenso en clase
    σi = m.σ[i]
    frio = m.frio[i]
    
    vecinosCaliente, vecinosFrio = vecinos(m,i)
    sumC = 0
    sumF = 0
    for i in 1:4
        sumC += vecinosCaliente[i]
        sumF += vecinosFrio[i]
    end
    return -2J*σi*sumC, -2J*frio*sumF
end   

function energiaTotal(m::MicroEstados)
    sumC = 0
    sumF = 0
    energias = []
    eCaliente = zeros(m.N^2)
    eFria = zeros(m.N^2)
    for i in 1:m.N^2
        push!(energias,energia_i(m,i))
        eCaliente[i] , eFria[i] = energias[i]
        #eCaliente[i] = energias[i]
        sumC += eCaliente[i]
        sumF += eFria[i]
    end
    
    return sumC,sumF #Por qué dividiamos entre dos?
end


#TEMPERATURA CRÍTICA

#Igualando la ecuación trascendental a cero y quedándonos con el lado izquierdo
tempCritica(T) = 2*tanh(1/T)^2 - 1 

function Dcentrada(f,a,h)
    derivada = (f(a+h) - f(a-h))/2h
    return derivada
end

function newtonRhapson(f,x_inicial,epsilon,n)
    h=0.0001
    if n==1
        return x_inicial
    elseif n > 1
        iteracion = newtonRhapson(f,x_inicial,epsilon,n-1)-f(newtonRhapson(f,x_inicial,epsilon,n-1))/Dcentrada(f,newtonRhapson(f,x_inicial,epsilon,n-1),h)
    end
    return iteracion
end


"""
Generamos una función para voltear el spín; se va a utilizar para
pasar de una configuración C_j a una C_k después de cierto número de pasos.

Parámetros

m := MicroEstado
i := Sitio de la configuración de la lista de spines generada.

"""
volteaSpinC(m::MicroEstados,i::Int) = m.σ[i] *= -1 #De esta manera se voltean los ±1
volteaSpinF(m::MicroEstados,i::Int) = m.frio[i] *= -1

"""
Generamos una función para cambiar de Micro Estado, si estamos en uno C_j y queremos
cambiar a otro C_k, proponemos ese cambio en lugar de pasar por un número n de
micro estados, esperando que de esta manera C_j y C_k no estén correlacionados.

Parámetros

m := MicroEstado


"""
function eleccionDeMicroEstado(m::MicroEstados)
    i = rand(1:m.N)
    ΔEcaliente, ΔEfria = energia_i(m,i)
    
    return ΔEcaliente, ΔEfria, i
end

"""
Generamos el paso Metropolis-MonteCarlo que consiste en ubicarse en un Micro Estado
arbitrario para voltearle el spin y calcular su energía de interacción con sus vecinos.
El micro estado se acepta si cumple que el mínimo entre e^(-β*ΔE) y 1 es mayor a un 
número aleatorio entre 0 y 1. Las pruebas en la simulación se harán para una temperatura
fría y otra caliente.

Parámteros

m := MicroEstado
T := Temperatura
conf := configuración caliente o fria

"""

function Metropolis(m::MicroEstados,T::Float64,conf::String)
    aceptado = false
    β = 1/T
    while aceptado == false
        
        ΔEcaliente , ΔEfria, i = eleccionDeMicroEstado(m)
        if conf == "caliente"
            parametro = min(1.,exp(-β*ΔEcaliente))
            if rand() < parametro
                aceptado = true
                volteaSpinC(m,i)
                return energia_i(m,i)[1]
            end
        
        elseif conf == "fria"
            parametro = min(1.,exp(-β*ΔEfria))
            if rand() < parametro
                aceptado = true
                volteaSpinF(m,i)
                return energia_i(m,i)[2]
            end
        else
            return "Escoge entre caliente o fria"
        end
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
    
    energiaCaliente = zeros(n)
    energiaFria = zeros(n)
    energiaCaliente[1] , energiaFria[1] = energiaTotal(m)
    

    if conf == "caliente"
        for i in 1:n-1 
            ΔE = Metropolis(m,β,conf) 
            energiaCaliente[i+1] = energiaCaliente[i] + ΔE
        end
        return energiaCaliente
   
    elseif conf == "fria"
        for i in 1:n-1
            ΔE = Metropolis(m,β,conf)
            energiaFria[i+1] = energiaFria[i] + ΔE
        end
        return energiaFria
    end
end




