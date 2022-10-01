module MySparse

export CRSMatrix

mutable struct CRSMatrix
    addres::Vector{Int64}
    columns::Vector{Int64}
    values::Vector{Float64}

    function CRSMatrix(addres::Vector{Int64}, columns::Vector{Int64}, values::Vector{Float64})
        addres[end] != length(values) + 1 && throw(ArgumentError("Последний индекс вектора `addres` Должен равняться числу элементов вектора `values` + 1"))
        length(columns) != length(values) && throw(ArgumentError("Размерность вектора `columns` должна равняться размерности вектора `values`"))
        for i in 1:length(addres)-1
            addres[i] >= addres[i+1] || addres[i] < 0 && throw(ArgumentError("`addres` должен быть возрастающим вектром положительных чисел"))
        end
        all((x->x>0).(columns)) || throw(ArgumentError("`columns` должен быть положительным вектором"))
        new(addres, columns, values)    
    end
end

function Base.show(io::IO, c::CRSMatrix)
    println(io, "CRSMatrix:")
    println(io, "└addres: $(c.addres)")
    println(io, "└columns: $(c.columns)")
    println(io, "└values: $(c.values)")
end

function Base.:*(c::CRSMatrix, vector::Vector)
    length(vector) == max(c.columns...) || throw(error("Размерности матрицы и вектора различны. Умножение невозможно"))
    result = Float64[]
    for i in 1:length(c.addres)-1
        push!(result, sum([c.values[j] * vector[c.columns[j]] for j in c.addres[i]:c.addres[i+1]-1]))
    end
    return result
end

end #module MySparse

#TODO сборка матрицы; доспуп к элементам по идексам; сложение матриц; 