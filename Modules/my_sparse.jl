module MySparse

export CSRMatrix, solve

mutable struct CSRMatrix
    addres::Vector{Int64}
    columns::Vector{Int64}
    values::Vector{Float64}
    rows::Int64
    cols::Int64

    function CSRMatrix(addres::Vector{Int64}, columns::Vector{Int64}, values::Vector{Float64})
        addres[end] != length(values) + 1 && throw(ArgumentError("Последний индекс вектора `addres` Должен равняться числу элементов вектора `values` + 1"))
        length(columns) != length(values) && throw(ArgumentError("Размерность вектора `columns` должна равняться размерности вектора `values`"))
        for i in 1:length(addres)-1
            addres[i] > addres[i+1] || addres[i] < 0 && throw(ArgumentError("`addres` должен быть возрастающим вектром положительных чисел"))
        end
        all((x->x>0).(columns)) || throw(ArgumentError("`columns` должен быть положительным вектором"))
        rows = length(addres) - 1
        cols = max(columns...)
        new(addres, columns, values, rows, cols)    
    end
end

function Base.show(io::IO, c::CSRMatrix)
    println(io, "CSRMatrix:")
    println(io, "└addres: $(c.addres)")
    println(io, "└columns: $(c.columns)")
    println(io, "└values: $(c.values)")
    println(io, "└rows: $(c.rows)")
    println(io, "└cols: $(c.cols)")
end

function Base.:*(c::CSRMatrix, vector::Vector{<:Real})::Vector{<:Real}
    length(vector) == c.cols || throw(error("Размерности матрицы и вектора различны. Умножение невозможно"))
    result = Float64[]
    for i in 1:length(c.addres)-1
        push!(result, sum([c.values[j] * vector[c.columns[j]] for j in c.addres[i]:c.addres[i+1]-1]))
    end
    return result
end


function Base.getindex(c::CSRMatrix, i::Int64, j::Int64)::Real
    1 ≤ i ≤ c.rows || throw(error("Индекс строки выходит за пределы размерности матрицы"))
    1 ≤ j ≤ c.cols || throw(error("Индекс столбца выходит за пределы размерности матрицы"))
    (ind1, ind2) = c.addres[i], c.addres[i+1]

    res = 0.0
    for ind in ind1:(ind2 - 1)
        if j == c.columns[ind]
            res = c.values[ind]
        end
    end

    return res
end

function Base.:+(m1::CSRMatrix, m2::CSRMatrix)::CSRMatrix
    m1.cols != m2.cols  || m1.rows != m2.rows && throw(error("Размеры матриц различны, сложение невозможно"))

    addres = [1]
    columns = Int64[]
    values = Float64[]

    iter1 = 1
    iter2 = 1
    for i in 1:m1.rows
        count1 = m1.addres[i+1]-m1.addres[i]
        count2 = m2.addres[i+1]-m2.addres[i]
        d = Dict()
        for _ in 1:count1
            d[m1.columns[iter1]] = m1.values[iter1]
            iter1 += 1
        end
        for _ in 1:count2
            d[m2.columns[iter2]] = get(d, m2.columns[iter2], 0.0) + m2.values[iter2]
            iter2 += 1
        end
        push!(addres, addres[end] + length(d))
        for (col, val) in sort(collect(d), by=x->x[1])
            push!(columns, col)
            push!(values, val)
        end
    end
    return CSRMatrix(addres, columns, values)
end

Base.:-(m1::CSRMatrix, m2::CSRMatrix)::CSRMatrix = m1 + m2 * (-1)

function Base.:(==)(m1::CSRMatrix, m2::CSRMatrix)::Bool
    m1.addres != m2.addres && return false 
    m1.columns != m2.columns && return false 
    m1.values != m2.values && return false 
    m1.rows != m2.rows && return false 
    m1.cols != m2.cols && return false 
    return true
end

function Base.:*(c::CSRMatrix, number::Real)::CSRMatrix
    return CSRMatrix(c.addres, c.columns, c.values*number)
end

function solve(A::CSRMatrix, f::Vector{<:Real}; solver::Symbol=:Jacobi, ω::Float64=1.95)::Vector{<:Real}
    if solver == :Jacobi
        return _solve_Jacobi(A, f)
    elseif solver == :Seidel
        return _solve_Seidel(A, f)
    elseif solver == :SOR
        return _solve_SOR(A, f, ω)
    else
        throw(ArgumentError("Неизвестный тип решателя \":$solver\". Доступные: :Jacobi, :Seidel, :SOR"))
    end
end

function _solve_Jacobi(A::CSRMatrix, f::Vector{<:Real})::Vector{<:Real}
    [1.0] #TODO
end

function _solve_Seidel(A::CSRMatrix, f::Vector{<:Real})::Vector{<:Real}
    [1.0] #TODO
end

function _solve_SOR(A::CSRMatrix, f::Vector{<:Real}, ω::Float64)::Vector{<:Real}
    [1.0] #TODO
end

end #module MySparse

#TODO сборка матрицы;