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

function solve(A::CSRMatrix, f::Vector{<:Real}; solver::Symbol=:Jacobi, ω::Float64=1.95, ε::Float64=1.0e-3, max_iter::Int64=100)::Vector{<:Real}
    #TODO проверка на размерности
    if solver == :Jacobi
        return _solve_Jacobi(A, f, ε, max_iter)
    elseif solver == :Seidel
        return _solve_Seidel(A, f)
    elseif solver == :SOR
        return _solve_SOR(A, f, ω)
    else
        throw(ArgumentError("Неизвестный тип решателя \":$solver\". Доступные: :Jacobi, :Seidel, :SOR"))
    end
end

function _solve_Jacobi(A::CSRMatrix, f::Vector{<:Real}, ε::Float64, max_iter::Int64)::Vector{<:Real}
    _check_Jacobi(A)
    
    u_old = zeros(Float64, A.rows)
    u = zeros(Float64, A.rows)

    iter = 0
    while true
        for i in 1:A.rows
            ind1 = A.addres[i]
            ind2 = A.addres[i+1]-1
            d = Dict()
            for j in ind1:ind2
                d[A.columns[j]] = A.values[j]
            end
            s = 0.0
            for (col, val) in d
                s -= col != i ? val * u_old[col] : 0.0
            end
            s+=f[i]
            u[i] = s / get(d, i, 0.0)
        end
        
        if sum(abs.(u - u_old)) < ε
            break
        end
        u_old .= u
        
        iter += 1
        if iter >= max_iter
            @warn "Достигнуто максимальное число итераций $max_iter"
            break
        end
    end

    return u
end

function _check_Jacobi(A::CSRMatrix)::Nothing
    bad_rows = Vector{Int64}()
    for i in 1:A.rows
        ind1 = A.addres[i]
        ind2 = A.addres[i+1]-1
        d = Dict()
        for j in ind1:ind2
            d[A.columns[j]] = A.values[j]
        end
        s = 0.0
        for (col, val) in d
            s += col != i ? abs(val) : 0.0
        end
        if s > abs(get(d, i, 0.0))
            push!(bad_rows, i)
        end
    end

    if !isempty(bad_rows)
        @warn "Невыполнено достаточное условие сходимости с строках: $bad_rows"
    end
end

function _solve_Seidel(A::CSRMatrix, f::Vector{<:Real})::Vector{<:Real}
    _check_Seidel(A)
    
    u_old = zeros(Float64, A.rows)
    u = zeros(Float64, A.rows)

    iter = 0
    while true
        for i in 1:A.rows
            ind1 = A.addres[i]
            ind2 = A.addres[i+1]-1
            d = Dict()
            for j in ind1:ind2
                d[A.columns[j]] = A.values[j]
            end
            s = 0.0
            for (col, val) in d
                s -= col > i ? val * u[col] : 0.0
                s += col < i ? val * u[col] : 0.0
            end
            s+=f[i]
            u[i] = s / get(d, i, 0.0)
        end
        
        if sum(abs.(u - u_old)) < ε
            break
        end
        u_old .= u
        
        iter += 1
        if iter >= max_iter
            @warn "Достигнуто максимальное число итераций $max_iter"
            break
        end
    end

    return u
end

function _check_Seidel(A::CSRMatrix)::Nothing
    bad_rows = Vector{Int64}()
    for i in 1:A.rows
        ind1 = A.addres[i]
        ind2 = A.addres[i+1]-1
        d = Dict()
        for j in ind1:ind2
            d[A.columns[j]] = A.values[j]
        end
        if get(d, i, 0.0) < 0.0
            push!(bad_rows, i)
        end
    end

    if !isempty(bad_rows)
        @warn "Невыполнено достаточное условие сходимости с строках: $bad_rows"
    end
end

function _solve_SOR(A::CSRMatrix, f::Vector{<:Real}, ω::Float64)::Vector{<:Real}
    [1.0] #TODO
end

end #module MySparse

#TODO сборка матрицы;