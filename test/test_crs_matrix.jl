@from "../Modules/my_sparse.jl" using MySparse

@testset verbose = true "Умножение CRS матрицы на вектор" begin
    # {{1, 0, 0}, {0, 2, 0}, {5, 3, 4}}.{1, 2, 3} == {1, 4, 23}
    addres = [1, 2, 3, 6]
    values = [1.0, 2.0, 5.0, 3.0, 4.0]
    columns = [1, 2, 1, 2, 3]
    vector = [1, 2, 3]
    @test CRSMatrix(addres, columns, values) * vector == [1, 4, 23]

    # {{1.0, 0.0, 7.0, 0.0}, {2.0, 3.0, 0.0, 0.0}, {9.0, 0.0, 0.0, 6.0},{5.0, 4.0, 3.0, 1.0}}.{7.0, 0.0, 1.0, 2.0} == {14.0, 14.0, 75.0, 40.0}
    addres = [1, 3, 5, 7, 11]
    values = [1.0, 7.0, 2.0, 3.0, 9.0, 6.0, 5.0, 4.0, 3.0, 1.0]
    columns = [1, 3, 1, 2, 1, 4, 1, 2, 3, 4]
    vector = [7.0, 0, 1, 2]
    @test CRSMatrix(addres, columns, values) * vector == [14.0, 14.0, 75.0, 40.0]

    expected_message = "Размерности матрицы и вектора различны. Умножение невозможно"
    @test_throws ErrorException(expected_message) begin
        addres = [1, 2, 3, 6]
        values = [1.0, 2.0, 5.0, 3.0, 4.0]
        columns = [1, 2, 1, 2, 3]
        vector = [1, 2, 3, 4]
        CRSMatrix(addres, columns, values) * vector
    end

    expected_message = "Последний индекс вектора `addres` Должен равняться числу элементов вектора `values` + 1"
    @test_throws ArgumentError(expected_message) begin
        addres = [1, 2, 3, 8]
        values = [1.0, 2.0, 5.0, 3.0, 4.0]
        columns = [1, 2, 1, 2, 3]
        m = CRSMatrix(addres, columns, values)
    end

    expected_message = "Размерность вектора `columns` должна равняться размерности вектора `values`"
    @test_throws ArgumentError(expected_message) begin
        addres = [1, 2, 3, 6]
        values = [1.0, 2.0, 5.0, 3.0, 4.0]
        columns = [1, 2, 1]
        m = CRSMatrix(addres, columns, values)
    end

    expected_message = "`addres` должен быть возрастающим вектром положительных чисел"
    @test_throws ArgumentError(expected_message) begin
        addres = [-1, 2, 3, 6]
        values = [1.0, 2.0, 5.0, 3.0, 4.0]
        columns = [1, 2, 1, 2, 3]
        m = CRSMatrix(addres, columns, values)
    end

    expected_message = "`addres` должен быть возрастающим вектром положительных чисел"
    @test_throws ArgumentError(expected_message) begin
        addres = [-1, 6]
        values = [1.0, 2.0, 5.0, 3.0, 4.0]
        columns = [1, 2, 1, 2, 3]
        m = CRSMatrix(addres, columns, values)
    end

    expected_message = "`columns` должен быть положительным вектором"
    @test_throws ArgumentError(expected_message) begin
        addres = [1, 2, 3, 6]
        values = [1.0, 2.0, 5.0, 3.0, 4.0]
        columns = [-1, 2, 1, -2, 3]
        m = CRSMatrix(addres, columns, values)
    end

end
