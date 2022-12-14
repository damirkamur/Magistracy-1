@from "../Modules/my_sparse.jl" using MySparse
@from "../Modules/tools.jl" using Tools

@testset verbose = true "Умножение CRS матрицы на вектор" begin
    # {{1, 0, 0}, {0, 2, 0}, {5, 3, 4}}.{1, 2, 3} == {1, 4, 23}
    addres = [1, 2, 3, 6]
    values = [1.0, 2.0, 5.0, 3.0, 4.0]
    columns = [1, 2, 1, 2, 3]
    vector = [1, 2, 3]
    @test CSRMatrix(addres, columns, values) * vector == [1, 4, 23]

    # {{1.0, 0.0, 7.0, 0.0}, {2.0, 3.0, 0.0, 0.0}, {9.0, 0.0, 0.0, 6.0},{5.0, 4.0, 3.0, 1.0}}.{7.0, 0.0, 1.0, 2.0} == {14.0, 14.0, 75.0, 40.0}
    addres = [1, 3, 5, 7, 11]
    values = [1.0, 7.0, 2.0, 3.0, 9.0, 6.0, 5.0, 4.0, 3.0, 1.0]
    columns = [1, 3, 1, 2, 1, 4, 1, 2, 3, 4]
    vector = [7.0, 0, 1, 2]
    @test CSRMatrix(addres, columns, values) * vector == [14.0, 14.0, 75.0, 40.0]

    # {{1.0, 0.0, 2.0}, {0.0, 0.0, 0.0}, {0.0, 4.0, 4.0}}.{1.0, 2.0, 3.0} == {7.0, 0.0, 20.0}
    addres = [1, 3, 3, 5]
    values = [1.0, 2.0, 4.0, 4.0]
    columns = [1, 3, 2, 3]
    vector = [1.0, 2.0, 3.0]
    @test CSRMatrix(addres, columns, values) * vector == [7.0, 0.0, 20.0]

    expected_message = "Размерности матрицы и вектора различны. Умножение невозможно"
    @test_throws ErrorException(expected_message) begin
        addres = [1, 2, 3, 6]
        values = [1.0, 2.0, 5.0, 3.0, 4.0]
        columns = [1, 2, 1, 2, 3]
        vector = [1, 2, 3, 4]
        CSRMatrix(addres, columns, values) * vector
    end

    expected_message = "Последний индекс вектора `addres` Должен равняться числу элементов вектора `values` + 1"
    @test_throws ArgumentError(expected_message) begin
        addres = [1, 2, 3, 8]
        values = [1.0, 2.0, 5.0, 3.0, 4.0]
        columns = [1, 2, 1, 2, 3]
        m = CSRMatrix(addres, columns, values)
    end

    expected_message = "Размерность вектора `columns` должна равняться размерности вектора `values`"
    @test_throws ArgumentError(expected_message) begin
        addres = [1, 2, 3, 6]
        values = [1.0, 2.0, 5.0, 3.0, 4.0]
        columns = [1, 2, 1]
        m = CSRMatrix(addres, columns, values)
    end

    expected_message = "`addres` должен быть возрастающим вектром положительных чисел"
    @test_throws ArgumentError(expected_message) begin
        addres = [-1, 2, 3, 6]
        values = [1.0, 2.0, 5.0, 3.0, 4.0]
        columns = [1, 2, 1, 2, 3]
        m = CSRMatrix(addres, columns, values)
    end

    expected_message = "`addres` должен быть возрастающим вектром положительных чисел"
    @test_throws ArgumentError(expected_message) begin
        addres = [-1, 6]
        values = [1.0, 2.0, 5.0, 3.0, 4.0]
        columns = [1, 2, 1, 2, 3]
        m = CSRMatrix(addres, columns, values)
    end

    expected_message = "`columns` должен быть положительным вектором"
    @test_throws ArgumentError(expected_message) begin
        addres = [1, 2, 3, 6]
        values = [1.0, 2.0, 5.0, 3.0, 4.0]
        columns = [-1, 2, 1, -2, 3]
        m = CSRMatrix(addres, columns, values)
    end

end

@testset verbose = true "Взятие элемента по индексам" begin
    # {{1, 0, 0}, {0, 0, 0}, {5, 3, 4}}
    addres = [1, 2, 2, 5]
    values = [1.0, 5.0, 3.0, 4.0]
    columns = [1, 1, 2, 3]
    m = CSRMatrix(addres, columns, values)
    @test m[1, 1] == 1.0
    @test m[1, 2] == 0.0
    @test m[1, 3] == 0.0
    @test m[2, 1] == 0.0
    @test m[2, 2] == 0.0
    @test m[2, 3] == 0.0
    @test m[3, 1] == 5.0
    @test m[3, 2] == 3.0
    @test m[3, 3] == 4.0

    expected_message = "Индекс строки выходит за пределы размерности матрицы"
    @test_throws ErrorException(expected_message) begin
        v = m[4, 3]
    end

    expected_message = "Индекс столбца выходит за пределы размерности матрицы"
    @test_throws ErrorException(expected_message) begin
        v = m[2, 0]
    end
end

@testset verbose = true "Сложение/Вычитание матриц" begin
    #{{1, 0}, {0, 1}} + {{0, 0}, {1, 0}}} == {{1, 0}, {1, 1}}
    addres = [1, 2, 3]
    values = [1.0, 1.0]
    columns = [1, 2]
    m1 = CSRMatrix(addres, columns, values)
    addres = [1, 1, 2]
    values = [1.0]
    columns = [1]
    m2 = CSRMatrix(addres, columns, values)
    addres = [1, 2, 4]
    values = [1.0, 1.0, 1.0]
    columns = [1, 1, 2]
    m3 = CSRMatrix(addres, columns, values)
    @test m1 + m2 == m3

    #{{7, 2, 15}, {0, 1, 0}, {3, 0, 81}, {0, 0, -1}} + {{1, 0, 0}, {0, 5, 0}, {4, 5, 0}, {0, 0, 0}} == {{8, 2, 15}, {0, 6, 0}, {7, 5, 81}, {0, 0, -1}}
    addres = [1, 4, 5, 7, 8]
    values = [7.0, 2.0, 15.0, 1.0, 3.0, 81.0, -1.0]
    columns = [1, 2, 3, 2, 1, 3, 3]
    m1 = CSRMatrix(addres, columns, values)
    addres = [1, 2, 3, 5, 5]
    values = [1.0, 5.0, 4.0, 5.0]
    columns = [1, 2, 1, 2]
    m2 = CSRMatrix(addres, columns, values)
    addres = [1, 4, 5, 8, 9]
    values = [8.0, 2.0, 15.0, 6.0, 7.0, 5.0, 81.0, -1.0]
    columns = [1, 2, 3, 2, 1, 2, 3, 3]
    m3 = CSRMatrix(addres, columns, values)
    @test m1 + m2 == m3

    #{{7, 2, 15}, {0, 1, 0}, {3, 0, 81}, {0, 0, -1}} - {{1, 0, 0}, {0, 5, 0}, {4, 5, 0}, {0, 0, 0}} == {{6, 2, 15}, {0, -4, 0}, {-1, -5, 81}, {0, 0, -1}}
    addres = [1, 4, 5, 7, 8]
    values = [7.0, 2.0, 15.0, 1.0, 3.0, 81.0, -1.0]
    columns = [1, 2, 3, 2, 1, 3, 3]
    m1 = CSRMatrix(addres, columns, values)
    addres = [1, 2, 3, 5, 5]
    values = [1.0, 5.0, 4.0, 5.0]
    columns = [1, 2, 1, 2]
    m2 = CSRMatrix(addres, columns, values)
    addres = [1, 4, 5, 8, 9]
    values = [6.0, 2.0, 15.0, -4.0, -1.0, -5.0, 81.0, -1.0]
    columns = [1, 2, 3, 2, 1, 2, 3, 3]
    m3 = CSRMatrix(addres, columns, values)
    @test m1 - m2 == m3

    expected_message = "Размеры матриц различны, сложение невозможно"
    @test_throws ErrorException(expected_message) begin
        addres = [1, 3, 4, 5, 5]
        values = [7.0, 2.0, 1.0, 3.0]
        columns = [1, 2, 2, 1]
        m1 = CSRMatrix(addres, columns, values)
        addres = [1, 2, 3, 5]
        values = [1.0, 5.0, 4.0, 5.0]
        columns = [1, 2, 1, 2]
        m2 = CSRMatrix(addres, columns, values)
        m3 = m1 + m2
    end

end

@testset verbose = true "Умножение/Деление матрицы на скаляр" begin
    #{{1, 0}, {0, 1}} * 2 = {{2, 0}, {0, 2}}
    addres = [1, 2, 3]
    values = [1.0, 1.0]
    columns = [1, 2]
    m1 = CSRMatrix(addres, columns, values)
    addres = [1, 2, 3]
    values = [2.0, 2.0]
    columns = [1, 2]
    m2 = CSRMatrix(addres, columns, values)
    number = 2.0
    @test m2 == m1 * number

    #{{6, 3, 15}, {0, 9, 0}, {3, 0, 81}} * 2/3 = {{4, 2, 10},{0, 6, 0},{2, 0, 54}}
    addres = [1, 4, 5, 7]
    values = [6.0, 3.0, 15.0, 9.0, 3.0, 81.0]
    columns = [1, 2, 3, 2, 1, 3]
    m1 = CSRMatrix(addres, columns, values)
    addres = [1, 4, 5, 7]
    values = [4.0, 2.0, 10.0, 6.0, 2.0, 54.0]
    columns = [1, 2, 3, 2, 1, 3]
    m2 = CSRMatrix(addres, columns, values)
    number = 2.0 / 3.0
    @test m2 == m1 * number

    #{{6, 4, 15}, {0, 9, 0}, {3, 0, 80}} / 2 = {{3, 2, 7.5},{0, 4.5, 0},{1.5, 0, 40}}
    addres = [1, 4, 5, 7]
    values = [6.0, 4.0, 15.0, 9.0, 3.0, 80.0]
    columns = [1, 2, 3, 2, 1, 3]
    m1 = CSRMatrix(addres, columns, values)
    addres = [1, 4, 5, 7]
    values = [3.0, 2.0, 7.5, 4.5, 1.5, 40.0]
    columns = [1, 2, 3, 2, 1, 3]
    m2 = CSRMatrix(addres, columns, values)
    number = 2.0
    @test m2 == m1 / number
end

@testset verbose = true "Решение системы уравнений" begin
    @testset "Ошибки" begin
        # {{0.0, 1.0, 2.0},{1.0, 0.0, 0.0}}.{..., ..., ...} == {1.0, 2.0, 3.0}
        addres = [1, 3, 4]
        values = [1.0, 2.0, 1.0]
        columns = [2, 3, 1]
        A = CSRMatrix(addres, columns, values)
        expected_message = "Матрица имеет некорректную размерность для расчета СЛАУ: 2x3"
        @test_throws ArgumentError(expected_message) begin
            solve(A, [1.0, 2.0, 3.0])
        end

        # {{0.0, 1.0},{1.0, 0.0}}.{..., ..., ...} == {1.0, 2.0, 3.0}
        addres = [1, 2, 3]
        values = [1.0, 1.0]
        columns = [2, 1]
        A = CSRMatrix(addres, columns, values)
        expected_message = "Размерность матрицы 2 не совпадает с размерностью правой части 3"
        @test_throws ArgumentError(expected_message) begin
            solve(A, [1.0, 2.0, 3.0])
        end

        expected_message = "Неизвестный тип решателя \":Solver\". Доступные: :Jacobi, :Seidel, :SOR"
        @test_throws ArgumentError(expected_message) begin
            addres = [1, 2, 3]
            values = [1.0, 1.0]
            columns = [1, 2]
            m1 = CSRMatrix(addres, columns, values)
            solve(m1, [1.0, 2.0]; solver = :Solver)
        end
    end

    @testset "Решения" begin
        #{{5.0, 0.0, 0.0}, {1.0, 2.0, 0.0}, {0.0, 0.0, 6.0}}.{1.0, 2.0, 3.0} == {5.0, 5.0, 18.0}
        addres = [1, 2, 4, 5]
        values = [5.0, 1.0, 2.0, 6.0]
        columns = [1, 1, 2, 3]
        A = CSRMatrix(addres, columns, values)
        res = [1.0, 2.0, 3.0]
        f = [5.0, 5.0, 18.0]
        @test L2_norm(A * solve(A, f) - f) < 10 - 5
        @test L2_norm(A * solve(A, f; solver = :Seidel) - f) < 10 - 5
        @test L2_norm(A * solve(A, f; solver = :SOR, norm = :L2) - f) < 10 - 5

        #{{3.0, 2.3, 0.0, 0.0}, {0.0, 3.5, 0.0, 0.0}, {0.0, 0.0, 8.0, 7.0}, {-1.0, 0.0, 3.0, 15.0}}.{2., 0., 9., -1.} = {6., 0., 65., 10.}
        addres = [1, 3, 4, 6, 9]
        values = [3.0, 2.3, 3.5, 8.0, 7.0, -1.0, 3.0, 15.0]
        columns = [1, 2, 2, 3, 4, 1, 3, 4]
        A = CSRMatrix(addres, columns, values)
        res = [2.0, 0.0, 9.0, -1.0]
        f = [6.0, 0.0, 65.0, 10.0]
        @test L2_norm(A * solve(A, f; ε = 1.0e-7) - f) < 10 - 5
        @test L2_norm(A * solve(A, f; solver = :Seidel, ε = 1.0e-7) - f) < 10 - 5
        @test L2_norm(A * solve(A, f; solver = :SOR, norm = :L2) - f) < 10 - 5

        # {{0.0, 2.0, 3.0}, {0.0, 0.0, 0.0}, {0.0, 6.0, 0.0}}.{1.0, 2.0, 3.0} == {13.0, 0.0, 12.0}
        addres = [1, 3, 3, 4]
        values = [2.0, 3.0, 6.0]
        columns = [2, 3, 2]
        A = CSRMatrix(addres, columns, values)
        res = [1.0, 2.0, 3.0]
        f = [13.0, 0.0, 12.0]
        @test_logs (:warn, "Невыполнено достаточное условие сходимости с строках: [1, 3]") (
            :warn,
            "Достигнуто максимальное число итераций 1000",
        ) solve(A, f)
        @test_logs (
            :warn,
            "Невыполнено достаточное условие сходимости с строках: [1, 2, 3]",
        ) (:warn, "Достигнуто максимальное число итераций 1000") solve(
            A,
            f;
            solver = :Seidel,
        )
    end
end

@testset verbose = true "Преобразование из удобного формата" begin
    # {{1, 0, 0}, {0, 2, 0}, {5, 3, 4}}
    addres = [1, 2, 3, 6]
    values = [1.0, 2.0, 5.0, 3.0, 4.0]
    columns = [1, 2, 1, 2, 3]
    row = [1, 2, 3, 3, 3]
    col = [1, 2, 1, 2, 3]
    val = [1.0, 2.0, 5.0, 3.0, 4.0]
    @test CSRMatrix(addres, columns, values) == CSRMatrix_from_RCV(row, col, val)

    # {{1.0, 0.0, 7.0, 0.0}, {2.0, 3.0, 0.0, 0.0}, {9.0, 0.0, 0.0, 6.0},{5.0, 4.0, 3.0, 1.0}}
    addres = [1, 3, 5, 7, 11]
    values = [1.0, 7.0, 2.0, 3.0, 9.0, 6.0, 5.0, 4.0, 3.0, 1.0]
    columns = [1, 3, 1, 2, 1, 4, 1, 2, 3, 4]
    row = [4, 4, 4, 4, 3, 3, 2, 2, 1, 1]
    col = [4, 3, 2, 1, 4, 1, 2, 1, 3, 1]
    val = [1.0, 3.0, 4.0, 5.0, 6.0, 9.0, 3.0, 2.0, 7.0, 1.0]
    @test CSRMatrix(addres, columns, values) == CSRMatrix_from_RCV(row, col, val)

    # {{1.0, 0.0, 2.0}, {0.0, 0.0, 0.0}, {0.0, 4.0, 4.0}}
    addres = [1, 3, 3, 5]
    values = [1.0, 2.0, 4.0, 4.0]
    columns = [1, 3, 2, 3]
    row = [1, 1, 3, 3]
    col = [3, 1, 3, 2]
    val = [2.0, 1.0, 4.0, 4.0]
    @test CSRMatrix(addres, columns, values) == CSRMatrix_from_RCV(row, col, val)

    # {{1.0, 0.0, 2.0}, {0.0, 0.0, 1.0}}
    addres = [1, 3, 4]
    values = [1.0, 2.0, 1.0]
    columns = [1, 3, 3]
    row = [1, 1, 2]
    col = [1, 3, 3]
    val = [1.0, 2.0, 1.0]
    @test CSRMatrix(addres, columns, values) == CSRMatrix_from_RCV(row, col, val)

    # {{0.0, 0.0, 0.0},{0.0, 1.0, 2.0}, {0.0, 3.0, 4.0}}
    addres = [1, 1, 3, 5]
    values = [1.0, 2.0, 3.0, 4.0]
    columns = [2, 3, 2, 3]
    row = [2, 2, 3, 3]
    col = [2, 3, 2, 3]
    val = [1.0, 2.0, 3.0, 4.0]
    @test CSRMatrix(addres, columns, values) == CSRMatrix_from_RCV(row, col, val)

    # {{1.0, 2.0},{3.0, 4.0},{0.0, 0.0}}}
    addres = [1, 3, 5, 5]
    values = [1.0, 2.0, 3.0, 4.0]
    columns = [1, 2, 1, 2]
    row = [1, 1, 2, 2]
    col = [1, 2, 1, 2]
    val = [1.0, 2.0, 3.0, 4.0]
    @test CSRMatrix(addres, columns, values) ==
          CSRMatrix_from_RCV(row, col, val; end_empty_rows = 1)

    # {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},{0.0, 1.0, 0.0},{0.0, 2.0, 0.0},{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}
    addres = [1, 1, 1, 2, 3, 3, 3]
    values = [1.0, 2.0]
    columns = [2, 2]
    row = [3, 4]
    col = [2, 2]
    val = [1.0, 2.0]
    @test CSRMatrix(addres, columns, values) ==
          CSRMatrix_from_RCV(row, col, val; end_empty_rows = 2)

    expected_message = "Размерности векторов row, col и val должны быть одинаковыми"
    @test_throws ArgumentError(expected_message) begin
        row = [3, 4, 6]
        col = [2, 2]
        val = [1.0, 2.0]
        m = CSRMatrix_from_RCV(row, col, val)
    end
end

@testset verbose = true "Разложение на L-D-U матрицы" begin
    # {{1, 0, 3}, {0, 0, 0}, {5, 3, 4}}
    addres = [1, 3, 3, 6]
    values = [1.0, 3.0, 5.0, 3.0, 4.0]
    columns = [1, 3, 1, 2, 3]
    m = CSRMatrix(addres, columns, values)
    (l1, d1, u1) = get_LDU(m)

    @testset verbose = true "Матрица L" begin
        l = get_L(m)
        @test l[1, 1] == 0.0
        @test l[1, 2] == 0.0
        @test l[2, 1] == 0.0
        @test l[2, 2] == 0.0
        @test l[3, 1] == 5.0
        @test l[3, 2] == 3.0
        @test l1[1, 1] == 0.0
        @test l1[1, 2] == 0.0
        @test l1[2, 1] == 0.0
        @test l1[2, 2] == 0.0
        @test l1[3, 1] == 5.0
        @test l1[3, 2] == 3.0
    end

    @testset verbose = true "Матрица D" begin
        d = get_D(m)
        @test d[1, 1] == 1.0
        @test d[1, 2] == 0.0
        @test d[1, 3] == 0.0
        @test d[2, 1] == 0.0
        @test d[2, 2] == 0.0
        @test d[2, 3] == 0.0
        @test d[3, 1] == 0.0
        @test d[3, 2] == 0.0
        @test d[3, 3] == 4.0
        @test d1[1, 1] == 1.0
        @test d1[1, 2] == 0.0
        @test d1[1, 3] == 0.0
        @test d1[2, 1] == 0.0
        @test d1[2, 2] == 0.0
        @test d1[2, 3] == 0.0
        @test d1[3, 1] == 0.0
        @test d1[3, 2] == 0.0
        @test d1[3, 3] == 4.0
    end

    @testset verbose = true "Матрица U" begin
        u = get_U(m)
        @test u[1, 1] == 0.0
        @test u[1, 2] == 0.0
        @test u[1, 3] == 3.0
        @test u[2, 1] == 0.0
        @test u[2, 2] == 0.0
        @test u[2, 3] == 0.0
        @test u[3, 1] == 0.0
        @test u[3, 2] == 0.0
        @test u[3, 3] == 0.0
        @test u1[1, 1] == 0.0
        @test u1[1, 2] == 0.0
        @test u1[1, 3] == 3.0
        @test u1[2, 1] == 0.0
        @test u1[2, 2] == 0.0
        @test u1[2, 3] == 0.0
        @test u1[3, 1] == 0.0
        @test u1[3, 2] == 0.0
        @test u1[3, 3] == 0.0
    end
end
