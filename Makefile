.PHONY: format, tests, docs

format:
	@echo "Форматирование проекта"
	julia -e "using JuliaFormatter;format(\"src\");format(\"test\")"

tests:
	@echo "Запуск тестирования проекта"
	julia --project test/runtest.jl
