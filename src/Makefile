.PHONY: all s21_matrix_oop.a

CC=g++
FLAGS= -g -Wall -Werror -Wextra
S21_MATRIX_CPP=*.cpp 
S21_MATRIX_H=*.h
S21_MATRIX_O=*.o
TESTS=test

ifeq ($(shell uname -s),Linux)
	TEST_FLAGS+=-lsubunit -lgcov
endif

all: s21_matrix_oop.a

s21_matrix_oop.a:
	$(CC) $(S21_MATRIX_CPP) -c $(FLAGS)
	ar rc s21_matrix_oop.a $(S21_MATRIX_O)
	ranlib s21_matrix_oop.a

s21_covered:
	$(CC) $(S21_MATRIX_CPP) -c $(FLAGS) --coverage
	ar rc s21_matrix_oop.a $(S21_MATRIX_O)
	ranlib s21_matrix_oop.a

test: clean s21_covered
	$(CC) -o $(TESTS).o $(TESTS).cpp -c $(FLAGS) --coverage 
	$(CC) $(TESTS).o s21_matrix_oop.a -o $(TESTS) -lcheck --coverage $(FLAGS) -lgtest
	./$(TESTS)

clean:
	rm -rf $(S21_MATRIX_O) s21_matrix_oop.a $(TESTS).o $(TESTS) *.gcda *.gcno *.gcov *.info *.html report

gcov_report: test
	lcov -t "s21_matrix_oop" -o s21_matrix_oop.info -c -d .
	genhtml -o report s21_matrix_oop.info
