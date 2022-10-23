CC=gcc -Wall -Wextra -Werror -std=c11
SRC=s21_matrix.c 
OBJ=$(SRC:%.c=%.o)
LIB=s21_matrix.a

all: 	$(LIB) test gcov_report

$(LIB):	$(OBJ)
		ar rc $(LIB) $^
		ranlib s21_matrix.a

%.o:	%.c s21_matrix.h
		$(CC) -c $<

test:
		$(CC) --coverage s21_test_matrix.c $(SRC) -lcheck -lm -lpthread -o test
		./test

gcov_report: 
		gcovr --html-details report.html

rebuild: clean all

clang:
	clang-format -n *.c

leaks:
	leaks -atExit -- ./test

clean :  
	rm -rf *.o test s21_matrix.a
	rm -rf *.c.gcov
	rm -rf *.gcno
	rm -rf *.gcda
	rm -rf *.css
	rm -rf *.functions.html
	rm -rf *.html
	rm -rf *.gcov