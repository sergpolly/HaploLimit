gcc -I./include testhaplo.c -L./lib -Wl,-rpath -Wl,./lib -lglpk -std=c99 -o ../bin/testhaplo