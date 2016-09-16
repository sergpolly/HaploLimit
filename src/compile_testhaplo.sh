gcc -I./include testhaplo.c -L./lib -Wl,-rpath -Wl,./lib -lglpk -std=c99 -D_POSIX_SOURCE -o ../bin/testhaplo
