# gcc -Ofast main_threads.cpp -lm -ljson11 -lstdc++ -pthread -lprofiler
# export OMP_NUM_THREADS=6
# gcc -I/usr/local/include -pthread -Wl,-rpath -Wl,/usr/local/lib -Wl,--enable-new-dtags -L/usr/local/lib -lmpi -Ofast main.cpp -lstdc++ -lprofiler -lm -ljson11 -fopenmp
gcc -g3 -march=native -Ofast -std=c++17 -I/usr/local/include -I/usr/include/eigen3 -pthread -lmpi -mfma -funroll-loops -mfpmath=sse -msse2 main_threads_omp.cpp -Wl,--start-group -lstdc++ -lprofiler -lm -ljson11 -fopenmp -lglog -lgflags -ldl -ldouble-conversion -lfolly -Wl,--end-group