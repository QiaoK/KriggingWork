CC=gcc
CFLAGS=-c -Wall -Wextra -O1
LIBS = -lm
CORE_COMPONENTS = krig_functions.o cluster_functions.o krig_cluster.o matrix_functions.o data_functions.o distributions.o random.o
IGRA_TEST_OBJS = IGRA_test.o $(CORE_COMPONENTS)
MATRIX_TEST_OBJS = matrix_test.o matrix_functions.o
SOCR_TEST_OBJS = SOCR_test.o $(CORE_COMPONENTS)
KRIG_TEST_OBJS = krig_test.o krig_functions.o cluster_functions.o matrix_functions.o
GENERATED_TEST_OBJS = generated_test.o spatial_temporal_generator.o $(CORE_COMPONENTS)
SOCR_REGRESSION_TEST_OBJS = SOCR_regression_test.o regression.o $(CORE_COMPONENTS)
VARIOGRAM_TEST_OBJS = variogram_test.o variogram_training.o $(CORE_COMPONENTS)
CORE_POINT_TEST_OBJS = core_point_test.o $(CORE_COMPONENTS)
IGRA_VARIOGRAM_TEST_OBJS = IGRA_variogram_test.o variogram_training.o $(CORE_COMPONENTS)

All: matrix_test IGRA_test SOCR_test krig_test SOCR_regression_test variogram_test
IGRA_test : $(IGRA_TEST_OBJS)
	$(CC) -o $@ $(IGRA_TEST_OBJS) $(LIBS)
matrix_test : $(MATRIX_TEST_OBJS)
	$(CC) -o $@ $(MATRIX_TEST_OBJS) $(LIBS)
SOCR_test : $(SOCR_TEST_OBJS)
	$(CC) -o $@ $(SOCR_TEST_OBJS) $(LIBS)
krig_test : $(KRIG_TEST_OBJS)
	$(CC) -o $@ $(KRIG_TEST_OBJS) $(LIBS)
generated_test : $(GENERATED_TEST_OBJS)
	$(CC) -o $@ $(GENERATED_TEST_OBJS) $(LIBS)
SOCR_regression_test : $(SOCR_REGRESSION_TEST_OBJS)
	$(CC) -o $@ $(SOCR_REGRESSION_TEST_OBJS) $(LIBS)
variogram_test : $(VARIOGRAM_TEST_OBJS)
	$(CC) -o $@ $(VARIOGRAM_TEST_OBJS) $(LIBS)
core_point_test : $(CORE_POINT_TEST_OBJS)
	$(CC) -o $@ $(CORE_POINT_TEST_OBJS) $(LIBS)
IGRA_variogram_test : $(IGRA_VARIOGRAM_TEST_OBJS)
	$(CC) -o $@ $(IGRA_VARIOGRAM_TEST_OBJS) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<  
clean:
	rm -rf *.o
	rm -rf cluster_test
	rm -rf matrix_test
	rm -rf IGRA_test
	rm -rf krig_test
	rm -rf SOCR_test
	rm -rf generated_test
	rm -rf SOCR_regression_test
	rm -rf variogram_test
	rm -rf core_point_test
	rm -rf cluster_test.exe
	rm -rf matrix_test.exe
	rm -rf IGRA_test.exe
	rm -rf krig_test.exe
	rm -rf SOCR_test.exe
	rm -rf generated_test.exe
	rm -rf SOCR_regression_test.exe
	rm -rf variogram_test.exe
	rm -rf core_point_test.exe
