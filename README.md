This software implements filtering-based clusterings algorithm for improving Kriging interpolation accuracy in https://dl.acm.org/doi/abs/10.1145/2983323.2983668.

Main parts: Matrix operations, Kriging interpolations, Clustering functions, and filtering-based clustering.

Input and output files are all text files. SOCR and IGRA datasets referenced in the paper are provided in this folder.

# System requirement

gcc compiler that supports C89.

# Compile

Use make file to compile code and test cases with command "make".


# Data types

Basic types and constants are defined in clustertype.h. For example, DTYPE means double (floating point). DWORD means long long int (integer).
Structure types are defined in cluster.h.

# Executables

All executables have built in inputs without arguments.

- ./IGRA_test	
 * 60 days data are read in. Chi-square values for Kriging with clustering and Kriging without clustering are computed and write to csv files for each time stamps. For example IGRA_compare_0.txt is comparison for time stamp 0.

- ./SOCR_test
 * Test for SOCR data. NMSE, chi-square and number of clusters will be printed. Timing information will also be printed. Each cluster will be written into an individual file that shows the objects in the cluster. For example, SOCR_0.txt contains all objects in cluster with id 0. Comparison of clustered statistics and unclustered statistics will be written into SOCR_compare.txt

- ./matrix_test
 * Unit test for matrix operations. Should pass all test cases without any alerts.

# Contact

Contact me at qiao.kang@eecs.northwestern.edu if you have questions.
