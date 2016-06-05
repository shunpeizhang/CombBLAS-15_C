# CMake generated Testfile for 
# Source directory: /home/allenzou/git/CombBLAS
# Build directory: /home/allenzou/git/CombBLAS
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(Multiplication_Test "mpirun" "-np" "4" "./ReleaseTests/MultTest" "TESTDATA/rmat_scale16_A.txt" "TESTDATA/rmat_scale16_B.txt" "TESTDATA/rmat_scale16_productAB.txt" "TESTDATA/x_65536_halfdense.txt" "TESTDATA/y_65536_halfdense.txt")
ADD_TEST(Reduction_Test "mpirun" "-np" "4" "./ReleaseTests/ReduceTest" "TESTDATA/sprand10000" "TESTDATA/sprand10000_sumcols" "TESTDATA/sprand10000_sumrows")
ADD_TEST(Iterator_Test "mpirun" "-np" "4" "./ReleaseTests/IteratorTest" "TESTDATA" "sprand10000")
ADD_TEST(Transpose_Test "mpirun" "-np" "4" "./ReleaseTests/TransposeTest" "TESTDATA" "betwinput_scale16" "betwinput_transposed_scale16")
ADD_TEST(Indexing_Test "mpirun" "-np" "4" "./ReleaseTests/IndexingTest" "TESTDATA" "B_100x100.txt" "B_10x30_Indexed.txt" "rand10outta100.txt" "rand30outta100.txt")
ADD_TEST(SpAsgn_Test "mpirun" "-np" "4" "./ReleaseTests/SpAsgnTest" "TESTDATA" "A_100x100.txt" "A_with20x30hole.txt" "dense_20x30matrix.txt" "A_wdenseblocks.txt" "20outta100.txt" "30outta100.txt")
ADD_TEST(GalerkinNew_Test "mpirun" "-np" "4" "./ReleaseTests/GalerkinNew" "TESTDATA/grid3d_k5.txt" "TESTDATA/offdiag_grid3d_k5.txt" "TESTDATA/diag_grid3d_k5.txt" "TESTDATA/restrict_T_grid3d_k5.txt")
ADD_TEST(FindSparse_Test "mpirun" "-np" "4" "./ReleaseTests/FindSparse" "TESTDATA" "findmatrix.txt")
ADD_TEST(BetwCent_Test "mpirun" "-np" "4" "./Applications/betwcent" "TESTDATA/SCALE16BTW-TRANSBOOL/" "10" "96")
ADD_TEST(TopDownBFS_Test "mpirun" "-np" "4" "./Applications/tdbfs" "Force" "17" "FastGen")
ADD_TEST(DirOptBFS_Test "mpirun" "-np" "4" "./Applications/dobfs" "17")
ADD_TEST(FBFS_Test "mpirun" "-np" "4" "./Applications/fbfs" "Gen" "16")
ADD_TEST(FMIS_Test "mpirun" "-np" "4" "./Applications/fmis" "17")
ADD_TEST(Mcl_Test "mpirun" "-np" "4" "./Applications/ppcl" "TESTDATA/SCALE168BTW-TRANSBOOL/" "2" "0.00001")
ADD_TEST(PPCL_Test "mpirun" "-np" "4" "./Applications/ppcl" "TESTDATA/SCALE168BTW-TRANSBOOL/")
SUBDIRS(ReleaseTests)
SUBDIRS(Applications)
SUBDIRS(graph500-1.2/generator)
