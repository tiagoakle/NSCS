#include <stdio.h>
#include <check.h>
#include <stdlib.h>
//Include the test headers
#include "test_linear_solver.h"
#include "test_smatvec.h"
#include "test_barriers.h"


int main(void)
{

 //The multithreading of open blas conflicts with 
 //the multithreading of check :(. This solves
 //the issue but the tests might become less representative 
 //of the real behavior of the program.
 openblas_set_num_threads(1);


 //smatvec tests
 int number_failed;
 Suite *s = matvec_suite();
 SRunner *sr = srunner_create (s);
 srunner_run_all (sr, CK_VERBOSE);
 number_failed = srunner_ntests_failed (sr);
 srunner_free (sr);

 //Barriers
 s = barriers_suite();
 sr = srunner_create(s);
 srunner_run_all(sr, CK_VERBOSE);
 number_failed += srunner_ntests_failed(sr);
 srunner_free(sr);

 //linear solver 
 s  = linear_solver_suite();
 sr = srunner_create (s);
 srunner_run_all (sr, CK_VERBOSE);
 number_failed = srunner_ntests_failed (sr);
 srunner_free (sr);

 return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
