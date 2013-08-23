#include <stdio.h>
#include <check.h>
#include <stdlib.h>
//Include the test headers
#include "test_linear_solver.h"
#include "test_smatvec.h"



int main(void)
{

  //smatvec tests
  int number_failed;
  Suite *s = matvec_suite();
  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
 
  //linear solver 
  s  = linear_solver_suite();
  sr = srunner_create (s);
  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
 
  //XXX:
  //For some reason I have not deduced
  //calling the test using check hangs.....
  //so this test was removed from the suite
  test_load_csv_and_solve_system();

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
