#ifndef H_LINE_SEARCH
#define H_LINE_SEARCH
//int linesearch_atd(state_t state ,parameters_t  pars, problem_t prob);
#ifdef __cplusplus
extern "C" {
#endif

int linesearch_atd(state_t* state , parameters_t  pars, problem_t prob);
int linesearch_centering(state_t* state , parameters_t  pars, problem_t prob);

#ifdef __cplusplus
}
#endif


#endif
