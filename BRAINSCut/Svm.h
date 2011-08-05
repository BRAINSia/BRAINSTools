#ifndef __Svm_h__
#define __Svm_h__
#include "Utilities.h"
struct svm_node
  {
  int index;
  double value;
  };
struct svm_problem
  {
  int l_svm;
  double *y_svm;
  struct svm_node * *x_svm;
  };
enum
  {
  C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR
  };   /* svm_type */
enum
  {
  LINEAR, POLY, RBF, SIGMOID
  };   /* kernel_type */
struct svm_parameter
  {
  int svm_type;
  int kernel_type;
  double degree;       /* for poly */
  double gamma;        /* for poly/rbf/sigmoid */
  double coef0;        /* for poly/sigmoid */
  double cache_size;   /* in MB */
  double eps;          /* stopping criteria */
  double C;            /* for C_SVC, EPSILON_SVR and NU_SVR */
  int nr_weight;       /* for C_SVC */
  int *weight_label;   /* for C_SVC */
  double *weight;      /* for C_SVC */
  double nu;           /* for NU_SVC, ONE_CLASS, and NU_SVR */
  double p;            /* for EPSILON_SVR */
  int shrinking;       /* use the shrinking heuristics */
  int probability;     /* do probability estimates */
  };
struct svm_model * svm_train(const struct svm_problem *prob, const struct svm_parameter *param);

int svm_save_model(const char *model_file_name, const struct svm_model *model);

struct svm_model * svm_load_model(const char *model_file_name);

int svm_get_svm_type(const struct svm_model *model);

double svm_get_svr_probability(const struct svm_model *model);

void svm_predict_values(const struct svm_model *model, const struct svm_node *x, double *dec_values);

double svm_predict_probability(const struct svm_model *model, const struct svm_node *x, double *prob_estimates);

void svm_destroy_model(struct svm_model *model);

void svm_destroy_param(struct svm_parameter *param);

const char * svm_check_parameter(const struct svm_problem *prob, const struct svm_parameter *param);

int svm_check_probability_model(const struct svm_model *model);

#endif
