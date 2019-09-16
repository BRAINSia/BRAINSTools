/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef _LIBLINEAR_H
#define _LIBLINEAR_H

#ifdef __cplusplus
extern "C"
{
#endif

  struct feature_node
  {
    int    index;
    double value;
  };

  struct problem
  {
    int                    l, n;
    int *                  y;
    struct feature_node ** x;
    double                 bias; /* < 0 if no bias term */
  };

  enum
  {
    L2R_LR,
    L2R_L2LOSS_SVC_DUAL,
    L2R_L2LOSS_SVC,
    L2R_L1LOSS_SVC_DUAL,
    MCSVM_CS,
    L1R_L2LOSS_SVC,
    L1R_LR,
    L2R_LR_DUAL
  }; /*
       solver_type
       */

  struct parameter
  {
    int solver_type;

    /* these are for training only */
    double   eps; /* stopping criteria */
    double   C;
    int      nr_weight;
    int *    weight_label;
    double * weight;
  };

  struct model
  {
    struct parameter param;
    int              nr_class; /* number of classes */
    int              nr_feature;
    double *         w;
    int *            label; /* label of each class */
    double           bias;
  };

  struct model *
  train(const struct problem * prob, const struct parameter * param);

  void
  cross_validation(const struct problem * prob, const struct parameter * param, int nr_fold, int * target);

  int
  predict_values(const struct model * model_, const struct feature_node * x, double * dec_values);

  int
  predict(const struct model * model_, const struct feature_node * x);

  int
  predict_probability(const struct model * model_, const struct feature_node * x, double * prob_estimates);

  int
  save_model(const char * model_file_name, const struct model * model_);

  struct model *
  load_model(const char * model_file_name);

  int
  get_nr_feature(const struct model * model_);

  int
  get_nr_class(const struct model * model_);

  void
  get_labels(const struct model * model_, int * label);

  void
  free_model_content(struct model * model_ptr);

  void
  free_and_destroy_model(struct model ** model_ptr_ptr);

  void
  destroy_param(struct parameter * param);

  const char *
  check_parameter(const struct problem * prob, const struct parameter * param);

  int
  check_probability_model(const struct model * model);

  void
  set_print_string_function(void (*print_func)(const char *));

#ifdef __cplusplus
}
#endif

#endif /* _LIBLINEAR_H */
