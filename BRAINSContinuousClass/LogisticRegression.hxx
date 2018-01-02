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
#include "LogisticRegression.h"
#include "linear.h"
#include <map>
#include <vector>

template <typename TSampleType>
LogisticRegressionSample<TSampleType>::LogisticRegressionSample(const unsigned int featureCount)
{
  this->m_sample = new std::vector<TSampleType>(featureCount + 1);
  this->m_label = 0;
  this->m_labelSet = false;
}

template <typename TSampleType>
LogisticRegressionSample<TSampleType>::LogisticRegressionSample(const LogisticRegressionSample & LRS)
{
  this->m_sample = LRS.m_sample;
  this->m_label = 0;
  this->m_labelSet = false;
}

template <typename TSampleType>
LogisticRegressionSample<TSampleType>::~LogisticRegressionSample()
{
  delete this->m_sample;
}

template <typename TSampleType>
double LogisticRegressionSample<TSampleType>::GetLabelProbability(unsigned int const & label)
{
  return this->m_predictedProbability[label];
}

template <typename TSampleType>
void LogisticRegressionSample<TSampleType>::SetLabelProbability(unsigned int const & label, double const & probability)
{
  this->m_predictedProbability[label] = probability;
}

template <typename TSampleType>
std::vector<TSampleType> const * LogisticRegressionSample<TSampleType>::GetSample() const
{
  return this->m_sample;
}

template <typename TSampleType>
void LogisticRegressionSample<TSampleType>::SetSample(std::vector<TSampleType> & sample)
{
  *this->m_sample = sample;
}

template <typename TSampleType>
void LogisticRegressionSample<TSampleType>::SetLabel(unsigned int const & label)
{
  this->m_label = label;
  this->m_labelSet = true;
}

template <typename TSampleType>
LogisticRegression<TSampleType>::LogisticRegression(const unsigned int featureCount, const unsigned int totalSamples) :
  m_totalSamples(totalSamples),
  m_featureCount(featureCount),
  m_classOneLabel(0),
  m_classTwoLabel(0),
  m_classOneLabelSet(false),
  m_classTwoLabelSet(false)
{
  this->m_sampleCount = 0;
  this->m_parameters.solver_type = L1R_LR;
  this->m_parameters.C = 1;
  this->m_parameters.eps = 0.01;
  this->m_parameters.nr_weight = 0;
  this->m_parameters.weight_label = nullptr;
  this->m_parameters.weight = nullptr;
  this->m_problem.bias = 1;
  this->m_problem.n = this->m_problem.bias + this->m_featureCount;
  this->m_problem.l = this->m_totalSamples;
  this->m_problem.y = new int[totalSamples];
  this->m_problem.x = new struct feature_node *[this->m_problem.l];
  this->m_featureNodes = new struct feature_node[this->m_problem.n * this->m_problem.l];
  this->m_model = nullptr;
}

template <typename TSampleType>
LogisticRegression<TSampleType>::LogisticRegression(const LogisticRegression & LR) :
  m_totalSamples(LR.m_totalSamples),
  m_featureCount(LR.m_featureCount),
  m_classOneLabel(0),
  m_classTwoLabel(0),
  m_classOneLabelSet(false),
  m_classTwoLabelSet(false)
{
  this->m_sampleCount = 0;
  this->m_parameters.solver_type = L1R_LR;
  this->m_parameters.C = 1;
  this->m_parameters.eps = 0.01;
  this->m_parameters.nr_weight = 0;
  this->m_parameters.weight_label = nullptr;
  this->m_parameters.weight = nullptr;
  this->m_problem.bias = 1;
  this->m_problem.n = LR.m_problem.bias + LR.m_featureCount;
  this->m_problem.l = LR.m_totalSamples;
  this->m_problem.y = new int[LR.m_totalSamples];
  this->m_problem.x = new struct feature_node *[LR.m_problem.l];
  this->m_featureNodes = new struct feature_node[LR.m_problem.n * LR.m_problem.l];
  this->m_model = nullptr;
}

template <typename TSampleType>
LogisticRegression<TSampleType>::~LogisticRegression()
{
  delete this->m_problem.y;
  delete this->m_problem.x;
  delete this->m_featureNodes;
  destroy_param(&this->m_parameters);
}

template <typename TSampleType>
void LogisticRegression<TSampleType>::AddLabeledSample(LogisticRegressionSample<TSampleType> const & labeledSample)
{
  assert(labeledSample.LabelIsSet() );

  std::vector<TSampleType> const * const samples = labeledSample.GetSample();

  struct feature_node * currentNode;

  this->m_problem.x[this->m_sampleCount] = &this->m_featureNodes[this->m_problem.n * this->m_sampleCount];
  this->m_problem.y[this->m_sampleCount]  = labeledSample.GetLabel();
  for( unsigned int i = 0; i < samples->size(); ++i )
    {
    currentNode = &this->m_featureNodes[this->m_problem.n * this->m_sampleCount + i];
    currentNode->index = i + 1;
    currentNode->value = samples->at(i);
    }
  currentNode = &this->m_featureNodes[this->m_problem.n * this->m_sampleCount + samples->size()];
  currentNode->index = -1;

  this->m_sampleCount++;
}

template <typename TSampleType>
void LogisticRegression<TSampleType>::TrainModel()
{
  assert(this->m_totalSamples >= this->m_sampleCount);
  if( this->m_totalSamples > this->m_sampleCount )
    {
    m_problem.l = m_sampleCount;
    }
  this->m_model = train(&this->m_problem, &this->m_parameters);
}

template <typename TSampleType>
void LogisticRegression<TSampleType>::SetClassOneLabel(const unsigned int classLabel)
{
  this->m_classOneLabel = classLabel;
  this->m_classOneLabelSet = true;
}

template <typename TSampleType>
void LogisticRegression<TSampleType>::SetClassTwoLabel(const unsigned int classLabel)
{
  this->m_classTwoLabel = classLabel;
  this->m_classTwoLabelSet = true;
}

template <typename TSampleType>
void LogisticRegression<TSampleType>::ClassifySample(LogisticRegressionSample<TSampleType> & labeledSample)
{
  // TODO Check if the class labels have been set and throw an exception if they haven't.
  assert(this->m_classTwoLabelSet && this->m_classOneLabelSet);

  std::vector<TSampleType> const * const samples = labeledSample.GetSample();
  struct feature_node *                  sampleToPredict = new struct feature_node[this->m_problem.n];
  for( unsigned int i = 0; i < samples->size(); ++i )
    {
    sampleToPredict[i].index = i + 1;
    sampleToPredict[i].value = samples->at(i);
    }
  sampleToPredict[samples->size()].index = -1;

  double predictedProbabilities[2];
  predict_probability(this->m_model, sampleToPredict, predictedProbabilities);

  if( this->m_classOneLabel > this->m_classTwoLabel )
    {
    labeledSample.SetLabelProbability(this->m_classOneLabel, predictedProbabilities[1]);
    labeledSample.SetLabelProbability(this->m_classTwoLabel, predictedProbabilities[0]);
    }
  else
    {
    labeledSample.SetLabelProbability(this->m_classOneLabel, predictedProbabilities[0]);
    labeledSample.SetLabelProbability(this->m_classTwoLabel, predictedProbabilities[1]);
    }

  delete [] sampleToPredict;
}
