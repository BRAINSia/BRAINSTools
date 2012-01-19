#include "LogisticRegression.h"
#include "linear.h"
#include <map>
#include <vector>

template class LogisticRegression<float>;
template class LogisticRegressionSample<float>;

template <typename TSampleType>
LogisticRegression<TSampleType>::LogisticRegression(const unsigned int featureCount, const unsigned int totalSamples)
{
  this->m_sampleCount = 0;
  this->m_totalSamples = totalSamples;
  this->m_parameters.solver_type = L1R_LR;
  this->m_parameters.C = 1;
  this->m_parameters.eps = 0.01;
  this->m_parameters.nr_weight = 0;
  this->m_parameters.weight_label = NULL;
  this->m_parameters.weight = NULL;
  this->m_featureCount = featureCount;
  this->m_problem.bias = 1;
  this->m_problem.n = this->m_problem.bias + this->m_featureCount;
  this->m_problem.l = this->m_totalSamples;
  this->m_problem.y = new int[totalSamples];
  this->m_problem.x = new struct feature_node *[this->m_problem.l];
  this->m_featureNodes = new struct feature_node[this->m_problem.n * this->m_problem.l];
}

template <typename TSampleType>
LogisticRegression<TSampleType>::~LogisticRegression()
{
  delete this->m_problem.y;
  delete this->m_problem.x;
  delete this->m_featureNodes;
}

template <typename TSampleType>
void LogisticRegression<TSampleType>::AddLabeledSample(LogisticRegressionSample<TSampleType> const & labeledSample)
{
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
  this->m_model = train(&this->m_problem, &this->m_parameters);
}

template <typename TSampleType>
void LogisticRegression<TSampleType>::ClassifySample(LogisticRegressionSample<TSampleType> & labeledSample)
{
  std::vector<TSampleType> const * const samples = labeledSample.GetSample();
  struct feature_node                    sampleToPredict[static_cast<int>(m_featureCount + this->m_problem.bias)];

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
}

template <typename TSampleType>
LogisticRegressionSample<TSampleType>::LogisticRegressionSample(const unsigned int featureCount)
{
  this->m_sample = new std::vector<TSampleType>[featureCount + 1];
}

template <typename TSampleType>
LogisticRegressionSample<TSampleType>::~LogisticRegressionSample()
{
  delete this->m_sample;
}

template <typename TSampleType>
double const LogisticRegressionSample<TSampleType>::GetLabelProbability(unsigned int const & label) const
{
  return this->m_predictedProbability->at(label);
}

template <typename TSampleType>
void LogisticRegressionSample<TSampleType>::SetLabelProbability(unsigned int const & label, double const & probability)
{
  this->m_predictedProbability->at(label) = probability;
}

template <typename TSampleType>
std::vector<TSampleType> const * LogisticRegressionSample<TSampleType>::GetSample() const
{
  return this->m_sample;
}

template <typename TSampleType>
void LogisticRegressionSample<TSampleType>::SetSample(std::vector<TSampleType> * sample)
{
  this->m_sample = sample;
}
