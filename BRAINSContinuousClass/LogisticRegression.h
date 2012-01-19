#ifndef LOGISTICREGRESSION_h
#define LOGISTICREGRESSION_h

#include "linear.h"
#include <vector>
#include <map>

template <typename TSampleType>
class LogisticRegressionSample
{
private:
  std::vector<TSampleType> *     m_sample;
  unsigned int                   m_label;
  std::map<unsigned int, double> m_predictedProbability[2];
  bool                           m_ReversedClassOrder;
public:
  LogisticRegressionSample(const unsigned int featureCount);
  ~LogisticRegressionSample();
  double const GetLabelProbability(unsigned int const &) const;

  void SetSample(std::vector<TSampleType> *);

  std::vector<TSampleType> const * GetSample() const;

  void SetLabelProbability(unsigned int const &, double const &);

  unsigned int const GetLabel() const
  {
    return this->m_label;
  };
};

template <typename TSampleType>
class LogisticRegression
{
private:
  unsigned int          m_sampleCount;
  unsigned int          m_totalSamples;
  struct parameter      m_parameters;
  struct problem        m_problem;
  struct model *        m_model;
  struct feature_node * m_featureNodes;
  unsigned int          m_featureCount;
  unsigned int          m_classOneLabel;
  unsigned int          m_classTwoLabel;
public:
  LogisticRegression(const unsigned int featureCount, const unsigned int sampleCount);
  ~LogisticRegression();
  void AddLabeledSample(LogisticRegressionSample<TSampleType> const & );

  void TrainModel();

  void SetClassOneLabel(const unsigned int);

  void SetClassTwoLabel(const unsigned int);

  void ClassifySample(LogisticRegressionSample<TSampleType> &);
};
#endif
