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
  std::map<unsigned int, double> m_predictedProbability;
  bool                           m_labelSet;
public:
  LogisticRegressionSample(const unsigned int featureCount);
  ~LogisticRegressionSample();
  double GetLabelProbability(unsigned int const &);

  void SetSample(std::vector<TSampleType> &);

  std::vector<TSampleType> const * GetSample() const;

  void SetLabelProbability(unsigned int const &, double const &);

  unsigned int GetLabel() const
  {
    return this->m_label;
  };
  void SetLabel(unsigned int const &);

  bool LabelIsSet() const
  {
    return this->m_labelSet;
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
  bool                  m_classOneLabelSet;
  bool                  m_classTwoLabelSet;
public:
  LogisticRegression(const unsigned int featureCount, const unsigned int sampleCount);
  ~LogisticRegression();
  void AddLabeledSample(LogisticRegressionSample<TSampleType> const & );

  void TrainModel();

  void SetClassOneLabel(const unsigned int);

  void SetClassTwoLabel(const unsigned int);

  void ClassifySample(LogisticRegressionSample<TSampleType> &);
};
#include "LogisticRegression.hxx"
#endif
