#ifndef DISSIMILARITY_FUNCTIONS_H
#define DISSIMILARITY_FUNCTIONS_H

#include <openipML/SimilarityFunction.h>
#include <openipDS/Image.h>
#include <openipML/DiscreteHistogram.h>
#include <openipML/KMeansClusteringDiscretizer.h>
#include <openipML/EqualWidthDiscretizer.h>
#include <openipML/EqualFrequencyDiscretizer.h>
#include <openipSC/matrixOperations.h>
#include <openipML/PoolAdjacentViolatorsAlgorithm.h>
#include <openipLL/thinning.h>

#include <vector>
#include <exception>
#include <stdexcept>
#include <algorithm>

namespace openip
{
  template<typename T>
  class L1NormDSF: public DissimilarityFunction<T>
  {
  public:
      using DissimilarityFunction<T>::descriptor;
      using DissimilarityFunction<T>::initialized;
      using DissimilarityFunction<T>::apply1;
      using DissimilarityFunction<T>::apply2;
      using DissimilarityFunction<T>::lastTemplate;

      L1NormDSF();

      virtual double evaluate(Vector<T>& a, Vector<T>& b);
      
      virtual double evaluate(Vector<T>& a, Template<T>& b, int n);
      
      virtual void init(Template<T>& b);
      
      virtual DissimilarityFunction<T>* clone();
  };

  template<typename T>
  L1NormDSF<T>::L1NormDSF()
  : DissimilarityFunction<T>()
  {
    descriptor= string("L1NormDSF");
  }

  template<typename T>
  double L1NormDSF<T>::evaluate(Vector<T>& a, Vector<T>& b)
  {
    float result= 0;
    for ( unsigned int i= 0; i < a.size(); ++i )
      result+= fabs(a(i) - b(i));
    return result;
  }
  
  template<typename T>
  double L1NormDSF<T>::evaluate(Vector<T>& a, Template<T>& b, int n)
  {
    if ( !initialized || lastTemplate != &b )
    {
      init(b);
      lastTemplate= &b;
    }
    
    float result= 0;
    for ( unsigned int i= 0; i < b.size(); ++i )
      result+= fabs(a(n + b(i).first) - b(i).second);
    
    return result;
  }
  
  template<typename T>
  void L1NormDSF<T>::init(Template<T>& /*b*/)
  {
    initialized= 1;
  }
  
  template<typename T>
  DissimilarityFunction<T>* L1NormDSF<T>::clone()
  {
    return new L1NormDSF<T>();
  }
  
  template<typename T>
  class MedianOfAbsoluteDifferencesDSF: public DissimilarityFunction<T>
  {
  public:
      using DissimilarityFunction<T>::descriptor;
      using DissimilarityFunction<T>::initialized;
      using DissimilarityFunction<T>::apply1;
      using DissimilarityFunction<T>::apply2;
      using DissimilarityFunction<T>::lastTemplate;

      MedianOfAbsoluteDifferencesDSF();

      virtual double evaluate(Vector<T>& a, Vector<T>& b);
      
      virtual double evaluate(Vector<T>& a, Template<T>& b, int n);
      
      virtual void init(Template<T>& b);
      
      virtual DissimilarityFunction<T>* clone();
      
      Vector<double> diffs;
  };

  template<typename T>
  MedianOfAbsoluteDifferencesDSF<T>::MedianOfAbsoluteDifferencesDSF()
  : DissimilarityFunction<T>()
  {
    descriptor= string("MedianOfAbsoluteDifferencesDSF");
  }

  template<typename T>
  double MedianOfAbsoluteDifferencesDSF<T>::evaluate(Vector<T>& a, Vector<T>& b)
  {
    diffs.resize(a.size());
    
    for ( unsigned int i= 0; i < a.size(); ++i )
      diffs(i)= fabs(a(i) - b(i));
    
    return diffs.getMedian();
  }
  
  template<typename T>
  double MedianOfAbsoluteDifferencesDSF<T>::evaluate(Vector<T>& a, Template<T>& b, int n)
  {
    if ( !initialized || lastTemplate != &b )
    {
      init(b);
      lastTemplate= &b;
    }
    
    diffs.resize(b.size());
    
    for ( unsigned int i= 0; i < b.size(); ++i )
      diffs(i)= fabs(a(n + b(i).first) - b(i).second);
    
    return diffs.getMedian();
  }
  
  template<typename T>
  void MedianOfAbsoluteDifferencesDSF<T>::init(Template<T>& /*b*/)
  {
    initialized= 1;
  }
  
  template<typename T>
  DissimilarityFunction<T>* MedianOfAbsoluteDifferencesDSF<T>::clone()
  {
    return new MedianOfAbsoluteDifferencesDSF<T>();
  }
  
  template<typename T>
  class SquareL2NormDSF: public DissimilarityFunction<T>
  {
  public:
      using DissimilarityFunction<T>::descriptor;
      using DissimilarityFunction<T>::initialized;
      using DissimilarityFunction<T>::apply1;
      using DissimilarityFunction<T>::apply2;
      using DissimilarityFunction<T>::lastTemplate;

      SquareL2NormDSF();

      virtual double evaluate(Vector<T>& a, Vector<T>& b);
      
      virtual double evaluate(Vector<T>& a, Template<T>& b, int n);
      
      virtual void init(Template<T>& b);
      
      virtual DissimilarityFunction<T>* clone();
  };

  template<typename T>
  SquareL2NormDSF<T>::SquareL2NormDSF()
  : DissimilarityFunction<T>()
  {
    descriptor= string("SquareL2NormDSF");
  }

  template<typename T>
  double SquareL2NormDSF<T>::evaluate(Vector<T>& a, Vector<T>& b)
  {
    float result= 0;
    for ( unsigned int i= 0; i < a.size(); ++i )
      result+= (a(i) - b(i))*(a(i) - b(i));
    return result;
  }
  
  template<typename T>
  double SquareL2NormDSF<T>::evaluate(Vector<T>& a, Template<T>& b, int n)
  {
    if ( !initialized || lastTemplate != &b )
    {
      init(b);
      lastTemplate= &b;
    }
    
    float result= 0;
    for ( unsigned int i= 0; i < b.size(); ++i )
      result+= (a(n + b(i).first) - b(i).second)*(a(n + b(i).first) - b(i).second);
    
    return result;
  }
  
  template<typename T>
  void SquareL2NormDSF<T>::init(Template<T>& /*b*/)
  {
    initialized= 1;
  }
  
  template<typename T>
  DissimilarityFunction<T>* SquareL2NormDSF<T>::clone()
  {
    return new SquareL2NormDSF<T>();
  }
  
  template<typename T>
  class MedianOfSquareDifferencesDSF: public DissimilarityFunction<T>
  {
  public:
      using DissimilarityFunction<T>::descriptor;
      using DissimilarityFunction<T>::initialized;
      using DissimilarityFunction<T>::apply1;
      using DissimilarityFunction<T>::apply2;
      using DissimilarityFunction<T>::lastTemplate;

      MedianOfSquareDifferencesDSF();

      virtual double evaluate(Vector<T>& a, Vector<T>& b);
      
      virtual double evaluate(Vector<T>& a, Template<T>& b, int n);
      
      virtual void init(Template<T>& b);
      
      virtual DissimilarityFunction<T>* clone();
      
      Vector<double> diffs;
  };

  template<typename T>
  MedianOfSquareDifferencesDSF<T>::MedianOfSquareDifferencesDSF()
  : DissimilarityFunction<T>()
  {
    descriptor= string("MedianOfSquareDifferencesDSF");
  }

  template<typename T>
  double MedianOfSquareDifferencesDSF<T>::evaluate(Vector<T>& a, Vector<T>& b)
  {
    diffs.resize(a.size());
    
    for ( unsigned int i= 0; i < a.size(); ++i )
      diffs(i)= (a(i) - b(i))*(a(i) - b(i));
    
    return diffs.getMedian();
  }
  
  template<typename T>
  double MedianOfSquareDifferencesDSF<T>::evaluate(Vector<T>& a, Template<T>& b, int n)
  {
    if ( !initialized || lastTemplate != &b )
    {
      init(b);
      lastTemplate= &b;
    }
    
    diffs.resize(b.size());
    
    for ( unsigned int i= 0; i < b.size(); ++i )
      diffs(i)= (a(n + b(i).first) - b(i).second)*(a(n + b(i).first) - b(i).second);
    
    return diffs.getMedian();
  }
  
  template<typename T>
  void MedianOfSquareDifferencesDSF<T>::init(Template<T>& /*b*/)
  {
    initialized= 1;
  }
  
  template<typename T>
  DissimilarityFunction<T>* MedianOfSquareDifferencesDSF<T>::clone()
  {
    return new MedianOfSquareDifferencesDSF<T>();
  }
  
  template<typename T>
  class NormalizedSquareL2NormDSF: public DissimilarityFunction<T>
  {
  public:
      using DissimilarityFunction<T>::descriptor;
      using DissimilarityFunction<T>::initialized;
      using DissimilarityFunction<T>::apply1;
      using DissimilarityFunction<T>::apply2;
      using DissimilarityFunction<T>::lastTemplate;

      NormalizedSquareL2NormDSF();

      virtual double evaluate(Vector<T>& a, Vector<T>& b);
      
      virtual double evaluate(Vector<T>& a, Template<T>& b, int n);
      
      virtual void init(Template<T>& b);
      
      virtual DissimilarityFunction<T>* clone();
      
      float meanb;
      float stdb;
  };

  template<typename T>
  NormalizedSquareL2NormDSF<T>::NormalizedSquareL2NormDSF()
  : DissimilarityFunction<T>()
  {
    descriptor= string("NormalizedSquareL2NormDSF");
    meanb= stdb= 0;
  }

  template<typename T>
  double NormalizedSquareL2NormDSF<T>::evaluate(Vector<T>& a, Vector<T>& b)
  {
    float result= 0;
    float tmp, meana= a.getMean(), stda= a.getStandardDeviation();
    meanb= b.getMean(), stdb= b.getStandardDeviation();
    for ( unsigned int i= 0; i < a.size(); ++i )
    {
      tmp= (a(i) - meana)/stda - (b(i) - meanb)/stdb;
      result= tmp*tmp;
    }
    
    return result;
  }
  
  template<typename T>
  double NormalizedSquareL2NormDSF<T>::evaluate(Vector<T>& a, Template<T>& b, int n)
  {
    if ( !initialized || lastTemplate != &b )
    {
      init(b);
      lastTemplate= &b;
    }
    
    float result= 0;
    float tmp, meana= 0, stda= 0;
    
    for ( unsigned int i= 0; i < b.size(); ++i )
    {
      meana+= a(n + b(i).first);
      stda+= a(n + b(i).first)*a(n + b(i).first);
    }
    meana/= b.size();
    stda/= b.size();
    stda= sqrt(stda - meana*meana);
    
    for ( unsigned int i= 0; i < b.size(); ++i )
    {
      tmp= (a(n + b(i).first) - meana)/stda - (b(i).second - meanb)/stdb;
      result= tmp*tmp;
    }
    
    return result;
  }
  
  template<typename T>
  void NormalizedSquareL2NormDSF<T>::init(Template<T>& b)
  {
    initialized= 1;
    Vector<float> weights;
    b.getWeights(weights);
    meanb= weights.getMean();
    stdb= weights.getStandardDeviation();
  }
  
  template<typename T>
  DissimilarityFunction<T>* NormalizedSquareL2NormDSF<T>::clone()
  {
    return new NormalizedSquareL2NormDSF<T>();
  }
  
  template<typename T>
  class IncrementalSignDistanceDSF: public DissimilarityFunction<T>
  {
  public:
      using DissimilarityFunction<T>::descriptor;
      using DissimilarityFunction<T>::initialized;
      using DissimilarityFunction<T>::apply1;
      using DissimilarityFunction<T>::apply2;
      using DissimilarityFunction<T>::lastTemplate;

      IncrementalSignDistanceDSF();

      virtual double evaluate(Vector<T>& a, Vector<T>& b);
      
      virtual double evaluate(Vector<T>& a, Template<T>& b, int n);
      
      virtual void init(Template<T>& b);
      
      virtual DissimilarityFunction<T>* clone();
  };

  template<typename T>
  IncrementalSignDistanceDSF<T>::IncrementalSignDistanceDSF()
  : DissimilarityFunction<T>()
  {
    descriptor= string("IncrementalSignDistanceDSF");
  }

  template<typename T>
  double IncrementalSignDistanceDSF<T>::evaluate(Vector<T>& a, Vector<T>& b)
  {
    float result= 0;
    
    for ( int i= 0; i < int(a.size())-1; ++i )
      if ( (a(i) > a(i+1) && b(i) <= b(i+1)) || (a(i) <= a(i+1) && b(i) > b(i+1)) )
	result+= 1.0f;
    
    return result;
  }
  
  template<typename T>
  double IncrementalSignDistanceDSF<T>::evaluate(Vector<T>& a, Template<T>& b, int n)
  {
    if ( !initialized || lastTemplate != &b )
    {
      init(b);
      lastTemplate= &b;
    }
    
    float result= 0;
    
    for ( int i= 0; i < int(b.size())-1; ++i )
      if ( (a(n + b(i).first) > a(n + b(i+1).first) && b(i).second <= b(i+1).second) || (a(n + b(i).first) <= a(n + b(i+1).first) && b(i).second > b(i+1).second) )
	result+= 1.0f;
    
    return result;
  }
  
  template<typename T>
  void IncrementalSignDistanceDSF<T>::init(Template<T>& /*b*/)
  {
    initialized= 1;
  }
  
  template<typename T>
  DissimilarityFunction<T>* IncrementalSignDistanceDSF<T>::clone()
  {
    return new IncrementalSignDistanceDSF<T>();
  }
  
  template<typename T>
  class IntensityRatioVarianceDSF: public DissimilarityFunction<T>
  {
  public:
      using DissimilarityFunction<T>::descriptor;
      using DissimilarityFunction<T>::initialized;
      using DissimilarityFunction<T>::apply1;
      using DissimilarityFunction<T>::apply2;
      using DissimilarityFunction<T>::lastTemplate;

      IntensityRatioVarianceDSF();

      virtual double evaluate(Vector<T>& a, Vector<T>& b);
      
      virtual double evaluate(Vector<T>& a, Template<T>& b, int n);
      
      virtual void init(Template<T>& b);
      
      virtual DissimilarityFunction<T>* clone();
      
      Vector<float> tmp;
  };

  template<typename T>
  IntensityRatioVarianceDSF<T>::IntensityRatioVarianceDSF()
  : DissimilarityFunction<T>()
  {
    descriptor= string("IntensityRatioVarianceDSF");
  }

  template<typename T>
  double IntensityRatioVarianceDSF<T>::evaluate(Vector<T>& a, Vector<T>& b)
  {
    float epsilon= 1;
    
    tmp.resize(a.size());
    
    for ( unsigned int i= 0; i < a.size(); ++i )
      tmp(i)= (a(i) + epsilon)/(b(i) + epsilon);
    
    return tmp.getVariance();
  }
  
  template<typename T>
  double IntensityRatioVarianceDSF<T>::evaluate(Vector<T>& a, Template<T>& b, int n)
  {
    if ( !initialized || lastTemplate != &b )
    {
      init(b);
      lastTemplate= &b;
    }
    
    float epsilon= 1;

    tmp.resize(b.size());
    
    for ( unsigned int i= 0; i < b.size(); ++i )
      tmp(i)= (a(n + b(i).first) + epsilon)/(b(i).second + epsilon);
    
    return tmp.getVariance();
  }
  
  template<typename T>
  void IntensityRatioVarianceDSF<T>::init(Template<T>& /*b*/)
  {
    initialized= 1;
  }
  
  template<typename T>
  DissimilarityFunction<T>* IntensityRatioVarianceDSF<T>::clone()
  {
    return new IntensityRatioVarianceDSF<T>();
  }
  
  template<typename T>
  class RankDistanceDSF: public DissimilarityFunction<T>
  {
  public:
      using DissimilarityFunction<T>::descriptor;
      using DissimilarityFunction<T>::initialized;
      using DissimilarityFunction<T>::apply1;
      using DissimilarityFunction<T>::lastTemplate;

      RankDistanceDSF();

      virtual double evaluate(Vector<T>& a, Vector<T>& b);
      
      virtual double evaluate(Vector<T>& a, Template<T>& b, int n);
      
      virtual void init(Template<T>& b);
      
      void rank(Vector<T>& input, Vector<int>& ranked);
      
      double computeDist();
      
      virtual DissimilarityFunction<T>* clone();
      
      //virtual void apply2(Image<T>& input, Template2<T>& temp, Image<T>& output, Image<unsigned char>* roi= NULL, Image<unsigned char>* support= NULL);
      
      Vector<int> rankeda;
      Vector<int> rankedb;
      
      Vector<T> tmp;
  };
  
  template<typename T>
  RankDistanceDSF<T>::RankDistanceDSF()
  : DissimilarityFunction<T>()
  {
    stringstream ss;
    ss << string("RankDistanceDSF");
    descriptor= ss.str();
  }

  template<typename T>
  double RankDistanceDSF<T>::evaluate(Vector<T>& a, Vector<T>& b)
  {
    rank(a, rankeda);
    rank(b, rankedb);
    
    return computeDist();
  }
  
  template<typename T>
  double RankDistanceDSF<T>::evaluate(Vector<T>& a, Template<T>& b, int n)
  {
    if ( !initialized || lastTemplate != &b )
    {
      init(b);
      lastTemplate= &b;
    }
        
    tmp.resize(b.size());
    for ( unsigned int i= 0; i < b.size(); ++i )
      tmp(i)= a(n + b(i).first);
    
    rank(tmp, rankeda);
	
    return computeDist();
  }
  
  template<typename T>
  void RankDistanceDSF<T>::init(Template<T>& b)
  {
    Vector<T> weights;
    b.getWeights(weights);
    
    rank(weights, rankedb);
    
    initialized= 1;
  }
  
  template<typename T>
  void RankDistanceDSF<T>::rank(Vector<T>& b, Vector<int>& ranked)
  {
    Vector<IndexWeightPair<T> > data;
    for ( unsigned int i= 0; i < b.size(); ++i )
      data.push_back(IndexWeightPair<T>(i, b(i)));
    sort(data.begin(), data.end());
    
    ranked.resize(b.size());
    
    for ( unsigned int i= 0; i < data.size(); ++i )
      ranked(data(i).index)= i;
  }
  
  template<typename T>
  double RankDistanceDSF<T>::computeDist()
  {
    float result= 0;
    for ( unsigned int i= 0; i < rankeda.size(); ++i )
      result+= fabs((rankeda(i) - rankedb(i)));
    result/= rankeda.size();
    
    return result;
  }
  
  template<typename T>
  DissimilarityFunction<T>* RankDistanceDSF<T>::clone()
  {
    return new RankDistanceDSF<T>();
  }
  
  //template<typename T>
  //void RankDistanceDSF<T>::apply2(Image<T>& input, Template2<T>& temp, Image<T>& output, Image<unsigned char>* /*roi*/, Image<unsigned char>* /*support*/)
  /*{
    temp.updateStride(input.columns);
    
    Image<T> filtered;
    filtered.resizeImage(input);
    filtered= 0;
    GaussianFilter2<T, T> gf(1, 7);
    gf.updateStride(input.columns);
    gf.apply(input, filtered);
    
    int start= -temp.getMin();
    int end= int(input.n) - temp.getMax();
    
#pragma omp parallel for
    for ( int i= start; i < end; ++i )
      output(i)= evaluate(filtered, temp, i);
  }*/
  
  template<typename T>
  class JointEntropyDSF: public DissimilarityFunction<T>
  {
  public:
      using DissimilarityFunction<T>::descriptor;
      using DissimilarityFunction<T>::initialized;
      using DissimilarityFunction<T>::apply1;
      using DissimilarityFunction<T>::apply2;
      using DissimilarityFunction<T>::lastTemplate;

      JointEntropyDSF();

      virtual double evaluate(Vector<T>& a, Vector<T>& b);
      
      virtual double evaluate(Vector<T>& a, Template<T>& b, int n);
      
      virtual void init(Template<T>& b);
      
      virtual DissimilarityFunction<T>* clone();
      
      DiscreteJointHistogram jpd;
      
      Vector<float> weights;
      Vector<float> values;
  };
  
  template<typename T>
  JointEntropyDSF<T>::JointEntropyDSF()
  : DissimilarityFunction<T>()
  {
    stringstream ss;
    ss << string("JointEntropyDSF ");
    descriptor= ss.str();
  }

  template<typename T>
  double JointEntropyDSF<T>::evaluate(Vector<T>& a, Vector<T>& b)
  {
    jpd.computeHistogram(a, b, 100);
    jpd.normalize();
    
    return jpd.getShannonJointEntropy();
  }
  
  template<typename T>
  double JointEntropyDSF<T>::evaluate(Vector<T>& a, Template<T>& b, int n)
  {
    if ( !initialized || lastTemplate != &b )
    {
      init(b);
      lastTemplate= &b;
    }
    
    values.resize(weights.size());
    for ( unsigned int i= 0; i < b.size(); ++i )
      values(i)= a(n + b(i).first);
    
    return evaluate(values, weights);
  }
  
  template<typename T>
  void JointEntropyDSF<T>::init(Template<T>& b)
  {
    b.getWeights(weights);
    initialized= 1;
  }
  
  template<typename T>
  DissimilarityFunction<T>* JointEntropyDSF<T>::clone()
  {
    return new JointEntropyDSF<T>();
  }
  
  template<typename T>
  class ExclusiveFInformationDSF: public DissimilarityFunction<T>
  {
  public:
      using DissimilarityFunction<T>::descriptor;
      using DissimilarityFunction<T>::initialized;
      using DissimilarityFunction<T>::apply1;
      using DissimilarityFunction<T>::apply2;
      using DissimilarityFunction<T>::lastTemplate;

      ExclusiveFInformationDSF();

      virtual double evaluate(Vector<T>& a, Vector<T>& b);
      
      virtual double evaluate(Vector<T>& a, Template<T>& b, int n);
      
      virtual void init(Template<T>& b);
      
      virtual DissimilarityFunction<T>* clone();
      
      DiscreteJointHistogram jpd;
      
      Vector<float> weights;
      Vector<float> values;
  };
  
  template<typename T>
  ExclusiveFInformationDSF<T>::ExclusiveFInformationDSF()
  : DissimilarityFunction<T>()
  {
    stringstream ss;
    ss << string("ExclusiveFInformationDSF ");
    descriptor= ss.str();
  }

  template<typename T>
  double ExclusiveFInformationDSF<T>::evaluate(Vector<T>& a, Vector<T>& b)
  {
    jpd.computeHistogram(a, b, 100);
    jpd.normalize();
    
    return jpd.getShannonJointEntropy() - jpd.getShannonEntropy();
  }
  
  template<typename T>
  double ExclusiveFInformationDSF<T>::evaluate(Vector<T>& a, Template<T>& b, int n)
  {
    if ( !initialized || lastTemplate != &b )
    {
      init(b);
      lastTemplate= &b;
    }
    
    values.resize(weights.size());
    for ( unsigned int i= 0; i < b.size(); ++i )
      values(i)= a(n + b(i).first);
    
    return evaluate(values, weights);
  }
  
  template<typename T>
  void ExclusiveFInformationDSF<T>::init(Template<T>& b)
  {
    b.getWeights(weights);
    initialized= 1;
  }
  
  template<typename T>
  DissimilarityFunction<T>* ExclusiveFInformationDSF<T>::clone()
  {
    return new ExclusiveFInformationDSF<T>();
  }
}
#endif
