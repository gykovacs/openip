#include <openipLL/imageIO.h>
#include <openipDS/Image.h>
#include <openipDS/OptionTable.h>
#include <openipLL/imageCorrection.h>
#include <openipML/SimilarityFunctions.h>
#include <openipML/DissimilarityFunctions.h>
#include <openipLL/colorSpaces.h>
#include <openipML/Noises.h>

#include <algorithm>

using namespace openip;

void saveFalseColoredResults(Image<float>& result, char* filename)
{
  Image<float> r, g, b;
  r.resizeImage(result);
  g.resizeImage(result);
  b.resizeImage(result);
  
  unsigned char R, G, B;
  short h, s, v;
  for ( unsigned int i= 0; i < result.n; ++i )
  {
    h= 255-result(i);
    s= 255;
    v= 255;
    hsv2rgb(h, s, v, R, G, B);
    r(i)= R;
    g(i)= G;
    b(i)= B;
  }
  
  writeImage(filename, r, g, b);
}

int main(int argc, char** argv)
{
  srand(time(NULL));
    
  OptionTable ot;

  bool cc= false;
  bool pearsoncc= false;
  bool tanimoto= false;
  bool stochasticsc= false;
  bool deterministicsc= false;
  bool minratio= false;
  bool spearman= false;
  bool kendalltau= false;
  bool greatdev= false;
  bool ordinalm= false;
  bool corrratio= false;
  bool energyjpd= false;
  bool matsim= false;
  bool shannon= false;
  bool renyi= false;
  bool tsallis= false;
  bool finform= false;
  bool l1norm= false;
  bool medabsdiff= false;
  bool squarel2= false;
  bool medsquarediff= false;
  bool normsquarel2= false;
  bool incsigndist= false;
  bool intratiovar= false;
  bool rankdist= false;
  bool jointentropy= false;
  bool exclusivef= false;
  
  bool coloroutput= false;
  
  bool gaussian= false;
  float gaussianstdev= 5;
  
  ot.addOption(string(""), OPTION_SEPARATOR, NULL, 0, string("Similarity/Dissimilarity functions:"));
  ot.addOption(string("--cc"), OPTION_BOOL, (char*)&cc, 0, string("correlation coefficient"));
  ot.addOption(string("--pearsoncc"), OPTION_BOOL, (char*)&pearsoncc, 0, string("Pearson correlation coefficient"));
  ot.addOption(string("--tanimoto"), OPTION_BOOL, (char*)&tanimoto, 0, string("Tanimoto measure"));
  ot.addOption(string("--stochasticsc"), OPTION_BOOL, (char*)&stochasticsc, 0, string("Stochastic sign change"));
  ot.addOption(string("--deterministicsc"), OPTION_BOOL, (char*)&deterministicsc, 0, string("Deterministic sign change"));
  ot.addOption(string("--minratio"), OPTION_BOOL, (char*)&minratio, 0, string("Minimum ratio"));
  ot.addOption(string("--spearmanrho"), OPTION_BOOL, (char*)&spearman, 0, string("Spearman's rho"));
  ot.addOption(string("--kendalltau"), OPTION_BOOL, (char*)&kendalltau, 0, string("Kendall's tau"));
  ot.addOption(string("--greatdev"), OPTION_BOOL, (char*)&greatdev, 0, string("Greatest deviation"));
  ot.addOption(string("--ordinalm"), OPTION_BOOL, (char*)&ordinalm, 0, string("Ordinal measure"));
  ot.addOption(string("--corrratio"), OPTION_BOOL, (char*)&corrratio, 0, string("Correlation ratio"));
  ot.addOption(string("--energyjpd"), OPTION_BOOL, (char*)&energyjpd, 0, string("Energy of joint probability distribution"));
  ot.addOption(string("--matsim"), OPTION_BOOL, (char*)&matsim, 0, string("Material similarity"));
  ot.addOption(string("--shannon"), OPTION_BOOL, (char*)&shannon, 0, string("Shannon mutual information"));
  ot.addOption(string("--renyi"), OPTION_BOOL, (char*)&renyi, 0, string("Renyi mutual information"));
  ot.addOption(string("--tsallis"), OPTION_BOOL, (char*)&tsallis, 0, string("Tsallis mutual information"));
  ot.addOption(string("--finform"), OPTION_BOOL, (char*)&finform, 0, string("F-information measures"));
  ot.addOption(string("--l1norm"), OPTION_BOOL, (char*)&l1norm, 0, string("L1 norm"));
  ot.addOption(string("--medabsdiff"), OPTION_BOOL, (char*)&medabsdiff, 0, string("Median of absolute differences"));
  ot.addOption(string("--squarel2"), OPTION_BOOL, (char*)&squarel2, 0, string("Square L2 norm"));
  ot.addOption(string("--medsquarediff"), OPTION_BOOL, (char*)&medsquarediff, 0, string("Median of square differences"));
  ot.addOption(string("--normsquarel2"), OPTION_BOOL, (char*)&normsquarel2, 0, string("Normalized square L2 norm"));
  ot.addOption(string("--incsigndist"), OPTION_BOOL, (char*)&incsigndist, 0, string("Incremental sign distance"));
  ot.addOption(string("--intratiovar"), OPTION_BOOL, (char*)&intratiovar, 0, string("Intensity-ratio variance"));
  ot.addOption(string("--rankdist"), OPTION_BOOL, (char*)&rankdist, 0, string("Rank distance"));
  ot.addOption(string("--jointentropy"), OPTION_BOOL, (char*)&jointentropy, 0, string("Joint entropy"));
  ot.addOption(string("--exclusivef"), OPTION_BOOL, (char*)&exclusivef, 0, string("Exclusive F-information"));
  
  ot.addOption(string(""), OPTION_SEPARATOR, NULL, 0, string("Other options:"));
  ot.addOption(string("--coloroutput"), OPTION_BOOL, (char*)&coloroutput, 0, string("false colored output"));
  
  ot.addOption(string(""), OPTION_SEPARATOR, NULL, 0, string("Noises:"));
  ot.addOption(string("--gaussian"), OPTION_BOOL, (char*)&gaussian, 0, string("add Gaussian noise"));
  ot.addOption(string("--gaussian.stdev"), OPTION_FLOAT, (char*)&gaussianstdev, 1, string("standard deviation of Gaussian noise"));
  
  if ( ot.processArgs(&argc, argv) )
      return 0;

  DissimilarityFunction<float>* dsf= NULL;
  SimilarityFunction<float>* sf= NULL;
  
  tprintf("instantiating similarity function\n");
  if ( cc )
    sf= new CorrelationCoefficientSF<float>();
  else if ( pearsoncc )
    sf= new PearsonCorrelationCoefficientSF<float>();
  else if ( tanimoto )
    sf= new TanimotoSF<float>();
  else if ( stochasticsc )
    sf= new StochasticSignChangeSF<float>();
  else if ( deterministicsc )
    sf= new DeterministicSignChangeSF<float>();
  else if ( minratio )
    sf= new MinimumRatioSF<float>();
  else if ( spearman )
    sf= new SpearmanRhoSF<float>();
  else if ( kendalltau )
    sf= new KendallTauSF<float>();
  else if ( greatdev )
    sf= new GreatestDeviationSF<float>();
  else if ( ordinalm )
    sf= new OrdinalMeasureSF<float>();
  else if ( corrratio )
    sf= new CorrelationRatioSF<float>();
  else if ( energyjpd )
    sf= new EnergyOfJointProbabilityDistributionSF<float>();
  else if ( matsim )
    sf= new MaterialSF<float>();
  else if ( shannon )
    sf= new ShannonMutualInformationSF<float>();
  else if ( renyi )
    sf= new RenyiMutualInformationSF<float>();
  else if ( tsallis )
    sf= new TsallisMutualInformationSF<float>();
  else if ( finform )
    sf= new FInformationMeasureISF<float>();
  
  tprintf("instantiating dissimilarity function\n");
  if ( l1norm )
    dsf= new L1NormDSF<float>();
  else if ( medabsdiff )
    dsf= new MedianOfAbsoluteDifferencesDSF<float>();
  else if ( squarel2 )
    dsf= new SquareL2NormDSF<float>();
  else if ( medsquarediff )
    dsf= new MedianOfSquareDifferencesDSF<float>();
  else if ( normsquarel2 )
    dsf= new NormalizedSquareL2NormDSF<float>();
  else if ( incsigndist )
    dsf= new IncrementalSignDistanceDSF<float>();
  else if ( intratiovar )
    dsf= new IntensityRatioVarianceDSF<float>();
  else if ( rankdist )
    dsf= new RankDistanceDSF<float>();
  else if ( jointentropy )
    dsf= new JointEntropyDSF<float>();
  else if ( exclusivef )
    dsf= new ExclusiveFInformationDSF<float>();
  
  if ( sf != NULL )
    dsf= new DissimilarityFunction<float>(sf);
  
  tprintf("reading images\n");
  Image<float> input;
  Image<float> temp;
  Image<float> output;
  
  readImage(argv[1], input);
  readImage(argv[2], temp);
  
  Template2<float> temp0(temp);
  
  Border2 b= temp0.getProposedBorder();
  b.borderMode= BORDER_MODE_MIRRORED;
  input.setBorder(b);
  output.resizeImage(input);
  output= 0;
  
  tprintf("starting template matching\n");
  dsf->apply2(input, temp0, output);
  
  tprintf("template matching finished\n");
  
  output.removeBorder();
  
  tprintf("normalizing output\n");
  output.normalize(0, 255);
  
  tprintf("writing output\n");
  if ( coloroutput )
    saveFalseColoredResults(output, argv[3]);
  else
    writeImage(argv[3], output);
  
  return 0;
}




