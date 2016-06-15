//Author: G.Vianello (giacomov@slac.stanford.edu)

//This implements a real ModelInterface, which is used to bridge between a 3ML
//LikelihoodModel class (living in the python world) and plug-ins which are
//living in the C++ world. The plugins will only use the pyToCppModelInterface class
//to talk (without knowing it) to the python LikelihoodModel class.

#ifndef PYTOCPPMODEL_INTERFACE_CACHE_H
#define PYTOCPPMODEL_INTERFACE_CACHE_H

#include <boost/python.hpp>

#include <Python.h>

#include "ModelInterface.h"
#include <vector>
#include <map>
#include <string>
#include <boost/ptr_container/ptr_map.hpp>

using namespace boost::python;

namespace threeML {

typedef std::pair<float, float> SkyCoord;

typedef struct BoundingBox
{
  float lon_min;
  float lon_max;
  float lat_min;
  float lat_max;
} BoundingBox;

typedef std::vector<double> Spectrum;

// Using ptr_map from boost to avoid the copying that the standard
// map would make every time we insert an element. These are also
// auto pointers, so the moment the map is emptied, the memory
// will be freed

typedef boost::ptr_map<SkyCoord, Spectrum> ExtSrcCube;


class pyToCppModelInterfaceCache: public ModelInterface {

 public:

  pyToCppModelInterfaceCache() : m_nPtSources(0), m_nExtSources(0) {};

  void setExtSourceBoundaries(const int id, const float lon_min, const float lon_max,
                              const float lat_min, const float lat_max);

  void setExtSourceCube(const int id, const numeric::array& cube,
                        const numeric::array& lon, const numeric::array& lat);

  void setPtsSourceSpectrum(const int id, const numeric::array& spectrum);

  void setPtsSourcePosition(const int id, const float lon, const float lat);

  //Point source interface

  int getNumberOfPointSources() const;

  void getPointSourcePosition(int srcid, double *j2000_ra, double *j2000_dec) const;

  //Fluxes are differential fluxes in MeV^-1 cm^-1 s^-1
  std::vector<double> getPointSourceFluxes(int srcid, std::vector<double> energies) const;

  std::string getPointSourceName(int srcid) const;

  //Extended source interface

  int getNumberOfExtendedSources() const;

  std::vector<double> getExtendedSourceFluxes(int srcid, double j2000_ra, double j2000_dec,
                                              std::vector<double> energies) const;

  std::vector<double> getExtendedSourceFluxes_test(int srcid, double j2000_ra, double j2000_dec) const;

  std::string getExtendedSourceName(int srcid) const;

  bool isInsideAnyExtendedSource(double j2000_ra, double j2000_dec) const;

  void getExtendedSourceBoundaries(int srcid, double *j2000_ra_min,
                                   double *j2000_ra_max,
                                   double *j2000_dec_min,
                                   double *j2000_dec_max) const;

  void reset();

 private:

  int m_nPtSources, m_nExtSources;

  std::map<int, ExtSrcCube> m_extSources;

  std::map<int, BoundingBox> m_boundingBoxes;

  boost::ptr_map<int, Spectrum> m_ptsSources;

  std::map<int, SkyCoord> m_ptsSourcesPos;

};

}

#endif
