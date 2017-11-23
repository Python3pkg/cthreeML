#include <boost/python.hpp>
#include "ModelInterface.h"
#include "pyToCppModelInterfaceCache.h"

#include <boost/python/stl_iterator.hpp>
#include <iostream>
#include <stdexcept>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <numpy/noprefix.h>


namespace threeML {

template<typename T>
inline
std::vector<T> to_std_vector(const boost::python::object &iterable) {
  return std::vector<T>(boost::python::stl_input_iterator<T>(iterable),
                        boost::python::stl_input_iterator<T>());
}

void pyToCppModelInterfaceCache::setPtsSourceSpectrum(const int id, const numeric::array& spectrum)
{

  // These are both n_points in size

  PyArrayObject* spectrum__ = (PyArrayObject*) spectrum.ptr();

  unsigned int n_energies = (unsigned int) *(spectrum__->dimensions);

  double* spectrum_ = (double *) spectrum__->data;

  m_ptsSources[id].assign(spectrum_, spectrum_ + n_energies);

}

void pyToCppModelInterfaceCache::setPtsSourcePosition(const int id, const float lon, const float lat)
{

  m_ptsSourcesPos[id] = SkyCoord(lon, lat);

}

void pyToCppModelInterfaceCache::setExtSourceBoundaries(const int id, const float lon_min, const float lon_max,
                                                        const float lat_min, const float lat_max)
{

  BoundingBox this_boundingBox;

  this_boundingBox.lon_min = lon_min;
  this_boundingBox.lon_max = lon_max;
  this_boundingBox.lat_min = lat_min;
  this_boundingBox.lat_max = lat_max;

  m_boundingBoxes[id] = this_boundingBox;

}

void pyToCppModelInterfaceCache::setExtSourceCube(const int id, const numeric::array &cube,
                                                  const numeric::array &lon, const numeric::array &lat)
{

  // Add an extended source to the cache (a map)

  PyArrayObject* cube__ = (PyArrayObject*) cube.ptr();

  unsigned int n_points = (unsigned int) *(cube__->dimensions);
  unsigned int n_energies = (unsigned int) *(cube__->dimensions+1);

  // This is n_points * n_energies in size

  double* cube_ = (double *) cube__->data;

  // These are both n_points in size

  double* lon_ = (double *) ((PyArrayObject*) lon.ptr())->data;
  double* lat_ = (double *) ((PyArrayObject*) lat.ptr())->data;

  ExtSrcCube this_srcCube;

  for(unsigned int i=0; i < n_points; ++i)
  {

    float this_lon = (float) lon_[i];
    float this_lat = (float) lat_[i];

    SkyCoord this_coord(this_lon, this_lat);

    // A little of pointer algebra

    int this_data_start = i * n_energies;
    int this_data_stop = this_data_start + n_energies;

    // Note how this construct avoid the double copy which would
    // occur if we were to first create a vector and then
    // copy it into a standard map. Here instead the operator
    // [] creates an empty vector, which is then filled with
    // the assign method

    this_srcCube[this_coord].assign(cube_ + this_data_start, cube_ + this_data_stop);

  }

  // Finally add to the map

  m_extSources[id] = this_srcCube;

}

bool pyToCppModelInterfaceCache::isInsideAnyExtendedSource(double j2000_ra, double j2000_dec) const {
  return true;
}

int pyToCppModelInterfaceCache::getNumberOfPointSources() const {

  return m_ptsSourcesPos.size();

}

void pyToCppModelInterfaceCache::getPointSourcePosition(int srcid, double *j2000_ra, double *j2000_dec) const {

  try {

    SkyCoord this_position = m_ptsSourcesPos.at(srcid);
    *j2000_ra = this_position.first;
    *j2000_dec = this_position.second;

  } catch (...) {

    std::stringstream name;

    name << "Point source " << srcid << " not found in position cache";

    throw std::runtime_error(name.str());

  }

}


void pyToCppModelInterfaceCache::reset() {

  //Empty the cache

  m_extSources.clear();
  m_ptsSources.clear();
  /*
  m_boundingBoxes.clear();
  m_nExtSources = 0;
  m_nPtSources = 0;*/

}

std::vector<double>
pyToCppModelInterfaceCache::getPointSourceFluxes(int srcid, std::vector<double> energies) const {


  if (m_ptsSources.count(srcid) == 0)
  {
    // This happens during the construction of LikeHAWC because there is a energy reweighting
    // At that moment the cache is still empty, so let's just return the energy array which
    // has the right size.

    /*
    std::stringstream name;

    name << "Point source " << srcid << " not found in spectrum cache";


    throw std::runtime_error(name.str());*/

    std::cerr << "Point source " << srcid << " not found in spectrum cache" << std::endl;

    return energies;

  } else {

    return m_ptsSources.at(srcid);

  }



}

int pyToCppModelInterfaceCache::getNumberOfExtendedSources() const {

  // We use the size of bounding boxes and not of m_extSources because the latter
  // could get filled in a second moment
  return m_boundingBoxes.size();

}

std::vector<double>
pyToCppModelInterfaceCache::getExtendedSourceFluxes(int srcid, double j2000_ra, double j2000_dec,
                                               std::vector<double> energies) const {

  SkyCoord sky_pos(j2000_ra, j2000_dec);

  return m_extSources.at(srcid).at(sky_pos);

}

std::vector<double>
pyToCppModelInterfaceCache::getExtendedSourceFluxes_test(int srcid, double j2000_ra, double j2000_dec) const {

  SkyCoord sky_pos(j2000_ra, j2000_dec);

  return m_extSources.at(srcid).at(sky_pos);

}

std::string pyToCppModelInterfaceCache::getPointSourceName(int srcid) const {

  //TODO: implement a mechanism to actually keep the true name

  std::stringstream name;

  name << "Point source " << srcid;

  return name.str();
}

std::string pyToCppModelInterfaceCache::getExtendedSourceName(int srcid) const {

  //TODO: implement a mechanism to actually keep the true name

  std::stringstream name;

  name << "Extended source " << srcid;

  return name.str();
}

void pyToCppModelInterfaceCache::getExtendedSourceBoundaries(int srcid, double *j2000_ra_min,
                                                        double *j2000_ra_max,
                                                        double *j2000_dec_min,
                                                        double *j2000_dec_max) const {

  const BoundingBox *this_boundingBox = &m_boundingBoxes.at(srcid);

  *j2000_ra_min = this_boundingBox->lon_min;
  *j2000_ra_max = this_boundingBox->lon_max;
  *j2000_dec_min = this_boundingBox->lat_min;
  *j2000_dec_max = this_boundingBox->lat_max;

}

}

using namespace threeML;
using namespace boost::python;

//This is needed to wrap the interface (i.e., all methods are virtual)
//contained in ModelInterface.h
struct ModelInterfaceCacheWrap: ModelInterface, wrapper<ModelInterface> {
  int getNumberOfPointSources() const { return this->get_override("getNumberOfPointSources")(); }

  void getPointSourcePosition(int srcid, double *j2000_ra, double *j2000_dec) const {
    this->get_override("getPointSourcePosition")();
  }

  std::vector<double> getPointSourceFluxes(int srcid, std::vector<double> energies) const {
    return this->get_override("getPointSourceFluxes")();
  }

  std::string getPointSourceName(int srcid) const { return this->get_override("getPointSourceName")(); }

  int getNumberOfExtendedSources() const { return this->get_override("getNumberOfExtendedSources")(); }

  std::vector<double> getExtendedSourceFluxes(int srcid, double j2000_ra, double j2000_dec,
                                              std::vector<double> energies) const {
    return this->get_override("getExtendedSourceFluxes")();
  }

  std::string getExtendedSourceName(int srcid) const { return this->get_override("getExtendedSourceName")(); }

  bool isInsideAnyExtendedSource(double j2000_ra, double j2000_dec) const {
    return this->get_override("isInsideAnyExtendedSource")();
  }

  void getExtendedSourceBoundaries(int srcid, double *j2000_ra_min,
                                   double *j2000_ra_max,
                                   double *j2000_dec_min,
                                   double *j2000_dec_max) const { this->get_override("getExtendedSourceBoundaries")(); }
};

template<class T>
struct VecToList {
  static PyObject *convert(const std::vector<T> &vec) {
    boost::python::list *l = new boost::python::list();

    for (size_t i = 0; i < vec.size(); i++)
      (*l).append(vec[i]);

    return l->ptr();
  }
};


BOOST_PYTHON_MODULE (pyModelInterfaceCache) {
  //hello
  to_python_converter<std::vector<double, std::allocator<double> >, VecToList<double> >();

  class_<ModelInterfaceCacheWrap, boost::noncopyable>("ModelInterface")
      .def("getNumberOfPointSources", pure_virtual(&ModelInterface::getNumberOfPointSources))
      .def("getPointSourcePosition", pure_virtual(&ModelInterface::getPointSourcePosition))
      .def("getPointSourceFluxes", pure_virtual(&ModelInterface::getPointSourceFluxes))
      .def("getPointSourceName", pure_virtual(&ModelInterface::getPointSourceName))
      .def("getNumberOfExtendedSources", pure_virtual(&ModelInterface::getNumberOfExtendedSources))
      .def("getExtendedSourceFluxes", pure_virtual(&ModelInterface::getExtendedSourceFluxes))
      .def("getExtendedSourceName", pure_virtual(&ModelInterface::getExtendedSourceName))
      .def("isInsideAnyExtendedSource", pure_virtual(&ModelInterface::isInsideAnyExtendedSource))
      .def("getExtendedSourceBoundaries", pure_virtual(&ModelInterface::getExtendedSourceBoundaries));

  numeric::array::set_module_and_type("numpy", "ndarray");

  class_<pyToCppModelInterfaceCache, bases<ModelInterface> >("pyToCppModelInterfaceCache", init< >())
      .def("getNumberOfPointSources", &pyToCppModelInterfaceCache::getNumberOfPointSources)
      .def("getPointSourcePosition", &pyToCppModelInterfaceCache::getPointSourcePosition)
      .def("getPointSourceFluxes", &pyToCppModelInterfaceCache::getPointSourceFluxes)
      .def("getPointSourceName", &pyToCppModelInterfaceCache::getPointSourceName)
      .def("getNumberOfExtendedSources", &pyToCppModelInterfaceCache::getNumberOfExtendedSources)
      .def("getExtendedSourceFluxes", &pyToCppModelInterfaceCache::getExtendedSourceFluxes)
      .def("getExtendedSourceName", &pyToCppModelInterfaceCache::getExtendedSourceName)
      .def("isInsideAnyExtendedSource", &pyToCppModelInterfaceCache::isInsideAnyExtendedSource)
      .def("getExtendedSourceBoundaries", &pyToCppModelInterfaceCache::getExtendedSourceBoundaries)
      .def("setExtSourceBoundaries", &pyToCppModelInterfaceCache::setExtSourceBoundaries)
      .def("setExtSourceCube",&pyToCppModelInterfaceCache::setExtSourceCube)
      .def("getExtendedSourceFluxes_test",&pyToCppModelInterfaceCache::getExtendedSourceFluxes_test)
      .def("reset",&pyToCppModelInterfaceCache::reset)
      .def("setPtsSourceSpectrum",&pyToCppModelInterfaceCache::setPtsSourceSpectrum)
      .def("setPtsSourcePosition",&pyToCppModelInterfaceCache::setPtsSourcePosition);
}
