//
// Created by giacomov on 11/22/17.
//

#ifndef CTHREEML_BOOST_NUMPY_H_H
#define CTHREEML_BOOST_NUMPY_H_H

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

// Boost python numpy available in Boost 1.63+
// Boost python numeric removed in Boost 1.65+
#if BOOST_VERSION < 106500
#include <boost/python/numeric.hpp>
typedef typename boost::python::numeric::array NumpyArrayType;
#else
#include <boost/python/numpy.hpp>
typedef typename boost::python::numpy::ndarray NumpyArrayType;
#endif

#endif //CTHREEML_BOOST_NUMPY_H_H
