#ifndef MediumLoader_H
#define MediumLoader_H

#include <cmath>
#include <vector>
#include <string>
#include <H5Cpp.h>
#include <boost/multi_array.hpp>
#include <sstream>
#include <iomanip>
#include <iostream>

template <typename T> inline const H5::PredType& type();
template <> inline const H5::PredType& type<size_t>() { return H5::PredType::NATIVE_HSIZE; }
template <> inline const H5::PredType& type<double>() { return H5::PredType::NATIVE_DOUBLE; }
template <> inline const H5::PredType& type<int>() { return H5::PredType::NATIVE_INT; }

template <typename T>
void hdf5_add_scalar_attr(
  const H5::Group& gp, const std::string& name, const T& value) {
  const auto& datatype = type<T>();
  auto attr = gp.createAttribute(name.c_str(), datatype, H5::DataSpace{});
  attr.write(datatype, &value);
}

template <typename T>
void hdf5_read_scalar_attr(
  const H5::Group& gp, const std::string& name, T& value) {
  const auto& datatype = type<T>();
  auto attr = gp.openAttribute(name.c_str());
  attr.read(datatype, &value);
}

class MediumProfile {
  private:
    const double fmc_to_GeV_m1;
    const size_t N, _power_rank;
    std::vector<size_t> _shape, _spatial_dims;
    std::vector<double> _xl_limits, _xh_limits, _x_steps;
    boost::multi_array<double, 3> _Tgrid;
    boost::multi_array<double, 2> _buffer;
    double _dtau, _tau0;
    double _dx, _xl, _xh, _dy, _yl, _yh;
    int _Nt, _Nx, _Ny, _iXL, _iXH,_iYL, _iYH;
  public:
    MediumProfile();
    bool ReadHydroGridFile(std::string fname);
    bool GetTemp(double t, double x, double y, double & T);
};

class TABProfile {
  private:
    const double fmc_to_GeV_m1;
    const size_t N, _power_rank;
    std::vector<size_t> _shape;
    std::vector<double> _xl_limits, _xh_limits, _x_steps;
    boost::multi_array<double, 2> _TABgrid;
    double _dx, _dy, _x_min, _y_min, _x_max, _y_max;
    int _Nx, _Ny;
    double IntegratedTAB;
  public:
    TABProfile();
    bool ReadTABGridFile(std::string fname);
    double GetTAB(double x, double y);
    double GetIntegratedTAB() {return IntegratedTAB;};
};



#endif
