#include "MediumLoader.h"


MediumProfile::MediumProfile()
:fmc_to_GeV_m1(5.068), N(3), _power_rank(std::pow(2, N)){
}

bool MediumProfile::ReadHydroGridFile(std::string filename){
    
    // grid size information
    H5::H5File _file(filename, H5F_ACC_RDONLY);
    H5::Group _event(_file.openGroup("/Event"));
    hdf5_read_scalar_attr(_event, "Tau0", _tau0);
    hdf5_read_scalar_attr(_event, "dTau", _dtau);
    hdf5_read_scalar_attr(_event, "DX", _dx);
    hdf5_read_scalar_attr(_event, "DY", _dy);
    hdf5_read_scalar_attr(_event, "XL", _iXL);
    hdf5_read_scalar_attr(_event, "XH", _iXH);
    hdf5_read_scalar_attr(_event, "YL", _iYL);
    hdf5_read_scalar_attr(_event, "YH", _iYH);

    hsize_t  num_obj;
    H5Gget_num_objs(_event.getId(), &num_obj);
    _Nt = static_cast<int>(num_obj);
    
    _tau0 *=  fmc_to_GeV_m1;
    _dtau *=  fmc_to_GeV_m1; 
    _dx *=  fmc_to_GeV_m1;
    _dy *=  fmc_to_GeV_m1;
    _xl = _iXL * _dx;
    _xh = _iXH * _dx;
    _yl = _iYL * _dy;
    _yh = _iYH * _dy;
    

    
    _shape.resize(N);
    _shape[0] = _Nt; 
    _shape[1] = 2*_iXH+1; 
    _shape[2] =  2*_iYH+1;
    _spatial_dims.resize(N-1);
    for (auto i=0; i<N-1; ++i) _spatial_dims[i] = _shape[i+1]; 
    
    
    std::cout << "Grid size " << _shape[0] << " " << _shape[1] << " " << _shape[2] << std::endl;
    _Tgrid.resize(_shape);
    _buffer.resize(_spatial_dims);
    
    _xl_limits.resize(N); 
    _xh_limits.resize(N); 
    _x_steps.resize(N);

    _xl_limits[0] = _tau0; 
    _xh_limits[0] = _tau0 + (_Nt-1)*_dtau;
    _x_steps[0] = _dtau;

    _xl_limits[1] = _xl; 
    _xh_limits[1] = _xh;
    _x_steps[1] = _dx;

    _xl_limits[2] = _yl; 
    _xh_limits[2] = _yh; 
    _x_steps[2] = _dy; 
    
    
    // Load medium temperature grid
    H5::DataSet dataset;  
    hsize_t h5_spatial_dims[N]; 
    for (auto i=0; i<N-1; ++i) h5_spatial_dims[i]=_spatial_dims[i];   
    H5::DataSpace spatial_dataspace(N-1, h5_spatial_dims);
    for (int iTau=0; iTau<_Nt; ++iTau) {
        int istart = _buffer.num_elements()*iTau;
        std::stringstream FrameNumber;
        FrameNumber << std::setw(4) << std::setfill('0') << iTau;
        
        std::string FrameName = "/Event/Frame_"+FrameNumber.str();
        //std::cout << "reading " << FrameName << std::endl;
        dataset = _file.openDataSet(FrameName+"/Temp");
        dataset.read(_buffer.data(), H5::PredType::NATIVE_DOUBLE,
                     spatial_dataspace, dataset.getSpace());   
        // write it to the 3D array
        for(int i=0; i<_buffer.num_elements(); ++i) 
            _Tgrid.data()[i+istart] = _buffer.data()[i];
    }
    _file.close();
    return true;
}


bool MediumProfile::GetTemp(double t, double x, double y, double & T){

    double X[N]; X[0] = t; X[1] = x; X[2] = y;
    // check limits first
    for(int i=0; i<N; ++i) {
        if ( X[i] < _xl_limits[i] || X[i] > _xh_limits[i] ){
            T = 0.;
            return true;
        }
    }
   
    std::vector<size_t> start_index(N);
    std::vector<double> w(N);
    for(int i=0; i<N; ++i) {
        double var = (X[i]-_xl_limits[i])/_x_steps[i];
        var = std::min(std::max(var, 0.), _shape[i]-2.);
        size_t n = size_t(floor(var));
        double rx = var-n;
        w[i] = rx;
        start_index[i] = n;
    }

    std::vector<size_t> index(N);
    T = 0.;
    
    for(int i=0; i<_power_rank; ++i) {
        double W = 1.0;
        // loop over the 2^N corner of the hypercube 
        for (int j=0; j<N; ++j) {
            index[j] = start_index[j] + ((i & ( 1 << j )) >> j);
            W *= (index[j]==start_index[j])?(1.-w[j]):w[j];
        }
        // and perform linear interpolation (weighted average of the values at the corner)
        T = T + _Tgrid(index)*W;
    }
    return true;
}


TABProfile::TABProfile()
:fmc_to_GeV_m1(5.068), N(2), _power_rank(std::pow(2, N)){
}

bool TABProfile::ReadTABGridFile(std::string filename) {
    H5::H5File _file(filename, H5F_ACC_RDONLY);
    H5::Group _event(_file.openGroup("/event_0"));
    
    hdf5_read_scalar_attr(_event, "Nx", _Nx);
    hdf5_read_scalar_attr(_event, "Ny", _Ny);
    hdf5_read_scalar_attr(_event, "dxy", _dx);
    hdf5_read_scalar_attr(_event, "dxy", _dy);
    _dx *=  fmc_to_GeV_m1;
    _dy *=  fmc_to_GeV_m1;
    _x_min = -0.5*_Nx*_dx;
    _y_min = -0.5*_Ny*_dy;
    _x_max = 0.5*_Nx*_dx;
    _y_max = 0.5*_Ny*_dy;
    
    
    _xl_limits.resize(2);
    _xh_limits.resize(2);
    _x_steps.resize(2);
    
    _xl_limits[0] = _x_min; 
    _xh_limits[0] = _x_max;
    _x_steps[0] = _dx;

    _xl_limits[1] = _y_min; 
    _xh_limits[1] = _y_max; 
    _x_steps[1] = _dy; 
    
    
    _shape.resize(2);
    _shape[0] = _Nx;
    _shape[1] = _Ny;
    _TABgrid.resize(_shape);
    
    hsize_t dims[2]; dims[0] = _Nx; dims[1] = _Ny;
    H5::DataSet dataset;    
    H5::DataSpace dataspace(2, dims);
    
    dataset = _event.openDataSet("Ncoll_density");
    dataset.read(_TABgrid.data(), H5::PredType::NATIVE_DOUBLE, 
                 dataspace, dataset.getSpace());
                 
    // Integrated TAB (the number of binary coll, as is normalized in TRENTo)
    IntegratedTAB = 0.;
    for (int i=0; i<_Nx; i++){
        for (int j=0; j<_Ny; j++){
            IntegratedTAB += _TABgrid[i][j];
        }
    }
    IntegratedTAB *= _dx*_dy;
    
    return true;
}



double TABProfile::GetTAB(double x, double y){

    double X[N]; X[0] = x; X[1] = y;
    // check limits first
    for(int i=0; i<N; ++i) {
        if ( X[i] < _xl_limits[i] || X[i] > _xh_limits[i] )
            return 0.;
    }
   
    std::vector<size_t> start_index(N);
    std::vector<double> w(N);
    for(int i=0; i<N; ++i) {
        double var = (X[i]-_xl_limits[i])/_x_steps[i];
        var = std::min(std::max(var, 0.), _shape[i]-2.);
        size_t n = size_t(floor(var));
        double rx = var-n;
        w[i] = rx;
        start_index[i] = n;
    }

    std::vector<size_t> index(N);
    double res = 0.;
    
    for(int i=0; i<_power_rank; ++i) {
        double W = 1.0;
        // loop over the 2^N corner of the hypercube 
        for (int j=0; j<N; ++j) {
            index[j] = start_index[j] + ((i & ( 1 << j )) >> j);
            W *= (index[j]==start_index[j])?(1.-w[j]):w[j];
        }
        // and perform linear interpolation (weighted average of the values at the corner)
        res = res + _TABgrid(index)*W;
    }
    return res;
}


