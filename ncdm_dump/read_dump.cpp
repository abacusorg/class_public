#include <iostream>
#include "H5Cpp.h"

int main(int argc, char **argv) {
	//Open the HDF5 file
	H5::H5File f("multipoles_0.h5", H5F_ACC_RDONLY);

	//Open the three dataset (psi, wavenumbers, momentum bins)
	H5::DataSet dataset = f.openDataSet("Psi(k,q,l)");
	H5::DataSet k_dataset = f.openDataSet("wavenumber_bins (k)");
	H5::DataSet q_dataset = f.openDataSet("momentum_bins (q)");

	//Check that the data types are as expected
	H5T_class_t type = dataset.getTypeClass();
	H5T_class_t k_type = k_dataset.getTypeClass();
	H5T_class_t q_type = q_dataset.getTypeClass();

	//We expect a dataset consisting of floating point numbers
	if(type != H5T_FLOAT || k_type != H5T_FLOAT || q_type != H5T_FLOAT) {
		throw std::logic_error("Invalid data type in HDF5 file");
	}

	//We therefore only need one FloatType
	H5::FloatType flt = dataset.getFloatType();

	//Load the dataspaces
	H5::DataSpace dataspace = dataset.getSpace();
	H5::DataSpace k_dataspace = k_dataset.getSpace();
	H5::DataSpace q_dataspace = q_dataset.getSpace();

	//Load the psi dataspace and determine the dimensions
	int NDIMS = dataspace.getSimpleExtentNdims();
	hsize_t dimensions[NDIMS];
	dataspace.getSimpleExtentDims(dimensions, NULL);

	//The physical dimensions
	int Nk = dimensions[0]; //number of wavenumbers (k)
	int Nq = dimensions[1]; //number of momentum bins (q)
	int Nl = dimensions[2]; //number of multipoles (l)

	//Do the dimensions check out?
	hsize_t k_dimensions[1];
	hsize_t q_dimensions[1];
	k_dataspace.getSimpleExtentDims(k_dimensions, NULL);
	q_dataspace.getSimpleExtentDims(q_dimensions, NULL);
	if (k_dimensions[0] != Nk || q_dimensions[0] != Nq) {
		throw std::logic_error("Invalid dataset dimensions in HDF5 file");
	}

	//The data arrays
	float Psi[Nk*Nq*Nl];
	float kcol[Nk];
	float qcol[Nq];

	//Read the data into the array
	dataset.read(Psi, H5::PredType::NATIVE_FLOAT, dataspace);
	k_dataset.read(kcol, H5::PredType::NATIVE_FLOAT, k_dataspace);
	q_dataset.read(qcol, H5::PredType::NATIVE_FLOAT, q_dataspace);

	for (int i=0; i<Nq; i++) {
		std::cout << i << " " << qcol[i] << std::endl;
	}

	for (int i=0; i<Nk; i++) {
		std::cout << i << " " << kcol[i] << std::endl;
	}

	std::cout << "say wat " << H5T_INTEGER << std::endl;
}
