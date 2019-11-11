#include <iostream>
#include "class.h"
#include "H5Cpp.h"

void dumpMultipolesH5(std::string fname, size_t Nl, size_t Nq, size_t Nk, float *psi_arr, float *qcol, float *kcol);

/**
* dump_ncdm_perturb.cpp
* Copyright 2019 Willem Elbers <whe@willemelbers.com>
*
* This program runs CLASS as expected, but - assuming that
* ncdm species are evolved - dumps the full perturbation vector
* around redshift 40 (this is hard-coded for now but easily changed)
* into an HDF5 file. Specifically, the output is the Legendre coefficient
* $\Psi_l(k,q,tau)$ with tau corresponding to $z=40$ from the Boltzmann
* expansion. Here, we use the standard notation of Ma & Bertschinger (1995).
*/

int main(int argc, char **argv) {
	struct precision pr;        /* for precision parameters */
	struct background ba;       /* for cosmological background */
	struct thermo th;           /* for thermodynamics */
	struct perturbs pt;         /* for source functions */
	struct transfers tr;        /* for transfer functions */
	struct primordial pm;       /* for primordial spectra */
	struct spectra sp;          /* for output spectra */
	struct nonlinear nl;        /* for non-linear spectra */
	struct lensing le;          /* for lensed spectra */
	struct output op;           /* for output files */
	ErrorMsg errmsg;            /* for error messages */

	if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
	  printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
	  return _FAILURE_;
	}

	if (background_init(&pr,&ba) == _FAILURE_) {
      printf("\n\nError running background_init \n=>%s\n",ba.error_message);
      return _FAILURE_;
    }

    if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
      printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
      return _FAILURE_;
    }

    if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
      printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
      return _FAILURE_;
    }

	/** WHE */

	int Nn = ba.N_ncdm;
	int Nl = pr.l_max_ncdm;
	int Nq = pt.max_q_size_ncdm;
	int Nk = pt.k_size[pt.index_md_scalars];

	std::cout << "-> Exporting ncdm perturbations (one file per species)..." << std::endl;
	std::cout << "  ncdm dump redshift: " << pt.ncdm_dump_redshift << std::endl;
	std::cout << "  ncdm species: " << Nn << std::endl;
	std::cout << "  multipoles: " << Nl << std::endl;
	std::cout << "  momentum bins: " << Nq << std::endl;
	std::cout << "  wavenumbers: " << Nk << std::endl;

	//Get the wavenumbers
	float kcol[Nk];
	for (int ik=0; ik<Nk; ik++) {
		kcol[ik] = pt.k[pt.index_md_scalars][ik];
	}

	//For each ncdm species
	for (int in=0; in<Nn; in++) {
		//Get the multipoles
		float Psi[Nl*Nq*Nk];
		for (int ik=0; ik<Nk; ik++) {
			for (int iq=0; iq<Nq; iq++) {
				for (int il=0; il<Nl; il++) {
					Psi[il + iq*Nl + ik*Nl*Nq] = (float) pt.ncdm_y[il + iq*Nl + ik*Nl*Nq + in*Nl*Nq*Nk];
				}
			}
		}

		//And get the momentum bins
		int Nq_this = ba.q_size_ncdm[in];
		float qcol[Nq_this];
		for (int iq=0; iq<Nq_this; iq++) {
			qcol[iq] = ba.q_ncdm[in][iq];
		}

		dumpMultipolesH5("multipoles_" + std::to_string(in) + ".h5", Nl, Nq, Nk, Psi, qcol, kcol);
		std::cout << "  Dumped to " << "multipoles_" + std::to_string(in) + ".h5" << std::endl;
	}

	/** end WHE */

    if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
      printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
      return _FAILURE_;
    }

    if (nonlinear_init(&pr,&ba,&th,&pt,&pm,&nl) == _FAILURE_) {
      printf("\n\nError in nonlinear_init \n=>%s\n",nl.error_message);
      return _FAILURE_;
    }

    if (transfer_init(&pr,&ba,&th,&pt,&nl,&tr) == _FAILURE_) {
      printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
      return _FAILURE_;
    }

    if (spectra_init(&pr,&ba,&pt,&pm,&nl,&tr,&sp) == _FAILURE_) {
      printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
      return _FAILURE_;
    }

    if (lensing_init(&pr,&pt,&sp,&nl,&le) == _FAILURE_) {
      printf("\n\nError in lensing_init \n=>%s\n",le.error_message);
      return _FAILURE_;
    }

    if (output_init(&ba,&th,&pt,&pm,&tr,&sp,&nl,&le,&op) == _FAILURE_) {
      printf("\n\nError in output_init \n=>%s\n",op.error_message);
      return _FAILURE_;
    }



    /****** all calculations done, now free the structures ******/

    if (lensing_free(&le) == _FAILURE_) {
      printf("\n\nError in lensing_free \n=>%s\n",le.error_message);
      return _FAILURE_;
    }

    if (spectra_free(&sp) == _FAILURE_) {
      printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
      return _FAILURE_;
    }

    if (transfer_free(&tr) == _FAILURE_) {
      printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
      return _FAILURE_;
    }

    if (nonlinear_free(&nl) == _FAILURE_) {
      printf("\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
      return _FAILURE_;
    }

    if (primordial_free(&pm) == _FAILURE_) {
      printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
      return _FAILURE_;
    }

    if (perturb_free(&pt) == _FAILURE_) {
      printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
      return _FAILURE_;
    }

    if (thermodynamics_free(&th) == _FAILURE_) {
      printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
      return _FAILURE_;
    }

    if (background_free(&ba) == _FAILURE_) {
      printf("\n\nError in background_free \n=>%s\n",ba.error_message);
      return _FAILURE_;
    }

    return _SUCCESS_;
}

/** WHE */

/**
* Write the Legendre coefficients of the neutrino perturbation to the disk.
* The array is 3D: Arr[k,q,l] = \Psi_l(k,q) a function of wavenumber
* magnitude of the momentum vector, and multipole
*/
void dumpMultipolesH5(std::string fname, size_t Nl, size_t Nq, size_t Nk, float *psi_arr, float *qcol, float *kcol) {
	H5::H5File file(fname, H5F_ACC_TRUNC);

	//Write the 3D array to disk
	std::size_t NDIMS = 3;
	hsize_t dims[NDIMS] = {Nk, Nq, Nl};
	H5::DataSpace dataspace(NDIMS, dims);
	H5::DataSet dataset = file.createDataSet("Psi(k,q,l)", H5::PredType::NATIVE_FLOAT, dataspace);
	dataset.write(psi_arr, H5::PredType::NATIVE_FLOAT);

	//Write a description of the dimensions
	// size_t string_length = 37;
	// H5::StrType strtype(H5::PredType::C_S1, string_length);
	// std::vector<hsize_t> attributes = {1};
	// H5::DataSpace attr_space(1, attributes.data());
	// H5::Attribute attr = dataset.createAttribute("dimensions", strtype, attr_space);
	// attr.write(strtype, "wavenumber k, momentum q, multipole l");

	//We also write the momentum and wavenumber bins to disk
	hsize_t q_dimensions[1] = {Nq};
	hsize_t k_dimensions[1] = {Nk};
	H5::DataSpace q_dataspace(1, q_dimensions);
	H5::DataSpace k_dataspace(1, k_dimensions);
	H5::DataSet q_dataset = file.createDataSet("momentum_bins (q)", H5::PredType::NATIVE_FLOAT, q_dataspace);
	H5::DataSet k_dataset = file.createDataSet("wavenumber_bins (k)", H5::PredType::NATIVE_FLOAT, k_dataspace);

	q_dataset.write(qcol, H5::PredType::NATIVE_FLOAT);
	k_dataset.write(kcol, H5::PredType::NATIVE_FLOAT);
}

/** end WHE */
