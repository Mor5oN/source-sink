#ifndef MYOCYTE_H
#define MYOCYTE_H
#define CKIIOE 0 
#define plb_val 106 
#define LCCtotDyad 31.4 * 0.9 
#define LCCtotSL 0.0846		  
#define RyRtot 382.6		  
#define PP1_dyad 95.7		  
#define PP1_SL 0.57			  
#define PP2A_dyad 95.76		  
#define OA 0.0				  
#define PLBtot plb_val		  
#define Ligtot 0.0		 
#define LCCtotBA 0.025	 
#define RyRtotBA 0.135	 
#define PLBtotBA plb_val 
#define TnItotBA 70.0	 
#define IKstotBA 0.025	 
#define ICFTRtotBA 0.025 
#define PP1_PLBtot 0.89	 
#define IKurtotBA 0.025	 
#define PLMtotBA 48.0	 
#define Qpow ((Temp - 310.0) / 10.0)
template <unsigned int NMYOS, unsigned int NCable>
class myocyte
{
public:
	double y[181][NMYOS * NCable];	  
	int myo_location[NMYOS * NCable];
	int myo_cable_location[NMYOS * NCable];
	int myo_cable_id[NMYOS * NCable];
	double A_cell_lat[NMYOS * NCable];
	double A_cell_junc[2][NMYOS * NCable];
	double C_cell_lat[NMYOS * NCable];
	double C_cell_junc[2][NMYOS * NCable];
	double NUMofCell_connected[2][NMYOS * NCable];
	double ryo_f_Ina_lat[NMYOS * NCable];
	double ryo_f_k1_lat[NMYOS * NCable];
	double ryo_f_Ina_junc[NMYOS * NCable];
	double ryo_f_k1_junc[NMYOS * NCable];
	double O1_ChR2_myo[NMYOS * NCable];
	double O2_ChR2_myo[NMYOS * NCable];
	double C1_ChR2_myo[NMYOS * NCable];
	double C2_ChR2_myo[NMYOS * NCable];
	double p_ChR2_myo[NMYOS * NCable];
	double ydot[181][NMYOS * NCable];
	double d_O1_ChR2_myo[NMYOS * NCable];
	double d_O2_ChR2_myo[NMYOS * NCable];
	double d_C1_ChR2_myo[NMYOS * NCable];
	double d_C2_ChR2_myo[NMYOS * NCable];
	double d_p_ChR2_myo[NMYOS * NCable];
	double i_Na_lateral[NMYOS * NCable];
	double i_K1_lateral[NMYOS * NCable];
	double i_comp_lateral[NMYOS * NCable];
	double i_ChR2_myo[NMYOS * NCable];
	double I_tot_lateral[NMYOS * NCable];
	double I_Na_lateral_foroutput[NMYOS * NCable];
	double I_K1_lateral_foroutput[NMYOS * NCable];
	double m_myo_junc[2][NMYOS * NCable];
	double h_myo_junc[2][NMYOS * NCable];
	double j_myo_junc[2][NMYOS * NCable];
	double V_myo_junc[2][NMYOS * NCable];
	double Nao[2][NMYOS * NCable];
	double Ko[2][NMYOS * NCable];
	double d_m_myo_junc[2][NMYOS * NCable];
	double d_h_myo_junc[2][NMYOS * NCable];
	double d_j_myo_junc[2][NMYOS * NCable];
	double i_Na_myo_cl[2][NMYOS * NCable];
	double i_K1_myo_cl[2][NMYOS * NCable];
	double I_tot_myo_cl[2][NMYOS * NCable];
	double I_Na_tot_myo_cl[2][NMYOS * NCable];
	double I_K1_tot_myo_cl[2][NMYOS * NCable];
	__host__ myocyte();
	__host__ __device__ void eccODEfile(int id, double dt, double &Ilight_myo0,double in_time);
	__host__ __device__ void DAD_eccODEfile(int id, double dt, double &Ilight_myo0,double in_time);
	__host__ __device__ void rabbit_juncODEfile(int id, int j, double dt);
	__host__ __device__ void rabbit_covert_to_Current(int id, double Istim);
	__host__ __device__ void rabbit_update(int id, double dt);
	__host__ __device__ void i_Na_lat(int id, double dt, double &ena_junc, double &ena_sl, double &Fjunc, double &Fsl);
	__host__ __device__ void i_K1_lat(int id, double &Ko, double &ek, double &Kcoeff);
	__host__ __device__ double i_NaL_lat(int id, double dt, double &ena_junc, double &ena_sl, double &Fjunc, double &Fsl);
	__host__ __device__ double i_Nabk_lat(int id, double &ena_junc, double &ena_sl, double &Fjunc, double &Fsl);
	__host__ __device__ double i_NaK_lat(int id, double &Nao, double &Ko, double &Fjunc_nak, double &Fsl_nak, double &PLM_PKAp);
	__host__ __device__ double i_Kur_lat(int id, double dt, double &ek, double &Kcoeff, double &IKur_PKAp);
	__host__ __device__ double i_ss_lat(int id, double dt, double &ek, double &Kcoeff);
	__host__ __device__ double i_Kr_lat(int id, double dt, double &Ko, double &ek, double &Kcoeff);
	__host__ __device__ double i_Kp_lat(int id, double &ek, double &Fjunc, double &Fsl, double &Kcoeff);
	__host__ __device__ double i_to_lat(int id, double dt, double &ek, double &Kcoeff);
	__host__ __device__ void i_Cl_lat(int id, double &ecl, double &Fjunc, double &Fsl, double &I_ClCa, double &I_Clbk);
	__host__ __device__ void i_LCC_HH_lat(int id, double dt, double &Nao, double &Ko, double &Cao, double &Fjunc_CaL, double &Fsl_CaL, double &I_Ca_junc, double &I_Ca_sl, double &I_CaNa_junc, double &I_CaNa_sl, double &I_CaK);
	__host__ __device__ void i_LCC_Markov_lat(int id, double dt, double &Nao, double &Ko, double &Cao, double &Fjunc_CaL, double &Fsl_CaL, double &I_Ca_junc, double &I_Ca_sl, double &I_CaNa_junc, double &I_CaNa_sl, double &I_CaK);
	__host__ __device__ double i_NCX_lat(int id, double &Nao, double &Cao, double &Fjunc_ncx, double &Fsl_ncx, double &I_ncx_junc, double &I_ncx_sl);
	__host__ __device__ double i_pca_lat(int id, double &Fjunc, double &Fsl, double &I_pca_junc, double &I_pca_sl);
	__host__ __device__ double i_Cabk_lat(int id, double &Fjunc, double &Fsl, double &I_cabk_junc, double &I_cabk_sl, double &eca_junc, double &eca_sl);
	__host__ __device__ void i_RyR_lat(int id, double &J_SRCarel, double &J_SRleak);
	__host__ __device__ double i_SERCA_lat(int id);
	__host__ __device__ void Myofilament(int id, int &mechFlag);
	__host__ __device__ void Na_Ca_buffer(int id, double &Mgi, double &Vmyo, double &Vjunc, double &Vsl, double &J_CaB_cytosol, double &J_CaB_junction, double &J_CaB_sl);
	__host__ __device__ double i_ChR2(int id,double &Ilight);
};
#include "myocyte.cu"
#endif
