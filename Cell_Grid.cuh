#ifndef CELL_GRID_H
#define CELL_GRID_H
template <unsigned int NMYOS, unsigned int NCable>
class Cell_Grid
{
public:
    double V_Cell[NMYOS * NCable];
	double Cell_lat_current[NMYOS * NCable];
	double Cell_junc_current[2][NMYOS * NCable];           
	double Cell_Na_junc_current[2][NMYOS * NCable];        
	double Cell_K1_junc_current[2][NMYOS * NCable];
	double Cell_lat_Ina_current_foroutput[NMYOS * NCable];
	double Cell_lat_IK1_current_foroutput[NMYOS * NCable];
	double NUMofCells_in_grid[NMYOS * NCable]; 
	double NUMofCell_connected[2][NMYOS * NCable];         
	double V_cleft[(NMYOS - 1) * NCable];
	double Na_cleft[(NMYOS - 1) * NCable];
	double K_cleft[(NMYOS - 1) * NCable];
	double d_V_cleft[(NMYOS - 1) * NCable];
	double d_Na_cleft[(NMYOS - 1) * NCable];
	double d_K_cleft[(NMYOS - 1) * NCable];
	double I_Na_cleft[(NMYOS - 1) * NCable];
	double I_K_cleft[(NMYOS - 1) * NCable];
	double I_tot_cleft[(NMYOS - 1) * NCable];
	double I_Na_tot_left_cl[(NMYOS - 1) * NCable];          
	double I_Na_tot_right_cl[(NMYOS - 1) * NCable];
	double I_K_tot_left_cl[(NMYOS - 1) * NCable];
	double I_K_tot_right_cl[(NMYOS - 1) * NCable];
	double g_cleft[(NMYOS - 1) * NCable];                   
	double width[(NMYOS - 1) * NCable];
	double vcl[(NMYOS - 1) * NCable];
	double G_gap[(NMYOS - 1) * NCable];                     
	double I_gap[(NMYOS - 1) * NCable];
	int size_of_para_matrix;
	double inv_matrix[NCable][(NMYOS + NMYOS - 1)*(NMYOS + NMYOS - 1)];
	double constant[(NMYOS + NMYOS - 1) * NCable];
	double solution[(NMYOS + NMYOS - 1) * NCable];
    __host__ Cell_Grid();
	__device__ int myFloor(double x) {
		if (x >= 0.0) {
			return floor(x);
		} else {
			return ceil(x) - 1.0;
		}
	}
	__host__ __device__ void compute_I_cleft(int id);
	__host__ __device__ void compute_dion_dt(int idx, int idy);
	__host__ __device__ void cleft_update(int id, double dt);
	__host__ __device__ void arrange_constant(int idx, int idy);
};
#include"Cell_Grid.cu"
#endif
