#ifndef SPATIALTISSUE_H
#define SPATIALTISSUE_H
#include <fstream>
#include <iostream>
using namespace std;
#define BLOCK_SIZE_X	16		
#define BLOCK_SIZE_Y	8		
#define pi 3.141592653589793
#define R 8314.0	 
#define Frdy 96485.0 
#define Temp 310.0	 
#define FoRT (Frdy / R / Temp)
#define f_Ina 0.5
#define f_k1 0.9
#include "myocyte.cuh"
#include "Cell_Grid.cuh"
template <unsigned int NMYOS, unsigned int NCable>
class SpatialTissue
{
public:
	dim3 dimBlock_myo, dimGrid_myo, dimBlock_grid, dimGrid_grid;
	dim3 dimBlock_cleft, dimGrid_cleft, dimBlock_matrix, dimGrid_matrix;
	cudaStream_t MyoStream, GridStream;
	cudaStream_t memcpyStream_myo, memcpyStream_grid;
	int I_stim_location;
	int I_light_myo_location;
	myocyte<NMYOS, NCable>* rabbit_myocyte;
	myocyte<NMYOS, NCable>* d_rabbit_myocyte;
	Cell_Grid<NMYOS, NCable>* cable_grid;
	Cell_Grid<NMYOS, NCable>* d_cable_grid;
	double radius;
	double length;
	double A_m;
	double A_myo_j;
	double C_m;
	double V_Cell[NMYOS * NCable];
	double A_cell_lat[NMYOS * NCable];
	double A_cell_junc[2][NMYOS * NCable];
	double C_cell_lat[NMYOS * NCable];
	double C_cell_junc[2][NMYOS * NCable];
	double NUMofCells_in_grid[NMYOS * NCable]; 
	double NUMofCell_connected[2][NMYOS * NCable]; 
	double ryo_f_Ina_lat[NMYOS * NCable];
	double ryo_f_k1_lat[NMYOS * NCable];
	double ryo_f_Ina_junc[NMYOS * NCable];
	double ryo_f_k1_junc[NMYOS * NCable];
	int size_of_para_matrix; 
	SpatialTissue(int current_stim_location, int myo_light_location);
	~SpatialTissue();
	void arrange_cells_surface_and_Capacitance(double &rate_af_am, double &rate_aj_af);
	void arrange_cleft_parameter(double *g_cleft_control, double *width_of_cleft_myo, double *G_gap_myo);
	void deliver_parameter_to_class();
	void cpu_copy_gpu();
	void calculate_invert_armadillo();
	void calculate_invert_cublas();
	void invert(double** src, double** dst, int n, int batchSize);
	bool step(double dt, double &Istim, double &Ilight, double &Ilight_myo, double MAXDVDT, bool output_time_check,double in_time);
	void calculate_matrix();
};
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Compute_Myocyte_ODE(double dt, myocyte<NMYOS, NCable>* rabbit_myocyte,Cell_Grid<NMYOS, NCable>* cable_grid, double Ilight_myo, double Istim, int I_light_myo_location, int I_stim_location, double in_time);
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Compute_Cleft_ODE(Cell_Grid<NMYOS, NCable>* cable_grid);
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Arrange_Linear_Constant(Cell_Grid<NMYOS, NCable>* cable_grid);
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Calculate_Linear_dvdt(Cell_Grid<NMYOS, NCable>* cable_grid);
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Assignment_myocyte_dvdt(myocyte<NMYOS, NCable>* rabbit_myocyte, Cell_Grid<NMYOS, NCable>* cable_grid);
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Assignment_Cleft_dvdt(Cell_Grid<NMYOS, NCable>* cable_grid);
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Myocyte_Update(double dt, myocyte<NMYOS, NCable>* rabbit_myocyte){
	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int id  = idy * NMYOS + idx;
	if (idx >= NMYOS || idy >= NCable)
		return;
	rabbit_myocyte -> rabbit_update(id, dt);
}
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Cleft_update(double dt, Cell_Grid<NMYOS, NCable>* cable_grid){
	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int id  = idy * (NMYOS - 1) + idx;
	if (idx >= (NMYOS - 1) || idy >= NCable)
		return;
	cable_grid -> compute_dion_dt(idx, idy);
	cable_grid -> cleft_update(id, dt);
}
#include "SpatialTissue.cu"
#endif
