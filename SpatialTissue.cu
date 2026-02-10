#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
	if (code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}
#define cudacall(call)                                                                                                          \
    do                                                                                                                          \
    {                                                                                                                           \
        cudaError_t err = (call);                                                                                               \
        if(cudaSuccess != err)                                                                                                  \
        {                                                                                                                       \
            fprintf(stderr,"CUDA Error:\nFile = %s\nLine = %d\nReason = %s\n", __FILE__, __LINE__, cudaGetErrorString(err));    \
            cudaDeviceReset();                                                                                                  \
            exit(EXIT_FAILURE);                                                                                                 \
        }                                                                                                                       \
    }                                                                                                                           \
    while (0)
#define cublascall(call)                                                                                        \
    do                                                                                                          \
    {                                                                                                           \
        cublasStatus_t status = (call);                                                                         \
        if(CUBLAS_STATUS_SUCCESS != status)                                                                     \
        {                                                                                                       \
            fprintf(stderr,"CUBLAS Error:\nFile = %s\nLine = %d\nCode = %d\n", __FILE__, __LINE__, status);     \
            cudaDeviceReset();                                                                                  \
            exit(EXIT_FAILURE);                                                                                 \
        }                                                                                                       \
                                                                                                                \
    }                                                                                                           \
    while(0)
template <unsigned int NMYOS, unsigned int NCable>
SpatialTissue<NMYOS, NCable>::SpatialTissue(int current_stim_location)
{
	I_stim_location = current_stim_location;
	size_of_para_matrix = NMYOS + NMYOS - 1; 
	radius = 15.0;	
	length = 100.0; 
	dimBlock_myo.x = BLOCK_SIZE_X;
	dimBlock_grid.x = BLOCK_SIZE_X;
	dimBlock_cleft.x = BLOCK_SIZE_X;
	dimBlock_matrix.x = BLOCK_SIZE_X;
	dimBlock_myo.y = BLOCK_SIZE_Y;
	dimBlock_grid.y = BLOCK_SIZE_Y;
	dimBlock_cleft.y = BLOCK_SIZE_Y;
	dimBlock_matrix.y = BLOCK_SIZE_Y;
	dimGrid_myo.x = (NMYOS + dimBlock_myo.x - 1) / dimBlock_myo.x;
	dimGrid_grid.x = (NMYOS + dimBlock_grid.x - 1) / dimBlock_grid.x;
	dimGrid_cleft.x = (NMYOS - 1 + dimBlock_cleft.x - 1) / dimBlock_cleft.x;
	dimGrid_matrix.x = (size_of_para_matrix + dimBlock_matrix.x - 1) / dimBlock_matrix.x;
	dimGrid_myo.y = (NCable + dimBlock_myo.y - 1) / dimBlock_myo.y;
	dimGrid_grid.y = (NCable + dimBlock_grid.y - 1) / dimBlock_grid.y;
	dimGrid_cleft.y = (NCable + dimBlock_cleft.y - 1) / dimBlock_cleft.y;
	dimGrid_matrix.y = (NCable + dimBlock_matrix.y - 1) / dimBlock_matrix.y;
	cudaStreamCreate(&MyoStream);
	cudaStreamCreate(&GridStream);
	cudaStreamCreate(&memcpyStream_myo);
	cudaStreamCreate(&memcpyStream_grid);
	rabbit_myocyte = new myocyte<NMYOS, NCable>();
	cable_grid = new Cell_Grid<NMYOS, NCable>();
	cudaMalloc((void**) &d_rabbit_myocyte, sizeof(myocyte<NMYOS, NCable>));
	cudaMalloc((void**) &d_cable_grid, sizeof(Cell_Grid<NMYOS, NCable>));
	int id;
	for (int idy = 0; idy < NCable; idy++){
		for (int idx = 0; idx < NMYOS; idx++)
		{
			id = idx + NMYOS * idy;
			rabbit_myocyte -> myo_location[id] = id;
			rabbit_myocyte -> myo_cable_location[id] = idx;
			rabbit_myocyte -> myo_cable_id[id] = idy;
			V_Cell[id] = rabbit_myocyte -> y[39][id];
		}
	}
}
template <unsigned int NMYOS, unsigned int NCable>
SpatialTissue<NMYOS, NCable>::~SpatialTissue()
{
	free(rabbit_myocyte);
	free(cable_grid);
	cudaStreamDestroy(MyoStream);
	cudaStreamDestroy(GridStream);
	cudaStreamDestroy(memcpyStream_myo);
	cudaStreamDestroy(memcpyStream_grid);
	cudaFree(d_rabbit_myocyte);
	cudaFree(d_cable_grid);
}
template <unsigned int NMYOS, unsigned int NCable>
bool SpatialTissue<NMYOS, NCable>::step(double dt, double &Istim, double MAXDVDT, bool output_time_check,double in_time)
{
	Compute_Myocyte_ODE<<<dimGrid_myo, dimBlock_myo, 0, MyoStream>>>(dt, d_rabbit_myocyte, d_cable_grid, Istim, I_stim_location,in_time);
	gpuErrchk(cudaDeviceSynchronize());
	Compute_Cleft_ODE<<<dimGrid_cleft, dimBlock_cleft, 0, GridStream>>>(d_cable_grid);
	Arrange_Linear_Constant<<<dimGrid_matrix, dimBlock_matrix, 0, GridStream>>>(d_cable_grid);
	Calculate_Linear_dvdt<<<dimGrid_matrix, dimBlock_matrix, 0, GridStream>>>(d_cable_grid);
	gpuErrchk(cudaDeviceSynchronize());
	if(output_time_check){
		gpuErrchk(cudaMemcpyAsync(rabbit_myocyte, d_rabbit_myocyte, sizeof(myocyte<NMYOS, NCable>), cudaMemcpyDeviceToHost, memcpyStream_myo));
		gpuErrchk(cudaMemcpyAsync(cable_grid, d_cable_grid, sizeof(Cell_Grid<NMYOS, NCable>), cudaMemcpyDeviceToHost, memcpyStream_grid));
	}
	Assignment_myocyte_dvdt<<<dimGrid_matrix, dimBlock_matrix, 0, MyoStream>>>(d_rabbit_myocyte, d_cable_grid);
	Assignment_Cleft_dvdt<<<dimGrid_matrix, dimBlock_matrix, 0, GridStream>>>(d_cable_grid);
	Myocyte_Update<<<dimGrid_myo, dimBlock_myo, 0, MyoStream>>>(dt, d_rabbit_myocyte);
	Cleft_update<<<dimGrid_cleft, dimBlock_cleft, 0, GridStream>>>(dt, d_cable_grid);
	return true;
}
template <unsigned int NMYOS, unsigned int NCable> __global__
__global__ void Compute_Myocyte_ODE(double dt, myocyte<NMYOS, NCable>* rabbit_myocyte, 
									Cell_Grid<NMYOS, NCable>* cable_grid, 
									double Istim, 
									int I_stim_location,
									double in_time){
	double Istim0;
	int location, left_cleft_location, right_cleft_location;
	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int id  = idy * NMYOS + idx;
	if (idx >= NMYOS || idy >= NCable)
		return;
	int st_start = 99 - (I_stim_location) / 2;
	int st_end   = 99 + (I_stim_location + 1) / 2;
		rabbit_myocyte -> eccODEfile(id, dt, in_time);
	for (int j = 0; j < 2; j++)
	{
		if (rabbit_myocyte -> myo_cable_location[id] == 0)
		{
			right_cleft_location = rabbit_myocyte -> myo_cable_location[id] + (NMYOS - 1) * idy;
			rabbit_myocyte -> V_myo_junc[0][id] = rabbit_myocyte -> y[39][id] - 0;
			rabbit_myocyte -> V_myo_junc[1][id] = rabbit_myocyte -> y[39][id] - cable_grid -> V_cleft[right_cleft_location];
			rabbit_myocyte -> Nao[0][id] = 140.0;
			rabbit_myocyte -> Ko[0][id] = 5.4;
			rabbit_myocyte -> Nao[1][id] = cable_grid -> Na_cleft[right_cleft_location];
			rabbit_myocyte -> Ko[1][id] = cable_grid -> K_cleft[right_cleft_location];
		}
		else if (rabbit_myocyte -> myo_cable_location[id] == (NMYOS - 1))
		{
			left_cleft_location = rabbit_myocyte -> myo_cable_location[id] - 1 + (NMYOS - 1) * idy;
			rabbit_myocyte -> V_myo_junc[0][id] = rabbit_myocyte -> y[39][id] - cable_grid -> V_cleft[left_cleft_location];
			rabbit_myocyte -> V_myo_junc[1][id] = rabbit_myocyte -> y[39][id] - 0;
			rabbit_myocyte -> Nao[0][id] = cable_grid -> Na_cleft[left_cleft_location];
			rabbit_myocyte -> Ko[0][id] = cable_grid -> K_cleft[left_cleft_location];
			rabbit_myocyte -> Nao[1][id] = 140.0;
			rabbit_myocyte -> Ko[1][id] = 5.4;
		}
		else
		{
			left_cleft_location = rabbit_myocyte -> myo_cable_location[id] - 1 + (NMYOS - 1) * idy;
			right_cleft_location = rabbit_myocyte -> myo_cable_location[id] + (NMYOS - 1) * idy;
			rabbit_myocyte -> V_myo_junc[0][id] = rabbit_myocyte -> y[39][id] - cable_grid -> V_cleft[left_cleft_location];
			rabbit_myocyte -> V_myo_junc[1][id] = rabbit_myocyte -> y[39][id] - cable_grid -> V_cleft[right_cleft_location];
			rabbit_myocyte -> Nao[0][id] = cable_grid -> Na_cleft[left_cleft_location];
			rabbit_myocyte -> Ko[0][id] = cable_grid -> K_cleft[left_cleft_location];
			rabbit_myocyte -> Nao[1][id] = cable_grid -> Na_cleft[right_cleft_location];
			rabbit_myocyte -> Ko[1][id] = cable_grid -> K_cleft[right_cleft_location];
		}
		rabbit_myocyte -> rabbit_juncODEfile(id, j, dt);
	}
	Istim0 = (idx >= st_start && idx <= st_end) ? Istim : 0.0;
	rabbit_myocyte -> rabbit_covert_to_Current(id, Istim0);
	location = rabbit_myocyte -> myo_location[id];
	cable_grid -> V_Cell[location] = rabbit_myocyte -> y[39][id];
	cable_grid -> Cell_lat_current[location] = rabbit_myocyte -> I_tot_lateral[id];
	cable_grid -> Cell_junc_current[0][location] = rabbit_myocyte -> I_tot_myo_cl[0][id];
	cable_grid -> Cell_junc_current[1][location] = rabbit_myocyte -> I_tot_myo_cl[1][id];
	cable_grid -> Cell_Na_junc_current[0][location] = rabbit_myocyte -> I_Na_tot_myo_cl[0][id];
	cable_grid -> Cell_Na_junc_current[1][location] = rabbit_myocyte -> I_Na_tot_myo_cl[1][id];
	cable_grid -> Cell_K1_junc_current[0][location] = rabbit_myocyte -> I_K1_tot_myo_cl[0][id];
	cable_grid -> Cell_K1_junc_current[1][location] = rabbit_myocyte -> I_K1_tot_myo_cl[1][id];
	cable_grid ->Cell_lat_Ina_current_foroutput[location] = rabbit_myocyte -> I_Na_lateral_foroutput[id];
	cable_grid ->Cell_lat_IK1_current_foroutput[location] = rabbit_myocyte -> I_K1_lateral_foroutput[id];
}
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Compute_Cleft_ODE(Cell_Grid<NMYOS, NCable>* cable_grid)
{
	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int id  = idy * (NMYOS - 1) + idx;
	const unsigned int id_cell  = idy * NMYOS + idx;
	if (idx >= (NMYOS - 1) || idy >= NCable)
		return;
	cable_grid -> compute_I_cleft(id);
	cable_grid -> I_gap[id] = cable_grid -> G_gap[id] * (cable_grid -> V_Cell[id_cell] - cable_grid -> V_Cell[id_cell + 1]);
}
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Arrange_Linear_Constant(Cell_Grid<NMYOS, NCable>* cable_grid)
{
	int size_of_para_matrix = NMYOS + NMYOS - 1;
	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	if (idx >= size_of_para_matrix || idy >= NCable)
		return;
	cable_grid -> arrange_constant(idx, idy);
}
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Calculate_Linear_dvdt(Cell_Grid<NMYOS, NCable>* cable_grid)
{
	int size_of_para_matrix = NMYOS + NMYOS - 1;
	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int id  = idy * size_of_para_matrix + idx;
	if (idx >= size_of_para_matrix || idy >= NCable)
		return;
	cable_grid -> solution[id] = 0;
	for (int j = 0; j < size_of_para_matrix; ++j) {
		cable_grid -> solution[id] += cable_grid -> inv_matrix[idy][idx * size_of_para_matrix + j] * cable_grid -> constant[j + idy * size_of_para_matrix];
	}
}
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Assignment_myocyte_dvdt(myocyte<NMYOS, NCable>* rabbit_myocyte, 
										Cell_Grid<NMYOS, NCable>* cable_grid)
{
	int size_of_para_matrix = NMYOS + NMYOS - 1;
	int id_of_cell_or_cleft;
	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int id  = idy * size_of_para_matrix + idx;
	if (idx >= size_of_para_matrix || idy >= NCable)
		return;
	if (idx % 2 == 0)
	{
		id_of_cell_or_cleft = cable_grid -> myFloor(idx / 2) + idy * NMYOS;
		rabbit_myocyte -> ydot[39][id_of_cell_or_cleft] = cable_grid -> solution[id];
	}
}
template <unsigned int NMYOS, unsigned int NCable>
__global__ void Assignment_Cleft_dvdt(Cell_Grid<NMYOS, NCable>* cable_grid)
{
	int id_of_cell_or_cleft;
	int size_of_para_matrix = NMYOS + NMYOS - 1;
	const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	const unsigned int id  = idy * size_of_para_matrix + idx;
	if (idx >= size_of_para_matrix || idy >= NCable)
		return;
	if (idx % 2 == 1)
	{
		id_of_cell_or_cleft = cable_grid -> myFloor(idx / 2) + idy * (NMYOS - 1);			
		cable_grid -> d_V_cleft[id_of_cell_or_cleft] = cable_grid -> solution[id];
	}
}
template <unsigned int NMYOS, unsigned int NCable>
void SpatialTissue<NMYOS, NCable>::calculate_invert_cublas(){
	const int mybatch = NCable;
	int id_cell;
	double *result_flat = (double *)malloc(mybatch * size_of_para_matrix * size_of_para_matrix * sizeof(double));
    double **results = (double **)malloc(mybatch * sizeof(double *));
	double **inputs = (double **)malloc(mybatch*sizeof(double *));
    for (int idy = 0; idy < mybatch; idy++){
      	results[idy] = result_flat + (idy * size_of_para_matrix * size_of_para_matrix);
      	inputs[idy] = result_flat + (idy * size_of_para_matrix * size_of_para_matrix);
	}	
	for (int idy = 0; idy < mybatch; idy++){
		for(int idx = 0 ;idx < size_of_para_matrix * size_of_para_matrix; idx++){
			inputs[idy][idx] = 0;
		}
	}
	for (int idy = 0; idy < mybatch; idy++){
		for (int idx = 0; idx < NMYOS + NMYOS - 1; idx++)
		{
			if(idx == 0){
				id_cell = 0 + NMYOS * idy;
				inputs[idy][idx * size_of_para_matrix + idx] = C_cell_lat[id_cell] + NUMofCell_connected[1][id_cell] * C_cell_junc[1][id_cell];
				inputs[idy][idx * size_of_para_matrix + idx + 1] = -NUMofCell_connected[1][id_cell] * C_cell_junc[1][id_cell];
			}
			else if(idx == NMYOS + NMYOS - 2){
				id_cell = NMYOS - 1 + NMYOS * idy;
				inputs[idy][idx * size_of_para_matrix + idx - 1] = 
					-NUMofCell_connected[0][id_cell] * C_cell_junc[0][id_cell];
				inputs[idy][idx * size_of_para_matrix + idx] = 
					C_cell_lat[id_cell] + NUMofCell_connected[0][id_cell] * C_cell_junc[0][id_cell]; 
			}
			else{
				if (idx % 2 == 1) 
				{
					id_cell = floor(idx / 2) + NMYOS * idy;
					inputs[idy][idx * size_of_para_matrix + idx] = -(C_cell_junc[1][id_cell] + C_cell_junc[0][id_cell + 1]);
					inputs[idy][idx * size_of_para_matrix + idx - 1] = C_cell_junc[1][id_cell];
					inputs[idy][idx * size_of_para_matrix + idx + 1] = C_cell_junc[0][id_cell + 1];
				}
				if (idx % 2 == 0) 
				{
					id_cell = floor(idx / 2) + NMYOS * idy;
					inputs[idy][idx * size_of_para_matrix + idx] = C_cell_lat[id_cell] + 
													NUMofCell_connected[0][id_cell] * C_cell_junc[0][id_cell] + 
													NUMofCell_connected[1][id_cell] * C_cell_junc[1][id_cell];
					inputs[idy][idx * size_of_para_matrix + idx - 1] = -NUMofCell_connected[0][id_cell] * C_cell_junc[0][id_cell];
					inputs[idy][idx * size_of_para_matrix + idx + 1] = -NUMofCell_connected[1][id_cell] * C_cell_junc[1][id_cell];
				}
			}
		}
	}
	FILE *fp;
	fp = fopen("Mat_Inv_out", "w");
	if (!fp) 
	{
		   fprintf(stderr, "Failed to open Mat_Inv_out.\n");
		   exit(1);
	}
	if(mybatch < 6 && size_of_para_matrix < 11)
	{
		for (int idy = 0; idy < mybatch; idy++)
		{
			if(mybatch == 1)
				fprintf(stdout, "Input Matrix, M :\n\n");
			else
				fprintf(stdout, "Input Matrix %d:\n\n", idy);
			for(int i = 0; i < size_of_para_matrix; i++)
			{
				for(int j = 0; j < size_of_para_matrix; j++)
				{	
					if(inputs[idy][i * size_of_para_matrix + j])
						fprintf(stdout,"%e ",inputs[idy][i * size_of_para_matrix + j]);
					else
						fprintf(stdout,"%e ",inputs[idy][i * size_of_para_matrix + j]);
				}
				fprintf(stdout,"\n");
			}
		}
	    fprintf(stdout,"\n\n");
	}
	else{ 
		printf("\nThe order of matrix is too large to display in terminal\n, Please open the file : Mat_Inv_out.txt located in the current folder. To see the output.\n\n");
		for (int idy = 0; idy < mybatch; idy++)
		{
			if(mybatch==1)
				fprintf(fp, "Input Matrix , M:\n\n");
			else
				fprintf(fp, "Input Matrix %d:\n\n", idy);
			for(int i = 0; i < size_of_para_matrix; i++)
			{
				for(int j = 0; j < size_of_para_matrix; j++)
			{
				if(inputs[idy][i * size_of_para_matrix + j])
					fprintf(fp, "%e ", inputs[idy][i * size_of_para_matrix + j]);
				else
					fprintf(fp, "%e ", inputs[idy][i * size_of_para_matrix + j]);
			}
				fprintf(fp, "\n");
			}
		}
		fprintf(fp, "\n\n");
	}
	invert(inputs, results, size_of_para_matrix, mybatch);
	for (int idy = 0; idy < mybatch; idy++){
		for(int i = 0; i<size_of_para_matrix; i++){
			for(int j = 0; j<size_of_para_matrix; j++){	
				cable_grid -> inv_matrix[idy][i * size_of_para_matrix + j] = results[idy][i * size_of_para_matrix + j];
			}
		}
	}
	if(mybatch < 6 && size_of_para_matrix < 11)
	{
		for (int idy = 0; idy < mybatch; idy++)
		{
			if(mybatch == 1)
				fprintf(stdout, "Inverse of the Input Matrix, Inv(M):\n\n");
			else
				fprintf(stdout, "Inverse Matrix %d:\n\n", idy);
			for(int i = 0; i < size_of_para_matrix; i++)
			{
				for(int j = 0; j < size_of_para_matrix; j++)
				{
					if(results[idy][i * size_of_para_matrix + j])
						fprintf(stdout,"%e ",results[idy][i * size_of_para_matrix + j]);
					else
						fprintf(stdout,"%e ",results[idy][i * size_of_para_matrix + j]);
				}
				fprintf(stdout,"\n");
			}
		}
		fprintf(stdout,"\n\n");
		for (int idy = 0; idy < mybatch; idy++)
		{
			if(mybatch == 1)
				fprintf(stdout, "Inverse of the Input Matrix record, Inv(M):\n\n");
			else
				fprintf(stdout, "Inverse of the Input Matrix record, Inv(M %d):\n\n", idy);
			for(int i = 0; i < size_of_para_matrix; i++)
			{
				for(int j = 0; j < size_of_para_matrix; j++)
				{
					fprintf(stdout,"%e ",cable_grid -> inv_matrix[idy][i * size_of_para_matrix + j]);
				}
				fprintf(stdout,"\n");
			}
		}
	}
	else 
	{
		for (int idy = 0; idy < mybatch; idy++)
		{
			if(mybatch == 1)
				fprintf(fp, "Inverse of the Input Matrix, Inv(M):\n\n");
			else
				fprintf(fp, "Inverse %d:\n\n", idy);
			for(int i = 0; i < size_of_para_matrix; i++)
			{
				for(int j = 0; j < size_of_para_matrix; j++)	
				{
					if(results[idy][i * size_of_para_matrix + j])
						fprintf(fp, "%e ", results[idy][i * size_of_para_matrix + j]);
					else
						fprintf(fp, "%e ", results[idy][i * size_of_para_matrix + j]);
				}
				fprintf(fp, "\n");
			}
		}
		for (int idy = 0; idy < mybatch; idy++)
		{
			if(mybatch == 1)
				fprintf(fp, "Inverse of the Input Matrix record, Inv(M):\n\n");
			else
				fprintf(fp, "Inverse of the Input Matrix record, Inv(M %d):\n\n", idy);
			for(int i = 0; i < size_of_para_matrix; i++)
			{
				for(int j = 0; j<size_of_para_matrix; j++)	
				{
					fprintf(fp, "%e ", cable_grid -> inv_matrix[idy][i * size_of_para_matrix + j]);
				}
				fprintf(fp, "\n");
			}
		}
	}
	fclose(fp);
	free(result_flat);
	free(results);
	free(inputs);
}
template <unsigned int NMYOS, unsigned int NCable>
void SpatialTissue<NMYOS, NCable>::invert(double** src, double** dst, int n, int batchSize)
{
    	cublasHandle_t handle;
    	cublascall(cublasCreate_v2(&handle));
    	int *P, *INFO;
    	cudacall(cudaMalloc(&P, n * batchSize * sizeof(int)));
    	cudacall(cudaMalloc(&INFO,  batchSize * sizeof(int)));
    	int lda = n;
    	double **A = (double **)malloc(batchSize*sizeof(double *));
    	double **A_d, *A_dflat;
    	cudacall(cudaMalloc(&A_d,batchSize*sizeof(double *)));
    	cudacall(cudaMalloc(&A_dflat, n*n*batchSize*sizeof(double)));
		A[0] = A_dflat;
    	for (int i = 1; i < batchSize; i++)
      		A[i] = A[i-1]+(n*n);
    	cudacall(cudaMemcpy(A_d,A,batchSize*sizeof(double *),cudaMemcpyHostToDevice));
 	for (int i = 0; i < batchSize; i++)
      		cudacall(cudaMemcpy(A_dflat+(i*n*n), src[i], n*n*sizeof(double), cudaMemcpyHostToDevice));
    	cublascall(cublasDgetrfBatched (handle,n,A_d,lda,P,INFO,batchSize));
    	int INFOh[batchSize];
    	cudacall(cudaMemcpy(INFOh,INFO,batchSize*sizeof(int),cudaMemcpyDeviceToHost));
    	for (int i = 0; i < batchSize; i++)
      		if(INFOh[i]  != 0)
      		{
        		fprintf(stderr, "Factorization of matrix %d Failed: Matrix may be singular\n", i);
        		cudaDeviceReset();
        		exit(EXIT_FAILURE);
      		}
    	double **C = (double **)malloc(batchSize*sizeof(double *));
    	double **C_d, *C_dflat;
    	cudacall(cudaMalloc(&C_d,batchSize*sizeof(double *)));
    	cudacall(cudaMalloc(&C_dflat, n*n*batchSize*sizeof(double)));
    	C[0] = C_dflat;
    	for (int i = 1; i < batchSize; i++)
      		C[i] = C[i-1] + (n*n);
    	cudacall(cudaMemcpy(C_d,C,batchSize*sizeof(double *),cudaMemcpyHostToDevice));
    	cublascall(cublasDgetriBatched(handle,n,(const double **)A_d,lda,P,C_d,lda,INFO,batchSize));
    	cudacall(cudaMemcpy(INFOh,INFO,batchSize*sizeof(int),cudaMemcpyDeviceToHost));
    	for (int i = 0; i < batchSize; i++)
	      	if(INFOh[i] != 0)
      		{
        		fprintf(stderr, "Inversion of matrix %d Failed: Matrix may be singular\n", i);
        		cudaDeviceReset();
      			exit(EXIT_FAILURE);
      		}
    	for (int i = 0; i < batchSize; i++)
      	cudacall(cudaMemcpy(dst[i], C_dflat + (i*n*n), n*n*sizeof(double), cudaMemcpyDeviceToHost));
	cudaFree(A_d); cudaFree(A_dflat); free(A);
	cudaFree(C_d); cudaFree(C_dflat); free(C);
    cudaFree(P); cudaFree(INFO); cublasDestroy_v2(handle);
}
template <unsigned int NMYOS, unsigned int NCable>
void SpatialTissue<NMYOS, NCable>::arrange_cells_surface_and_Capacitance(double &rate_af_am, double &rate_aj_af)
{
	A_m = 2 * pi * radius * length +
		  2 * pi * radius * radius; 
	A_myo_j = pi * radius * radius;
	int id;
	for (int idy = 0; idy < NCable; idy++){
		for (int idx = 0; idx < NMYOS; idx++)
		{
			id = idx + NMYOS * idy;
			NUMofCells_in_grid[id] = 1;
		}
	}
	for (int idy = 0; idy < NCable; idy++){
		for (int idx = 0; idx < NMYOS; idx++){
			id = idx + NMYOS * idy;
			if (idx == 0){
				NUMofCell_connected[0][id] = 0;
				NUMofCell_connected[1][id] = NUMofCells_in_grid[id + 1];
				A_cell_lat[id] = A_m - NUMofCell_connected[1][id] * A_myo_j;
				A_cell_junc[0][id] = 0;
				A_cell_junc[1][id] = A_myo_j;
			}
			else if (idx == (NMYOS - 1)){	
				NUMofCell_connected[0][id] = NUMofCells_in_grid[id - 1];
				NUMofCell_connected[1][id] = 0;
				A_cell_lat[id] = A_m - NUMofCell_connected[0][id] * A_myo_j;
				A_cell_junc[0][id] = A_myo_j;
				A_cell_junc[1][id] = 0;
			}
			else{
				A_cell_lat[id] = A_m;
				NUMofCell_connected[0][id] = NUMofCells_in_grid[id - 1];
				A_cell_lat[id] = A_cell_lat[id] - NUMofCell_connected[0][id] * A_myo_j;
				A_cell_junc[0][id] = A_myo_j;
				NUMofCell_connected[1][id] = NUMofCells_in_grid[id + 1];
				A_cell_lat[id] = A_cell_lat[id] - NUMofCell_connected[1][id] * A_myo_j;
				A_cell_junc[1][id] = A_myo_j;
			}
		}
	}
	for (int idy = 0; idy < NCable; idy++){
		C_m = 1e-8;
		for (int idx = 0; idx < NMYOS; idx++)
		{
			id = idx + NMYOS * idy;
			C_cell_lat[id] = A_cell_lat[id] * C_m;
			C_cell_junc[0][id] = A_cell_junc[0][id] * C_m;
			C_cell_junc[1][id] = A_cell_junc[1][id] * C_m;
		}
		for (int idx = 0; idx < NMYOS; idx++)
		{
			id = idx + NMYOS * idy;
			ryo_f_Ina_lat[id] = (1 - f_Ina) * (A_cell_lat[id] + A_cell_junc[0][id] * NUMofCell_connected[0][id] + A_cell_junc[1][id] * NUMofCell_connected[1][id]) / A_cell_lat[id];																 
			ryo_f_k1_lat[id] = (1 - f_k1) * (A_cell_lat[id] + A_cell_junc[0][id] * NUMofCell_connected[0][id] + A_cell_junc[1][id] * NUMofCell_connected[1][id]) / A_cell_lat[id];																	 
			ryo_f_Ina_junc[id] = f_Ina * (A_cell_lat[id] + A_cell_junc[0][id] * NUMofCell_connected[0][id] + A_cell_junc[1][id] * NUMofCell_connected[1][id]) / (A_cell_junc[0][id] * NUMofCell_connected[0][id] + A_cell_junc[1][id] * NUMofCell_connected[1][id]); 
			ryo_f_k1_junc[id] = f_k1 * (A_cell_lat[id] + A_cell_junc[0][id] * NUMofCell_connected[0][id] + A_cell_junc[1][id] * NUMofCell_connected[1][id]) / (A_cell_junc[0][id] * NUMofCell_connected[0][id] + A_cell_junc[1][id] * NUMofCell_connected[1][id]);
		}
	}
}
template <unsigned int NMYOS, unsigned int NCable>
void SpatialTissue<NMYOS, NCable>::arrange_cleft_parameter(double *g_cleft_control, double *width_of_cleft_myo, double *G_gap_myo)
{
	int id_cell, id_cleft;
	for (int idy = 0; idy < NCable; idy++){
		for (int idx = 0; idx < NMYOS - 1; idx++)
		{
			id_cell = idx + NMYOS * idy;
			id_cleft = idx + (NMYOS - 1) * idy;
			cable_grid -> G_gap[id_cleft] = G_gap_myo[idy];
			cable_grid -> width[id_cleft] = 0.001 * width_of_cleft_myo[idy];														   
			cable_grid -> vcl[id_cleft] = cable_grid -> width[id_cleft] * A_cell_junc[1][id_cell];									   
			cable_grid -> g_cleft[id_cleft] = 8 * pi * cable_grid -> vcl[id_cleft] / (A_cell_junc[1][id_cell] * 150) / g_cleft_control[idy]; 
		}
	}
	deliver_parameter_to_class();
}
template <unsigned int NMYOS, unsigned int NCable>
void SpatialTissue<NMYOS, NCable>::deliver_parameter_to_class(){
	int id, location;
	for (int idy = 0; idy < NCable; idy++){
		for (int idx = 0; idx < NMYOS; idx++)
		{
			id = idx + NMYOS * idy;
			location = rabbit_myocyte -> myo_location[id];
			V_Cell[location] = rabbit_myocyte -> y[39][id];
			rabbit_myocyte -> A_cell_lat[id] = A_cell_lat[location];
			rabbit_myocyte -> C_cell_lat[id] = C_cell_lat[location];
			rabbit_myocyte -> ryo_f_Ina_lat[id] = ryo_f_Ina_lat[location];
			rabbit_myocyte -> ryo_f_k1_lat[id] = ryo_f_k1_lat[location];
			rabbit_myocyte -> ryo_f_Ina_junc[id] = ryo_f_Ina_junc[location];
			rabbit_myocyte -> ryo_f_k1_junc[id] = ryo_f_k1_junc[location];
			for (int j = 0; j < 2; j++){
				rabbit_myocyte -> A_cell_junc[j][id] = A_cell_junc[j][location];
				rabbit_myocyte -> C_cell_junc[j][id] = C_cell_junc[j][location];
				rabbit_myocyte -> NUMofCell_connected[j][id] = NUMofCell_connected[j][location];
			}
		}
		for (int idx = 0; idx < NMYOS; idx++)
		{
			id = idx + NMYOS * idy;
			cable_grid -> V_Cell[id] = V_Cell[id];
			for (int j = 0; j < 2; j++)
			{
				cable_grid -> NUMofCell_connected[j][id] = NUMofCell_connected[j][id];
			}
		}
	}
}
template <unsigned int NMYOS, unsigned int NCable>
void SpatialTissue<NMYOS, NCable>::cpu_copy_gpu()
{
	cudaMemcpy(d_rabbit_myocyte, rabbit_myocyte, sizeof(myocyte<NMYOS, NCable>), cudaMemcpyHostToDevice);
	cudaMemcpy(d_cable_grid, cable_grid, sizeof(Cell_Grid<NMYOS, NCable>), cudaMemcpyHostToDevice);
}
