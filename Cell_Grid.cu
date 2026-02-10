template <unsigned int NMYOS, unsigned int NCable>
__host__ Cell_Grid<NMYOS, NCable>::Cell_Grid()
{
	for(int i = 0; i < (NMYOS - 1) * NCable; ++i){
		I_Na_cleft[i] = 0;
		I_K_cleft[i] = 0;
		I_tot_cleft[i] = 0;
		I_Na_tot_left_cl[i] = 0;
		I_Na_tot_right_cl[i] = 0;
		I_K_tot_left_cl[i] = 0;
		I_K_tot_right_cl[i] = 0;
	}
	for (int i = 0; i < (NMYOS - 1) * NCable; ++i)
	{
		V_cleft[i] = 0.0;         
		Na_cleft[i] = 140;        
		K_cleft[i] = 5.4;         
	}
	size_of_para_matrix = NMYOS + NMYOS - 1; 
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void Cell_Grid<NMYOS, NCable>::compute_I_cleft(int id)
{
    double Ko = 5.4;    
    double Nao = 140.0; 
    if (fabs(V_cleft[id]) < 1e-5){
        I_Na_cleft[id] = 1e-1 * g_cleft[id] / FoRT * (Na_cleft[id] - Nao * exp(-V_cleft[id] * FoRT)) / (Nao + Ko);
        I_K_cleft[id] = 1e-1 * g_cleft[id] / FoRT * (K_cleft[id] - Ko * exp(-V_cleft[id] * FoRT)) / (Nao + Ko);
    }
    else{
        I_K_cleft[id] = 1e-1 * g_cleft[id] * V_cleft[id] * (K_cleft[id] - Ko * exp(-V_cleft[id] * FoRT)) / (Nao + Ko) / (1 - exp(-V_cleft[id] * FoRT));
        I_Na_cleft[id] = 1e-1 * g_cleft[id] * V_cleft[id] * (Na_cleft[id] - Nao * exp(-V_cleft[id] * FoRT)) / (Nao + Ko) / (1 - exp(-V_cleft[id] * FoRT));
    }
    I_tot_cleft[id] = I_Na_cleft[id] + I_K_cleft[id];
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void Cell_Grid<NMYOS, NCable>::compute_dion_dt(int idx, int idy)
{
	int id_cleft = idx + (NMYOS - 1) * idy;
	int id_cell = idx + NMYOS * idy;
	I_Na_tot_left_cl[id_cleft] = Cell_Na_junc_current[1][id_cell];
	I_K_tot_left_cl[id_cleft] = Cell_K1_junc_current[1][id_cell];
	I_Na_tot_right_cl[id_cleft] = Cell_Na_junc_current[0][id_cell + 1];
	I_K_tot_right_cl[id_cleft] = Cell_K1_junc_current[0][id_cell + 1];
    d_Na_cleft[id_cleft] = 1e6 * 1 / (Frdy / 1000 * vcl[id_cleft]) * (I_Na_tot_left_cl[id_cleft] + I_Na_tot_right_cl[id_cleft] - I_Na_cleft[id_cleft]);
    d_K_cleft[id_cleft] = 1e6 * 1 / (Frdy / 1000 * vcl[id_cleft]) * (I_K_tot_left_cl[id_cleft] + I_K_tot_right_cl[id_cleft] - I_K_cleft[id_cleft]);
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void Cell_Grid<NMYOS, NCable>::cleft_update(int id, double dt)
{
    V_cleft[id] += d_V_cleft[id] * dt;
    Na_cleft[id] += d_Na_cleft[id] * dt;
    K_cleft[id] += d_K_cleft[id] * dt;
    if (Na_cleft[id] < 1e-5){
        Na_cleft[id] = 1e-5;
    }
    if (K_cleft[id] < 1e-5){
        K_cleft[id] = 1e-5;
    }
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void Cell_Grid<NMYOS, NCable>::arrange_constant(int idx, int idy)
{
	int id = idx + size_of_para_matrix * idy;
	int id_cell, id_cleft;
	if(idx == 0){
		id_cell = 0 + NMYOS * idy;
		id_cleft = 0 + (NMYOS - 1) * idy;
		constant[id] = -(Cell_lat_current[id_cell] + NUMofCell_connected[1][id_cell] * I_gap[id_cleft] + NUMofCell_connected[1][id_cell] * Cell_junc_current[1][id_cell]);
	}
	else if(idx == NMYOS + NMYOS - 2){
		id_cell = NMYOS - 1 + NMYOS * idy;
		id_cleft = NMYOS - 2 + (NMYOS - 1) * idy;
		constant[id] = -(Cell_lat_current[id_cell] - NUMofCell_connected[0][id_cell] * I_gap[id_cleft] + NUMofCell_connected[0][id_cell] * Cell_junc_current[0][id_cell]);
	}
	else{
		if (idx % 2 == 1){ 
			id_cell = myFloor(idx / 2) + NMYOS * idy;
			id_cleft = myFloor(idx / 2) + (NMYOS - 1) * idy;
			constant[id] = -(Cell_junc_current[1][id_cell] + Cell_junc_current[0][id_cell + 1] - I_tot_cleft[id_cleft]);
		}
		if (idx % 2 == 0){ 
			id_cell = myFloor(idx / 2) + NMYOS * idy;
			id_cleft = myFloor(idx / 2) + (NMYOS - 1) * idy;
			constant[id] = -(Cell_lat_current[id_cell] -
							  NUMofCell_connected[0][id_cell] * I_gap[id_cleft - 1] +
							  NUMofCell_connected[1][id_cell] * I_gap[id_cleft] +
							  NUMofCell_connected[0][id_cell] * Cell_junc_current[0][id_cell] +
							  NUMofCell_connected[1][id_cell] * Cell_junc_current[1][id_cell]);
		}
	}
}
