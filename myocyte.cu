template <unsigned int NMYOS, unsigned int NCable>
__host__ myocyte<NMYOS, NCable>::myocyte()
{
	for (int i = 0; i < NMYOS * NCable; i++)
	{
		for (int k = 1; k < 181; k++)
		{
			ydot[k][i] = 0.0;
		}
		y[1][i] = 1.9346478e-3;
		y[2][i] = 9.8083497e-1;
		y[3][i] = 9.8725639e-1;
		y[4][i] = 7.0183977e-6;
		y[5][i] = 1.0006767e+0;
		y[6][i] = 2.7168607e-2;
		y[7][i] = 1.6167249e-2;
		y[8][i] = 2.0177245e-3;
		y[9][i] = 9.9060592e-1;
		y[10][i] = 2.0177245e-3;
		y[11][i] = 9.9919974e-1;
		y[12][i] = 1.1134534e-2;
		y[13][i] = 7.3659286e-3;
		y[14][i] = 7.0206164e-1;
		y[15][i] = 3.8152792e-6;
		y[16][i] = 1.6190847e-6;
		y[17][i] = 3.9936702e+0;
		y[18][i] = 8.7151906e-1;
		y[19][i] = 9.2550649e-3;
		y[20][i] = 1.1853455e-1;
		y[21][i] = 1.0068351e-2;
		y[22][i] = 2.5328101e-4;
		y[23][i] = 1.9806253e-3;
		y[24][i] = 1.3750194e-1;
		y[25][i] = 2.2053881e-3;
		y[26][i] = 2.1025492e-2;
		y[27][i] = 1.3050993e-2;
		y[28][i] = 1.2528448e-1;
		y[29][i] = 1.3885201e-1;
		y[30][i] = 1.1601603e+0;
		y[31][i] = 4.8972938e-1;
		y[32][i] = 1.1186322e+1;
		y[33][i] = 1.1186345e+1;
		y[34][i] = 1.1186354e+1;
		y[35][i] = 1.3499804e+2;
		y[36][i] = 5.0982368e-4;
		y[37][i] = 1.4117375e-4;
		y[38][i] = 8.8621161e-5;
		y[39][i] = -8.3644872e+1;
		y[40][i] = 9.4593705e-1;
		y[41][i] = 0.0000000e+0;
		y[42][i] = 0.0000000e+0;
		y[43][i] = 5.6569568e+4;
		y[44][i] = -3.3956882e+4;
		y[45][i] = -3.1108399e+5;
		y[46][i] = 2.8847655e+5;
		y[47][i] = 2.2242769e-1;
		y[48][i] = 1.0505728e-1;
		y[49][i] = 1.9190412e-3;
		y[50][i] = 4.1517773e-5;
		y[51][i] = 3.0335421e-1;
		y[52][i] = 5.6564909e-1;
		y[53][i] = 1.0133581e-2;
		y[54][i] = 7.0125806e-5;
		y[55][i] = 8.6199183e-8;
		y[56][i] = 1.4727892e-4;
		y[57][i] = 2.6384616e-6;
		y[58][i] = 1.8246494e-8;
		y[59][i] = 2.1700214e-11;
		y[60][i] = 9.3962162e-1;
		y[61][i] = 2.7035903e-5;
		y[62][i] = 8.1950862e-5;
		y[63][i] = 5.9995858e-4;
		y[64][i] = 4.9832842e-5;
		y[65][i] = 5.9616894e-2;
		y[66][i] = 9.3960602e-1;
		y[67][i] = 2.7036121e-5;
		y[68][i] = 8.1579440e-5;
		y[69][i] = 5.9731061e-4;
		y[70][i] = 4.9832657e-5;
		y[71][i] = 5.9616591e-2;
		y[72][i] = 9.4021367e-1;
		y[73][i] = 2.7051754e-5;
		y[74][i] = 3.6434118e-6;
		y[75][i] = 2.6377746e-5;
		y[76][i] = 4.9883211e-5;
		y[77][i] = 5.9676665e-2;
		y[78][i] = 9.4019428e-1;
		y[79][i] = 2.7051371e-5;
		y[80][i] = 3.9187813e-6;
		y[81][i] = 2.8353058e-5;
		y[82][i] = 4.9881691e-5;
		y[83][i] = 5.9674869e-2;
		y[84][i] = 7.3661132e-3;
		y[85][i] = 9.9091770e-1;
		y[86][i] = 7.3660239e-3;
		y[87][i] = 9.9491063e-1;
		y[88][i] = 2.0822669e-5;
		y[89][i] = 2.9415361e-7;
		y[90][i] = 6.9281230e-7;
		y[91][i] = 1.1799709e-6;
		y[92][i] = 1.0293854e+0;
		y[93][i] = 1.0352859e+0;
		y[94][i] = 3.8745665e+2;
		y[95][i] = 1.2268406e+1;
		y[96][i] = 7.8718090e-3;
		y[97][i] = 0.0000000e+0;
		y[98][i] = 0.0000000e+0;
		y[99][i] = 0.0000000e+0;
		y[100][i] = 6.6239556e-1;
		y[101][i] = 6.7745847e-2;
		y[102][i] = 2.2767130e-5;
		y[103][i] = 8.3401763e-9;
		y[104][i] = 3.4134683e-9;
		y[105][i] = 7.6777373e-5;
		y[106][i] = 2.5371197e-3;
		y[107][i] = 1.2993396e-2;
		y[108][i] = 3.6018771e+0;
		y[109][i] = 4.1783790e-2;
		y[110][i] = 6.5562737e-5;
		y[111][i] = 8.0301027e-9;
		y[112][i] = 2.4854171e+0;
		y[113][i] = 1.1814963e+1;
		y[114][i] = 4.0971457e-4;
		y[115][i] = 1.3172872e-5;
		y[116][i] = 5.5676725e-6;
		y[117][i] = 4.8679764e-8;
		y[118][i] = 5.3905510e-13;
		y[119][i] = 3.2358638e-9;
		y[120][i] = 5.3167136e-4;
		y[121][i] = 1.9031398e-6;
		y[122][i] = 4.8236139e-6;
		y[123][i] = 1.3774572e-3;
		y[124][i] = 4.1484253e-2;
		y[125][i] = 3.6452884e-5;
		y[126][i] = 5.2490766e-10;
		y[127][i] = 3.9357735e+0;
		y[128][i] = 1.3585880e+0;
		y[129][i] = 1.8561838e-5;
		y[130][i] = 7.2917272e-6;
		y[131][i] = 4.9658763e-8;
		y[132][i] = 8.0109619e-13;
		y[133][i] = 8.5618307e-18;
		y[134][i] = 5.3346515e-14;
		y[135][i] = 1.4081682e-4;
		y[136][i] = 5.0014716e-7;
		y[137][i] = 2.0996624e-7;
		y[138][i] = 2.0971104e-6;
		y[139][i] = 1.6454393e+1;
		y[140][i] = 1.7105696e+1;
		y[141][i] = 2.9735723e+2;
		y[142][i] = 7.7288182e+1;
		y[143][i] = 6.2280569e-1;
		y[144][i] = 7.9079682e-6;
		y[145][i] = -3.3652339e-35;
		y[146][i] = -2.1845103e-33;
		y[147][i] = 4.8024867e-4;
		y[148][i] = 4.4303902e-5;
		y[149][i] = 6.4822727e-4;
		y[150][i] = 9.6049733e-3;
		y[151][i] = 6.2100610e-4;
		y[152][i] = 1.0225979e-2;
		y[153][i] = 1.4157911e-3;
		y[154][i] = 2.2218056e-3;
		y[155][i] = 1.0225024e+0;
		y[156][i] = 8.0448326e-1;
		y[157][i] = 1.4168690e-1;
		y[158][i] = 4.4775461e-3;
		y[159][i] = 2.2904345e-1;
		y[160][i] = 8.5526409e-2;
		y[161][i] = 1.4351693e-1;
		y[162][i] = 5.1034659e-2;
		y[163][i] = 8.9883074e-3;
		y[164][i] = 2.8404573e-4;
		y[165][i] = 5.7688798e-2;
		y[166][i] = 2.1541445e-2;
		y[167][i] = 3.6147457e-2;
		y[168][i] = 7.2711569e-2;
		y[169][i] = 7.2800536e-2;
		y[170][i] = 8.4540392e+0;
		y[171][i] = 5.6034188e+0;
		y[172][i] = 5.4894151e-3;
		y[173][i] = 6.2664343e-3;
		y[174][i] = 2.7577265e-2;
		y[175][i] = 4.3888456e+0;
		y[176][i] = 1.5319292e-3;
		y[177][i] = 1.5319292e-3;
		y[178][i] = 1.8361112e-3;
		y[179][i] = 4.0592685e-3;
		y[180][i] = 1.0940869e-2;
		O1_ChR2_myo[i] = 0.0;
		O2_ChR2_myo[i] = 0.0;
		C1_ChR2_myo[i] = 1.0;
		C2_ChR2_myo[i] = 0.0;
		p_ChR2_myo[i] = 0.0;
		for (int j = 0; j < 2; j++)
		{
			m_myo_junc[j][i] = 0.00193464775226056; 
			h_myo_junc[j][i] = 0.980834968871309;	
			j_myo_junc[j][i] = 0.987256394518120;	
		}
		I_tot_lateral[i] = 0.0;
		I_Na_lateral_foroutput[i] = 0.0;
		I_K1_lateral_foroutput[i] = 0.0;
		for (int j = 0; j < 2; j++)
		{
			I_tot_myo_cl[j][i] = 0.0; 
			I_Na_tot_myo_cl[j][i] = 0.0;	
			I_K1_tot_myo_cl[j][i] = 0.0;	
		}
	}
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::eccODEfile(int id, double dt, double &Ilight_myo0,double in_time)
{
	double IKs_PKAp = y[87 + 6 + 45 + 6 + 34][id] / IKstotBA;
	double ICFTR_PKAp = y[87 + 6 + 45 + 6 + 35][id] / ICFTRtotBA;
	double IKur_PKAp = y[87 + 6 + 45 + 6 + 36][id] / IKurtotBA;
	double PLM_PKAp = y[87 + 6 + 45 + 6 + 27][id] / PLMtotBA;
	int ICa_MarkovFlag = 1;
	double Acell = 20.0e3; 
	double Cmem = 1.3810e-10; 
	double Fjunc = 17.0 / (17.0 + 31.0) * 7.0 / 17.0 + 31.0 / (17.0 + 31.0) * 2.0 / 31.0;
	double Fsl = 1.0 - Fjunc; 
	double Fjunc_nak = 1.6 * 17.0 / (1.6 * 17.0 + 31.0) * 7.0 / 17.0 + 31.0 / (1.6 * 17.0 + 31.0) * 2.0 / 31.0;
	double Fsl_nak = 1.0 - Fjunc_nak; 
	double Fjunc_ncx = Fjunc;
	double Fsl_ncx = 1.0 - Fjunc_ncx;
	double Fjunc_CaL = 0.9;
	double Fsl_CaL = 1.0 - Fjunc_CaL;
	double cellLength = 100.0;		  
	double cellRadius = 11;		  
	double junctionLength = 15.0e-3;  
	double junctionRadius = 160.0e-3; 
	double distSLcyto = 0.45;		  
	double distJuncSL = 0.5; 
	double DcaJuncSL = 1.64e-6;									   
	double DcaSLcyto = 1.22e-6;									   
	double DnaJuncSL = 1.09e-5;									   
	double DnaSLcyto = 1.79e-5;									   
	double Vcell = pi * pow(cellRadius, 2) * cellLength * 1.0e-15; 
	double Vmyo = 0.65 * Vcell;
	double Vsr = 0.035 * Vcell;
	double Vsl = 0.02 * Vcell;
	double Vjunc = 0.0539 * 0.01 * Vcell;
	double SAsl = Fsl * Acell;											
	double Njunc = (Fjunc * Acell) / (pi * pow(junctionRadius, 2));		
	double SAjunc = Njunc * pi * 2.0 * junctionLength * junctionRadius; 
	double J_ca_juncsl = 1 / 1.2134e12; 
	double J_ca_slmyo = 1 / 2.68510e11; 
	double J_na_juncsl = 1 / (1.6382e12 / 3 * 100); 
	double J_na_slmyo = 1 / (1.8308e10 / 3 * 100);  
	(void)ICFTR_PKAp;
	(void)IKur_PKAp;
	(void)distSLcyto;
	(void)distJuncSL;
	(void)DcaJuncSL;
	(void)DcaSLcyto;
	(void)DnaJuncSL;
	(void)DnaSLcyto;
	(void)SAsl;
	(void)SAjunc;
	(void)J_na_juncsl;
	(void)J_na_slmyo;
	double Cli = 15.0;	
	double Clo = 150.0; 
	double Ko = 5.4;	
	double Nao = 140.0; 
	double Cao = 1.8;	
	double Mgi = 1.0;	
	double ena_junc = (1.0 / FoRT) * log(Nao / y[32][id]);	   
	double ena_sl = (1.0 / FoRT) * log(Nao / y[33][id]);	   
	double ek = (1.0 / FoRT) * log(Ko / y[35][id]);			   
	double eca_junc = (1.0 / FoRT / 2) * log(Cao / y[36][id]); 
	double eca_sl = (1.0 / FoRT / 2) * log(Cao / y[37][id]);   
	double ecl = (1.0 / FoRT) * log(Cli / Clo);				   
	double Kcoeff = 1.0; 
	double Bmax_Csqn = 2.7;	 
	double koff_csqn = 65.0; 
	double kon_csqn = 100.0; 
	i_Na_lat(id, dt, ena_junc, ena_sl, Fjunc, Fsl);
	double I_Nal = i_NaL_lat(id, dt, ena_junc, ena_sl, Fjunc, Fsl);
	double I_nabk = i_Nabk_lat(id, ena_junc, ena_sl, Fjunc, Fsl);
	double I_nak = i_NaK_lat(id, Nao, Ko, Fjunc_nak, Fsl_nak, PLM_PKAp);
	i_K1_lat(id, Ko, ek, Kcoeff);
	double I_kur = 0;
	double I_ss = 0;
	double I_kr = i_Kr_lat(id, dt, Ko, ek, Kcoeff);
	double pNaK = 0.01833;
	double fracIKspo = 0.07344;  
	double fracIKsavail = (0.2*(IKs_PKAp / fracIKspo) + 0.8);
	double Xs05 = 1.5*(2.0 - IKs_PKAp / fracIKspo);
	double pcaks_junc = -log10(y[36][id]) + 3.0;
	double pcaks_sl = -log10(y[37][id]) + 3.0;
	double gks_junc = fracIKsavail*0.07*(0.057 + 0.19 / (1 + exp((-7.2 + pcaks_junc) / 0.6))); 
	double gks_sl = fracIKsavail*0.07*(0.057 + 0.19 / (1 + exp((-7.2 + pcaks_sl) / 0.6)));     
	double eks = (1 / FoRT)*log((Ko + pNaK*Nao) / (y[35][id] + pNaK*y[34][id]));
	double xsss = 1 / (1 + exp(-(y[39][id] - Xs05) / 16.7));   
	double tauxs = 1 / (7.19e-5*(y[39][id] + 30) / (1 - exp(-0.148*(y[39][id] + 30))) + 1.31e-4*(y[39][id] + 30) / (exp(0.0687*(y[39][id] + 30)) - 1));
	ydot[13][id] = (xsss - y[13][id]) / tauxs;
	double I_ks_junc = Fjunc * gks_junc * pow(y[13][id], 2) * (y[39][id] - eks); 
	double I_ks_sl = Fsl * gks_sl * pow(y[13][id], 2) * (y[39][id] - eks); 
	double I_ks = I_ks_junc + I_ks_sl;
	double I_kp = i_Kp_lat(id, ek, Fjunc, Fsl, Kcoeff);
	double I_to = i_to_lat(id, dt, ek, Kcoeff);
	double I_ClCa, I_Clbk;
	i_Cl_lat(id, ecl, Fjunc, Fsl, I_ClCa, I_Clbk);
	double I_Ca_junc, I_Ca_sl, I_CaNa_junc, I_CaNa_sl;
	double I_Ca, I_CaNa, I_CaK;
	if (ICa_MarkovFlag == 0)
	{
		i_LCC_HH_lat(id, dt, Nao, Ko, Cao, Fjunc_CaL, Fsl_CaL, I_Ca_junc, I_Ca_sl, I_CaNa_junc, I_CaNa_sl, I_CaK);
	}
	else
	{
		i_LCC_Markov_lat(id, dt, Nao, Ko, Cao, Fjunc_CaL, Fsl_CaL, I_Ca_junc, I_Ca_sl, I_CaNa_junc, I_CaNa_sl, I_CaK);
	}
	I_Ca = I_Ca_junc + I_Ca_sl;		  
	I_CaNa = I_CaNa_junc + I_CaNa_sl; 
	ydot[43][id] = -I_Ca * Cmem / (Vmyo * 2 * Frdy) * 1e3;
	double I_ncx_junc, I_ncx_sl;
	double I_ncx = i_NCX_lat(id, Nao, Cao, Fjunc_ncx, Fsl_ncx, I_ncx_junc, I_ncx_sl);
	ydot[45][id] = 2.0 * I_ncx * Cmem / (Vmyo * 2.0 * Frdy) * 1.0e3; 
	double I_pca_junc, I_pca_sl;
	double I_pca = i_pca_lat(id, Fjunc, Fsl, I_pca_junc, I_pca_sl);
	ydot[44][id] = -I_pca * Cmem / (Vmyo * 2 * Frdy) * 1e3;
	double I_cabk_junc, I_cabk_sl;
	double I_cabk = i_Cabk_lat(id, Fjunc, Fsl, I_cabk_junc, I_cabk_sl, eca_junc, eca_sl);
	ydot[46][id] = -I_cabk * Cmem / (Vmyo * 2.0 * Frdy) * 1.0e3;
	double Icftr = 0; 
	double J_SRCarel, J_SRleak;
	i_RyR_lat(id, J_SRCarel, J_SRleak);
	double J_serca = i_SERCA_lat(id);
	double J_CaB_cytosol, J_CaB_junction, J_CaB_sl;
	Na_Ca_buffer(id, Mgi, Vmyo, Vjunc, Vsl, J_CaB_cytosol, J_CaB_junction, J_CaB_sl);
	ydot[30][id] = kon_csqn * y[31][id] * (Bmax_Csqn - y[30][id]) - koff_csqn * y[30][id];	  
	ydot[31][id] = J_serca * Vmyo / Vsr - (J_SRleak * Vmyo / Vsr + J_SRCarel ) - ydot[30][id]; 
	ydot[32][id] = 0.0;
	ydot[33][id] = 0.0;
	ydot[34][id] = 0.0; 
	ydot[35][id] = 0.0;
	double I_Ca_tot_junc = I_Ca_junc + I_cabk_junc + I_pca_junc - 2.0 * I_ncx_junc;																									  
	double I_Ca_tot_sl = I_Ca_sl + I_cabk_sl + I_pca_sl - 2.0 * I_ncx_sl;																											  
	ydot[36][id] = -I_Ca_tot_junc * Cmem / (Vjunc * 2.0 * Frdy) + J_ca_juncsl / Vjunc * (y[37][id] - y[36][id]) - J_CaB_junction + (J_SRCarel)*Vsr / Vjunc + J_SRleak * Vmyo / Vjunc; 
	ydot[37][id] = -I_Ca_tot_sl * Cmem / (Vsl * 2.0 * Frdy) + J_ca_juncsl / Vsl * (y[36][id] - y[37][id]) + J_ca_slmyo / Vsl * (y[38][id] - y[37][id]) - J_CaB_sl;					  
	ydot[38][id] = -J_serca - J_CaB_cytosol + J_ca_slmyo / Vmyo * (y[37][id] - y[38][id]);																							  
	double junc_sl = J_ca_juncsl / Vsl * (y[36][id] - y[37][id]);
	double sl_junc = J_ca_juncsl / Vjunc * (y[37][id] - y[36][id]);
	double sl_myo = J_ca_slmyo / Vsl * (y[38][id] - y[37][id]);
	double myo_sl = J_ca_slmyo / Vmyo * (y[37][id] - y[38][id]);
	(void)junc_sl;
	(void)sl_junc;
	(void)sl_myo;
	(void)myo_sl;
	double i_Ca_lateral = I_Ca + I_CaNa + I_cabk + I_pca + I_ncx;
	double i_K_lateral = I_to + I_kr + I_ks + I_CaK + I_kp + I_kur + I_ss;
	double i_Cl_lateral = I_ClCa + I_Clbk + Icftr;
	double i_Na_lateral = I_nabk + I_Nal + I_nak;
	i_comp_lateral[id] = i_Na_lateral +
						 i_K_lateral +
						 i_Ca_lateral +
						 i_Cl_lateral ;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::DAD_eccODEfile(int id, double dt, double &Ilight_myo0,double in_time)
{
	double IKs_PKAp = y[87 + 6 + 45 + 6 + 34][id] / IKstotBA;
	double ICFTR_PKAp = y[87 + 6 + 45 + 6 + 35][id] / ICFTRtotBA;
	double IKur_PKAp = y[87 + 6 + 45 + 6 + 36][id] / IKurtotBA;
	double PLM_PKAp = y[87 + 6 + 45 + 6 + 27][id] / PLMtotBA;
	int ICa_MarkovFlag = 1;
	double Acell = 20.0e3; 
	double Cmem = Acell * 1.0e-14; 
	double Fjunc =0.11;
	double Fsl = 1.0 - Fjunc; 
	double Fjunc_nak = 1.6 * 17.0 / (1.6 * 17.0 + 31.0) * 7.0 / 17.0 + 31.0 / (1.6 * 17.0 + 31.0) * 2.0 / 31.0;
	double Fsl_nak = 1.0 - Fjunc_nak; 
	double Fjunc_ncx = Fjunc;
	double Fsl_ncx = 1.0 - Fjunc_ncx;
	double Fjunc_CaL = 0.9;
	double Fsl_CaL = 1.0 - Fjunc_CaL;
	double cellLength = 100.0;		  
	double cellRadius = 10.25;		  
	double junctionLength = 15.0e-3;  
	double junctionRadius = 160.0e-3; 
	double distSLcyto = 0.45;		  
	double distJuncSL = 0.3;									   
	double DcaJuncSL = 1.64e-6;									   
	double DcaSLcyto = 1.22e-6;									   
	double DnaJuncSL = 1.09e-5;									   
	double DnaSLcyto = 1.79e-5;									   
	double Vcell = pi * pow(cellRadius, 2) * cellLength * 1.0e-15; 
	double Vmyo = 0.65 * Vcell;
	double Vsr = 0.035 * Vcell;
	double Vsl = 0.02 * Vcell;
	double Vjunc = 0.0539 * 0.01 * Vcell;
	double SAsl = Fsl * Acell;											
	double Njunc = (Fjunc * Acell) / (pi * pow(junctionRadius, 2));		
	double SAjunc = Njunc * pi * 2.0 * junctionLength * junctionRadius; 
	double J_ca_juncsl = 1 / 1.2134e12; 
	double J_ca_slmyo = 1 / 2.68510e11; 
	double J_na_juncsl = 1 / (1.6382e12 / 3 * 100); 
	double J_na_slmyo = 1 / (1.8308e10 / 3 * 100);  
	(void)ICFTR_PKAp;
	(void)IKur_PKAp;
	(void)distSLcyto;
	(void)distJuncSL;
	(void)DcaJuncSL;
	(void)DcaSLcyto;
	(void)DnaJuncSL;
	(void)DnaSLcyto;
	(void)SAsl;
	(void)SAjunc;
	(void)J_na_juncsl;
	(void)J_na_slmyo;
	double Cli = 15.0;	
	double Clo = 150.0; 
	double Ko = 5.4;	
	double Nao = 140.0; 
	double Cao = 1.8;	
	double Mgi = 1.0;	
	double ena_junc = (1.0 / FoRT) * log(Nao / y[32][id]);	   
	double ena_sl = (1.0 / FoRT) * log(Nao / y[33][id]);	   
	double ek = (1.0 / FoRT) * log(Ko / y[35][id]);			   
	double eca_junc = (1.0 / FoRT / 2) * log(Cao / y[36][id]); 
	double eca_sl = (1.0 / FoRT / 2) * log(Cao / y[37][id]);   
	double ecl = (1.0 / FoRT) * log(Cli / Clo);				   
	double Kcoeff = 1.0; 
	double Bmax_Csqn = 2.7;	 
	double koff_csqn = 65.0; 
	double kon_csqn = 100.0; 
	i_Na_lat(id, dt, ena_junc, ena_sl, Fjunc, Fsl);
	double I_Nal = i_NaL_lat(id, dt, ena_junc, ena_sl, Fjunc, Fsl);
	double I_nabk = i_Nabk_lat(id, ena_junc, ena_sl, Fjunc, Fsl);
	double I_nak = i_NaK_lat(id, Nao, Ko, Fjunc_nak, Fsl_nak, PLM_PKAp);
	i_K1_lat(id, Ko, ek, Kcoeff);
	double I_kur = 0;
	double I_ss = 0;
	double I_kr = i_Kr_lat(id, dt, Ko, ek, Kcoeff);
	double pNaK = 0.01833;
	double fracIKspo = 0.07344;  
	double fracIKsavail = (0.2*(IKs_PKAp / fracIKspo) + 0.8);
	double Xs05 = 1.5*(2.0 - IKs_PKAp / fracIKspo);
	double pcaks_junc = -log10(y[36][id]) + 3.0;
	double pcaks_sl = -log10(y[37][id]) + 3.0;
	double gks_junc = fracIKsavail*0.07*(0.057 + 0.19 / (1 + exp((-7.2 + pcaks_junc) / 0.6))); 
	double gks_sl = fracIKsavail*0.07*(0.057 + 0.19 / (1 + exp((-7.2 + pcaks_sl) / 0.6)));     
	double eks = (1 / FoRT)*log((Ko + pNaK*Nao) / (y[35][id] + pNaK*y[34][id]));
	double xsss = 1 / (1 + exp(-(y[39][id] - Xs05) / 16.7));   
	double tauxs = 1 / (7.19e-5*(y[39][id] + 30) / (1 - exp(-0.148*(y[39][id] + 30))) + 1.31e-4*(y[39][id] + 30) / (exp(0.0687*(y[39][id] + 30)) - 1));
	ydot[13][id] = (xsss - y[13][id]) / tauxs;
	double I_ks_junc = Fjunc * gks_junc * pow(y[13][id], 2) * (y[39][id] - eks); 
	double I_ks_sl = Fsl * gks_sl * pow(y[13][id], 2) * (y[39][id] - eks); 
	double I_ks = I_ks_junc + I_ks_sl;
	double I_kp = i_Kp_lat(id, ek, Fjunc, Fsl, Kcoeff);
	double I_to = i_to_lat(id, dt, ek, Kcoeff);
	double I_ClCa, I_Clbk;
	i_Cl_lat(id, ecl, Fjunc, Fsl, I_ClCa, I_Clbk);
	double I_Ca_junc, I_Ca_sl, I_CaNa_junc, I_CaNa_sl;
	double I_Ca, I_CaNa, I_CaK;
	if (ICa_MarkovFlag == 0)
	{
		i_LCC_HH_lat(id, dt, Nao, Ko, Cao, Fjunc_CaL, Fsl_CaL, I_Ca_junc, I_Ca_sl, I_CaNa_junc, I_CaNa_sl, I_CaK);
	}
	else
	{
		i_LCC_Markov_lat(id, dt, Nao, Ko, Cao, Fjunc_CaL, Fsl_CaL, I_Ca_junc, I_Ca_sl, I_CaNa_junc, I_CaNa_sl, I_CaK);
	}
	I_Ca = I_Ca_junc + I_Ca_sl;		  
	I_CaNa = I_CaNa_junc + I_CaNa_sl; 
	ydot[43][id] = -I_Ca * Cmem / (Vmyo * 2 * Frdy) * 1e3;
	double I_ncx_junc, I_ncx_sl;
	double I_ncx = i_NCX_lat(id, Nao, Cao, Fjunc_ncx, Fsl_ncx, I_ncx_junc, I_ncx_sl);
	ydot[45][id] = 2.0 * I_ncx * Cmem / (Vmyo * 2.0 * Frdy) * 1.0e3; 
	double I_pca_junc, I_pca_sl;
	double I_pca = i_pca_lat(id, Fjunc, Fsl, I_pca_junc, I_pca_sl);
	ydot[44][id] = -I_pca * Cmem / (Vmyo * 2 * Frdy) * 1e3;
	double I_cabk_junc, I_cabk_sl;
	double I_cabk = i_Cabk_lat(id, Fjunc, Fsl, I_cabk_junc, I_cabk_sl, eca_junc, eca_sl);
	ydot[46][id] = -I_cabk * Cmem / (Vmyo * 2.0 * Frdy) * 1.0e3;
	double Icftr = 0; 
	double J_SRCarel, J_SRleak;
	i_RyR_lat(id, J_SRCarel, J_SRleak);
	double J_serca = i_SERCA_lat(id);
	double J_CaB_cytosol, J_CaB_junction, J_CaB_sl;
	Na_Ca_buffer(id, Mgi, Vmyo, Vjunc, Vsl, J_CaB_cytosol, J_CaB_junction, J_CaB_sl);
	double G_spon = 0.6;
	double beta = Vsl/Vsr;      
	double t0 = 425.0;      
	double tau1 = 10.0;     
	double tau2 = 30.0;     
	double g1 = 1.0 / (1.0 + exp(-(in_time - t0)/tau1));
	double g2 = 1.0 / (1.0 + exp((in_time - t0)/tau2));
	double J_spon = G_spon * (beta *y[31][id] - y[36][id]) * g1 * g2; 
	ydot[30][id] = kon_csqn * y[31][id] * (Bmax_Csqn - y[30][id]) - koff_csqn * y[30][id];	  
	ydot[31][id] = J_serca * Vmyo / Vsr - (J_SRleak * Vmyo / Vsr + J_SRCarel + J_spon) - ydot[30][id]; 
	ydot[32][id] = 0.0;
	ydot[33][id] = 0.0;
	ydot[34][id] = 0.0; 
	ydot[35][id] = 0.0;
	double I_Ca_tot_junc = I_Ca_junc + I_cabk_junc + I_pca_junc - 2.0 * I_ncx_junc;																									  
	double I_Ca_tot_sl = I_Ca_sl + I_cabk_sl + I_pca_sl - 2.0 * I_ncx_sl;																											  
	ydot[36][id] = -I_Ca_tot_junc * Cmem / (Vjunc * 2.0 * Frdy) + J_ca_juncsl / Vjunc * (y[37][id] - y[36][id]) - J_CaB_junction + (J_SRCarel + J_spon)*Vsr / Vjunc + J_SRleak * Vmyo / Vjunc; 
	ydot[37][id] = -I_Ca_tot_sl * Cmem / (Vsl * 2.0 * Frdy) + J_ca_juncsl / Vsl * (y[36][id] - y[37][id]) + J_ca_slmyo / Vsl * (y[38][id] - y[37][id]) - J_CaB_sl;					  
	ydot[38][id] = -J_serca - J_CaB_cytosol + J_ca_slmyo / Vmyo * (y[37][id] - y[38][id]);																							  
	double junc_sl = J_ca_juncsl / Vsl * (y[36][id] - y[37][id]);
	double sl_junc = J_ca_juncsl / Vjunc * (y[37][id] - y[36][id]);
	double sl_myo = J_ca_slmyo / Vsl * (y[38][id] - y[37][id]);
	double myo_sl = J_ca_slmyo / Vmyo * (y[37][id] - y[38][id]);
	(void)junc_sl;
	(void)sl_junc;
	(void)sl_myo;
	(void)myo_sl;
	double i_Ca_lateral = I_Ca + I_CaNa + I_cabk + I_pca + I_ncx;
	double i_K_lateral = I_to + I_kr + I_ks + I_CaK + I_kp + I_kur + I_ss;
	double i_Cl_lateral = I_ClCa + I_Clbk + Icftr;
	double i_Na_lateral = I_nabk + I_Nal + I_nak;
	i_comp_lateral[id] = i_Na_lateral +
						 i_K_lateral +
						 i_Ca_lateral +
						 i_Cl_lateral ;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::rabbit_juncODEfile(int id, int j, double dt)
{
	double Nai = y[34][id]; 
	double Ki = y[35][id];	
	double ena = (1.0 / FoRT) * log(Nao[j][id] / Nai);
	double ek = (1.0 / FoRT) * log(Ko[j][id] / Ki);
	;
	double GNa = 12.0; 
	double Kcoeff = 1.0; 
	double inashift = 0.0;
	double alphaCKII = 0.0;
	double am;
	if (fabs(0.1 * (V_myo_junc[j][id] - inashift + 47.13)) < 1e-4)
	{
		am = 0.32 / 0.1;
	}
	else
	{
		am = 0.32 * (V_myo_junc[j][id] - inashift + 47.13) / (1.0 - exp(-0.1 * (V_myo_junc[j][id] - inashift + 47.13)));
	}
	double bm = 0.08 * exp(-(V_myo_junc[j][id] - inashift) / 11.0);
	double ah, aj, bh, bj;
	if ((V_myo_junc[j][id] - inashift) >= -40.0)
	{
		ah = 0.0;
		aj = 0.0;
		bh = 1 / (0.13*(1 + exp(-((V_myo_junc[j][id] - inashift) + 10.66) / 11.1))); 
		bj = 0.3 * exp(-2.535e-7 * (V_myo_junc[j][id] - inashift)) / (1.0 + exp(-0.1 * ((V_myo_junc[j][id] - inashift) + 32.0)));
	}
	else
	{
		ah = 0.135 * exp((80.0 + (V_myo_junc[j][id] - inashift)) / -6.8);
		bh = 3.56*exp(0.079*(V_myo_junc[j][id] - inashift)) + 3.1e5*exp(0.35*(V_myo_junc[j][id] - inashift)); 
		aj = (1.0 + alphaCKII) * ((-1.2714e5 * exp(0.2444 * (V_myo_junc[j][id] - inashift)) - 3.474e-5 * exp(-0.04391 * (V_myo_junc[j][id] - inashift))) * ((V_myo_junc[j][id] - inashift) + 37.78) / (1.0 + exp(0.311 * ((V_myo_junc[j][id] - inashift) + 79.23))));
		bj = 0.1212 * exp(-0.01052 * (V_myo_junc[j][id] - inashift)) / (1.0 + exp(-0.1378 * ((V_myo_junc[j][id] - inashift) + 40.14)));
	}
	double taum = 1.0 / (am + bm);
	double tauh = 1.0 / (ah + bh);
	double tauj = 1.0 / (aj + bj);
	d_m_myo_junc[j][id] = (am / (am + bm) - ((am / (am + bm)) - m_myo_junc[j][id]) * exp(-dt / taum) - m_myo_junc[j][id]) / dt;
	d_h_myo_junc[j][id] = (ah / (ah + bh) - ((ah / (ah + bh)) - h_myo_junc[j][id]) * exp(-dt / tauh) - h_myo_junc[j][id]) / dt;
	d_j_myo_junc[j][id] = (aj / (aj + bj) - ((aj / (aj + bj)) - j_myo_junc[j][id]) * exp(-dt / tauj) - j_myo_junc[j][id]) / dt;
	double I_Na = GNa * pow(m_myo_junc[j][id], 3) * h_myo_junc[j][id] * j_myo_junc[j][id] * (V_myo_junc[j][id] - ena);
	double aki = 1.02 / (1 + exp(0.2385 * (V_myo_junc[j][id] - ek - 59.215)));
	double bki = (0.49124 * exp(0.08032 * (V_myo_junc[j][id] + 5.476 - ek)) + exp(0.06175 * (V_myo_junc[j][id] - ek - 594.31))) / (1.0 + exp(-0.5143 * (V_myo_junc[j][id] - ek + 4.753)));
	double kiss = aki / (aki + bki);
	double I_ki;
	I_ki = 0.9 * sqrt(Ko[j][id] / 5.4) * kiss * (V_myo_junc[j][id] - ek) * Kcoeff;
	i_Na_myo_cl[j][id] = I_Na; 
	i_K1_myo_cl[j][id] = I_ki;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::rabbit_covert_to_Current(int id, double Istim)
{
	I_tot_lateral[id] = 
	(ryo_f_Ina_lat[id] * i_Na_lateral[id] +
	ryo_f_k1_lat[id] * i_K1_lateral[id] - Istim) * C_cell_lat[id] +
	i_comp_lateral[id] * (C_cell_lat[id] + NUMofCell_connected[0][id] * C_cell_junc[0][id] + NUMofCell_connected[1][id] * C_cell_junc[1][id]);
	I_Na_lateral_foroutput[id] = ryo_f_Ina_lat[id] * i_Na_lateral[id] * C_cell_lat[id];
	I_K1_lateral_foroutput[id] = ryo_f_k1_lat[id] * i_K1_lateral[id] * C_cell_lat[id];
	for (int j = 0; j < 2; j++)
	{
		I_Na_tot_myo_cl[j][id] = ryo_f_Ina_junc[id] * i_Na_myo_cl[j][id] * C_cell_junc[j][id];
		I_K1_tot_myo_cl[j][id] = ryo_f_k1_junc[id] * i_K1_myo_cl[j][id] * C_cell_junc[j][id];
		I_tot_myo_cl[j][id] = I_Na_tot_myo_cl[j][id] + I_K1_tot_myo_cl[j][id];
	}
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::rabbit_update(int id, double dt)
{
	for (int i = 1; i < 181; i++)
	{
		y[i][id] += ydot[i][id] * dt;
	}
	O1_ChR2_myo[id] += d_O1_ChR2_myo[id] * dt;
	O2_ChR2_myo[id] += d_O2_ChR2_myo[id] * dt;
	C1_ChR2_myo[id] += d_C1_ChR2_myo[id] * dt;
	C2_ChR2_myo[id] += d_C2_ChR2_myo[id] * dt;
	p_ChR2_myo[id] += d_p_ChR2_myo[id] * dt;
	for (int j = 0; j < 2; j++)
	{
		m_myo_junc[j][id] += d_m_myo_junc[j][id] * dt;
		h_myo_junc[j][id] += d_h_myo_junc[j][id] * dt;
		j_myo_junc[j][id] += d_j_myo_junc[j][id] * dt;
	}
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::i_Na_lat(int id, double dt, double &ena_junc, double &ena_sl, double &Fjunc, double &Fsl)
{
	double GNa = 12.0; 
	double inashift = 0.0;
	double alphaCKII = 0.0;
	double am;
	if (fabs(0.1 * (y[39][id] - inashift + 47.13)) < 1e-4)
	{
		am = 0.32 / 0.1;
	}
	else
	{
		am = 0.32 * (y[39][id] - inashift + 47.13) / (1.0 - exp(-0.1 * (y[39][id] - inashift + 47.13)));
	}
	double bm = 0.08 * exp(-(y[39][id] - inashift) / 11.0);
	double ah, aj, bh, bj;
	if ((y[39][id] - inashift) >= -40.0)
	{
		ah = 0.0;
		aj = 0.0;
		bh = 1 / (0.13*(1 + exp(-((y[39][id] - inashift) + 10.66) / 11.1))); 
		bj = 0.3 * exp(-2.535e-7 * (y[39][id] - inashift)) / (1.0 + exp(-0.1 * ((y[39][id] - inashift) + 32.0)));
	}
	else
	{
		ah = 0.135 * exp((80.0 + (y[39][id] - inashift)) / -6.8);
		bh = 3.56*exp(0.079*(y[39][id] - inashift)) + 3.1e5*exp(0.35*(y[39][id] - inashift)); 
		aj = (1.0 + alphaCKII) * ((-1.2714e5 * exp(0.2444 * (y[39][id] - inashift)) - 3.474e-5 * exp(-0.04391 * (y[39][id] - inashift))) * ((y[39][id] - inashift) + 37.78) / (1.0 + exp(0.311 * ((y[39][id] - inashift) + 79.23))));
		bj = 0.1212 * exp(-0.01052 * (y[39][id] - inashift)) / (1.0 + exp(-0.1378 * ((y[39][id] - inashift) + 40.14)));
	}
	double tauh = 1.0 / (ah + bh);
	double tauj = 1.0 / (aj + bj);
	double taum = 1.0 / (am + bm);
	ydot[1][id] = (am / (am + bm) - ((am / (am + bm)) - y[1][id]) * exp(-dt / taum) - y[1][id]) / dt;
	ydot[2][id] = (ah / (ah + bh) - ((ah / (ah + bh)) - y[2][id]) * exp(-dt / tauh) - y[2][id]) / dt;
	ydot[3][id] = (aj / (aj + bj) - ((aj / (aj + bj)) - y[3][id]) * exp(-dt / tauj) - y[3][id]) / dt;
	double I_Na_junc = Fjunc * GNa * pow(y[1][id], 3) * y[2][id] * y[3][id] * (y[39][id] - ena_junc);
	double I_Na_sl = Fsl * GNa * pow(y[1][id], 3) * y[2][id] * y[3][id] * (y[39][id] - ena_sl);
	i_Na_lateral[id] = I_Na_junc + I_Na_sl;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::i_K1_lat(int id, double &Ko, double &ek, double &Kcoeff)
{
	double aki = 1.02 / (1 + exp(0.2385 * (y[39][id] - ek - 59.215)));
	double bki = (0.49124 * exp(0.08032 * (y[39][id] + 5.476 - ek)) + exp(0.06175 * (y[39][id] - ek - 594.31))) / (1.0 + exp(-0.5143 * (y[39][id] - ek + 4.753)));
	double kiss = aki / (aki + bki);
	double I_ki;
	I_ki = 0.9 * sqrt(Ko / 5.4) * kiss * (y[39][id] - ek) * Kcoeff;
	i_K1_lateral[id] = I_ki;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_NaL_lat(int id, double dt, double &ena_junc, double &ena_sl, double &Fjunc, double &Fsl)
{
	double deltGbarNal_CKII = 0.0;
	double GbarNal = 0.0065 * (1.0 + deltGbarNal_CKII) * 2.0; 
	double hlss = 1.0 / (1.0 + exp((y[39][id] + 91.0) / 6.1));
	double tauhl = 600.0; 
	ydot[47][id] = (hlss - (hlss - y[47][id]) * exp(-dt / tauhl) - y[47][id]) / dt; 
	double I_Nalj = Fjunc * GbarNal * pow(y[1][id], 3) * y[47][id] * (y[39][id] - ena_junc);
	double I_Nalsl = Fsl * GbarNal * pow(y[1][id], 3) * y[47][id] * (y[39][id] - ena_sl);
	double I_Nal = I_Nalj + I_Nalsl;
	return I_Nal;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_Nabk_lat(int id, double &ena_junc, double &ena_sl, double &Fjunc, double &Fsl)
{
	double GNaB = 0.297e-3; 
	double I_nabk_junc = Fjunc * GNaB * (y[39][id] - ena_junc);
	double I_nabk_sl = Fsl * GNaB * (y[39][id] - ena_sl);
	double I_nabk = I_nabk_junc + I_nabk_sl;
	return I_nabk;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_NaK_lat(int id, double &Nao, double &Ko, double &Fjunc_nak, double &Fsl_nak, double &PLM_PKAp)
{
	double IbarNaK = 1.90719; 
	double KmNaip = 11.0; 
	double KmKo = 1.5;	  
	double sigma = (exp(Nao / 67.3) - 1.0) / 7.0;
	double fnak = 1.0 / (1.0 + 0.1245 * exp(-0.1 * y[39][id] * FoRT) + 0.0365 * sigma * exp(-y[39][id] * FoRT));
	double fracPKA_PLMo = 0.116738;	  
	double fracPKA_PLMiso = 0.859251; 
	double kPKA_PLM = KmNaip * (1.0 - 0.7019) / (fracPKA_PLMiso / fracPKA_PLMo - 1.0); 
	double KmNaip_PKA = -kPKA_PLM + kPKA_PLM * (PLM_PKAp / fracPKA_PLMo);
	KmNaip = KmNaip - KmNaip_PKA;
	double I_nak_junc = Fjunc_nak * IbarNaK * fnak * Ko / (1.0 + pow((KmNaip / y[32][id]), 4)) / (Ko + KmKo);
	double I_nak_sl = Fsl_nak * IbarNaK * fnak * Ko / (1.0 + pow((KmNaip / y[32][id]), 4)) / (Ko + KmKo);
	double I_nak = I_nak_junc + I_nak_sl;
	return I_nak;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_Kur_lat(int id, double dt, double &ek, double &Kcoeff, double &IKur_PKAp)
{
	double Gkur1 = 1.1 * 0.16; 
	double Gkur2 = 0.14;	   
	double xurss = 1.0 / (1.0 + exp(-(y[39][id] + 15.0) / 14.0));
	double yurss = 1.0 / (1.0 + exp((y[39][id] + 48.0) / 6.2));
	double tauxur = 0.95 + 0.05 * exp(-0.08 * y[39][id]); 
	ydot[84][id] = (xurss - (xurss - y[84][id]) * exp(-dt / tauxur) - y[84][id]) / dt;
	double tauxur2 = 1.0 + 7.0 / (1.0 + exp(-(y[39][id] + 45.0) / 8.0)) + 20.0 * exp(-((y[39][id] + 35.0) / 10.0) * ((y[39][id] + 35.0) / 10.0));
	double tauyur1 = 400.0 + 900.0 * exp(-((y[39][id] + 55.0) / 16.0) * ((y[39][id] + 55.0) / 16.0)) - 250.0 / (1.0 + exp(-(y[39][id] + 60.0) / 8.0)); 
	ydot[85][id] = (yurss - (yurss - y[85][id]) * exp(-dt / tauyur1) - y[85][id]) / dt;
	double tauyur2 = 400.0 + 900.0 * exp(-((y[39][id] + 55.0) / 16.0) * ((y[39][id] + 55.0) / 16.0)) + 550.0 / (1.0 + exp(-(y[39][id] + 60.0) / 8.0)); 
	ydot[87][id] = (yurss - (yurss - y[87][id]) * exp(-dt / tauyur2) - y[87][id]) / dt;
	double fracIKurp0 = 0.437635;	
	double fracIKurpISO = 0.718207; 
	double a_Kur = (1.20 - 1.0) / (fracIKurpISO / fracIKurp0 - 1.0);
	double fracIKuravail = (1.0 - a_Kur) + a_Kur * (IKur_PKAp / fracIKurp0); 
	double I_kur1 = Kcoeff * fracIKuravail * Gkur1 * y[84][id] * y[85][id] * (y[39][id] - ek); 
	double I_kur2 = Kcoeff * Gkur2 * y[84][id] * y[87][id] * (y[39][id] - ek); 
	double I_kur = I_kur1 + I_kur2;
	return I_kur;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_ss_lat(int id, double dt, double &ek, double &Kcoeff)
{
	double Gss = 0.15; 
	double xurss = 1.0 / (1.0 + exp(-(y[39][id] + 15.0) / 14.0));
	double xssss = xurss; 
	double tauxss = 70.0 * exp(-((y[39][id] + 43.0) / 30.0) * ((y[39][id] + 43.0) / 30.0)) + 14.0; 
	ydot[86][id] = (xssss - (xssss - y[86][id]) * exp(-dt / tauxss) - y[86][id]) / dt;
	double I_ss = Kcoeff * Gss * y[86][id] * (y[39][id] - ek); 
	return I_ss;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_Kr_lat(int id, double dt, double &Ko, double &ek, double &Kcoeff)
{
	double gkr = 0.03 * sqrt(Ko / 5.4);
	double xrss = 1.0 / (1.0 + exp(-(y[39][id] + 50.0) / 7.5));
	double tauxr;
	if (fabs(y[39][id] + 7.0) < 1e-4)
	{
		tauxr = 1.0 / (1.38e-3 / 0.123 + 6.1e-4 * (y[39][id] + 10.0) / (exp(0.145 * (y[39][id] + 10.0)) - 1.0));
	}
	else if (fabs(y[39][id] + 10.0) < 1e-4)
	{
		tauxr = 1.0 / (1.38e-3 * (y[39][id] + 7.0) / (1.0 - exp(-0.123 * (y[39][id] + 7.0))) + 6.1e-4 / 0.145);
	}
	else
	{
		tauxr = 1.0 / (1.38e-3 * (y[39][id] + 7.0) / (1.0 - exp(-0.123 * (y[39][id] + 7.0))) + 6.1e-4 * (y[39][id] + 10.0) / (exp(0.145 * (y[39][id] + 10.0)) - 1.0));
	}
	ydot[12][id] = (xrss - (xrss - y[12][id]) * exp(-dt / tauxr) - y[12][id]) / dt;
	double rkr = 1.0 / (1.0 + exp((y[39][id] + 33.0) / 22.4));
	double I_kr = Kcoeff * gkr * y[12][id] * rkr * (y[39][id] - ek);
	return I_kr;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_Kp_lat(int id, double &ek, double &Fjunc, double &Fsl, double &Kcoeff)
{
	double gkp = 0.001;
	double kp_kp = 1.0 / (1.0 + exp(7.488 - y[39][id] / 5.98));
	double I_kp_junc = Kcoeff * Fjunc * gkp * kp_kp * (y[39][id] - ek);
	double I_kp_sl = Kcoeff * Fsl * gkp * kp_kp * (y[39][id] - ek);
	double I_kp = I_kp_junc + I_kp_sl;
	return I_kp;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_to_lat(int id, double dt, double &ek, double &Kcoeff)
{
	double GtoSlow = 0.06;  
	double GtoFast = 0.02; 
	double xtoss = 1.0 / (1.0 + exp(-(y[39][id] + 3.0) / 15.0));
	double ytoss = 1.0 / (1.0 + exp((y[39][id] + 33.5) / 10.0));
	double rtoss = 1.0 / (1.0 + exp((y[39][id] + 33.5) / 10.0)); 
	double tauxtos = 9.0 / (1.0 + exp(-(y[39][id] + 3.0) / 15.0)) + 0.5;
	double tauytos = 3000.0 / (1.0 + exp((y[39][id] + 60.0) / 10.0)) + 30.0;
	double taurtos = tauytos; 
	ydot[8][id] = (xtoss - (xtoss - y[8][id]) * exp(-dt / tauxtos) - y[8][id]) / dt;
	ydot[9][id] = (ytoss - (ytoss - y[9][id]) * exp(-dt / tauytos) - y[9][id]) / dt;
	ydot[40][id] = (rtoss - (rtoss - y[40][id]) * exp(-dt / taurtos) - y[40][id]) / dt;
	double I_tos = GtoSlow * y[8][id] * (y[9][id] + 0.5 * y[40][id]) * (y[39][id] - ek);
	double xtofs = xtoss; 
	double ytofs = ytoss; 
	double tauxtof = 3.5 * exp(-((y[39][id]) / 30.0) * ((y[39][id]) / 30.0)) + 1.5;
	double tauytof = 20.0 / (1.0 + exp((y[39][id] + 33.5) / 10.0)) + 20.0;
	ydot[10][id] = (xtofs - (xtofs - y[10][id]) * exp(-dt / tauxtof) - y[10][id]) / dt;
	ydot[11][id] = (ytofs - (ytofs - y[11][id]) * exp(-dt / tauytof) - y[11][id]) / dt;
	double I_tof = Kcoeff * GtoFast * y[10][id] * y[11][id] * (y[39][id] - ek);
	double I_to = I_tos + I_tof;
	return I_to;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::i_Cl_lat(int id, double &ecl, double &Fjunc, double &Fsl, double &I_ClCa, double &I_Clbk)
{
	double GClCa = 0.109625;  
	double GClB = 9.0e-3;	  
	double KdClCa = 100.0e-3; 
	double I_ClCa_junc = Fjunc * GClCa / (1.0 + KdClCa / y[36][id]) * (y[39][id] - ecl);
	double I_ClCa_sl = Fsl * GClCa / (1.0 + KdClCa / y[37][id]) * (y[39][id] - ecl);
	I_ClCa = I_ClCa_junc + I_ClCa_sl;
	I_Clbk = GClB * (y[39][id] - ecl);
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::i_LCC_HH_lat(int id, double dt, double &Nao, double &Ko, double &Cao, double &Fjunc_CaL, double &Fsl_CaL, double &I_Ca_junc, double &I_Ca_sl, double &I_CaNa_junc, double &I_CaNa_sl, double &I_CaK)
{
	double K_Ica = 1;		 
	double pNa = K_Ica * 1.5e-8; 
	double pCa = K_Ica * 5.4e-4; 
	double pK = K_Ica * 2.7e-7;	 
	double Q10CaL = 1.8;
	double dss = 1.0 / (1.0 + exp(-(y[39][id] + 14.5) / 6.0));
	double taud = dss * (1.0 - exp(-(y[39][id] + 14.5) / 6.0)) / (0.035 * (y[39][id] + 14.5));
	double fss = 1.0 / (1.0 + exp((y[39][id] + 35.06) / 3.6)) + 0.6 / (1.0 + exp((50 - y[39][id]) / 20));
	double tauf = 1.0 / (0.0197 * exp(-(0.0337 * (y[39][id] + 14.5)) * (0.0337 * (y[39][id] + 14.5))) + 0.02);
	ydot[4][id] = (dss - y[4][id]) / taud;
	ydot[5][id] = (fss - y[5][id]) / (tauf);
	ydot[6][id] = (1.7) * y[36][id] * (1 - y[6][id]) - 11.9e-3 * y[6][id]; 
	ydot[7][id] = 1.7 * y[37][id] * (1 - y[7][id]) - 11.9e-3 * y[7][id];   
	double ibarca_j;
	if (fabs(y[39][id] * 2.0 * FoRT) < 1e-4)
	{
		ibarca_j = pCa * 2.0 * Frdy * (0.341 * y[36][id] * exp(2.0 * y[39][id] * FoRT) - 0.341 * Cao);
	}
	else
	{
		ibarca_j = pCa * 4.0 * (y[39][id] * Frdy * FoRT) * (0.341 * y[36][id] * exp(2.0 * y[39][id] * FoRT) - 0.341 * Cao) / (exp(2.0 * y[39][id] * FoRT) - 1.0);
	}
	double ibarca_sl;
	if (fabs(y[39][id] * 2.0 * FoRT) < 1e-4)
	{
		ibarca_sl = pCa * 2.0 * Frdy * (0.341 * y[37][id] * exp(2.0 * y[39][id] * FoRT) - 0.341 * Cao);
	}
	else
	{
		ibarca_sl = pCa * 4.0 * (y[39][id] * Frdy * FoRT) * (0.341 * y[37][id] * exp(2.0 * y[39][id] * FoRT) - 0.341 * Cao) / (exp(2.0 * y[39][id] * FoRT) - 1.0);
	}
	double ibark;
	if (fabs(y[39][id] * FoRT) < 1e-4)
	{
		ibark = pK * Frdy * (0.75 * y[35][id] * exp(y[39][id] * FoRT) - 0.75 * Ko);
	}
	else
	{
		ibark = pK * (y[39][id] * Frdy * FoRT) * (0.75 * y[35][id] * exp(y[39][id] * FoRT) - 0.75 * Ko) / (exp(y[39][id] * FoRT) - 1.0);
	}
	double ibarna_j;
	if (fabs(y[39][id] * FoRT) < 1e-4)
	{
		ibarna_j = pNa * Frdy * (0.75 * y[32][id] * exp(y[39][id] * FoRT) - 0.75 * Nao);
	}
	else
	{
		ibarna_j = pNa * (y[39][id] * Frdy * FoRT) * (0.75 * y[32][id] * exp(y[39][id] * FoRT) - 0.75 * Nao) / (exp(y[39][id] * FoRT) - 1.0);
	}
	double ibarna_sl;
	if (fabs(y[39][id] * FoRT) < 1e-4)
	{
		ibarna_sl = pNa * Frdy * (0.75 * y[33][id] * exp(y[39][id] * FoRT) - 0.75 * Nao);
	}
	else
	{
		ibarna_sl = pNa * (y[39][id] * Frdy * FoRT) * (0.75 * y[33][id] * exp(y[39][id] * FoRT) - 0.75 * Nao) / (exp(y[39][id] * FoRT) - 1.0);
	}
	double I_Ca_junc1 = (Fjunc_CaL * ibarca_j * y[4][id] * y[5][id] * (1.0 - y[6][id]) * pow(Q10CaL, Qpow)) * 0.45;
	double I_Ca_sl1 = (Fsl_CaL * ibarca_sl * y[4][id] * y[5][id] * (1.0 - y[7][id]) * pow(Q10CaL, Qpow)) * 0.45;
	double I_CaK1 = (ibark * y[4][id] * y[5][id] * (Fjunc_CaL * (1.0 - y[6][id]) + Fsl_CaL * (1.0 - y[7][id])) * pow(Q10CaL, Qpow)) * 0.45;
	double I_CaNa_junc1 = (Fjunc_CaL * ibarna_j * y[4][id] * y[5][id] * (1.0 - y[6][id]) * pow(Q10CaL, Qpow)) * 0.45;
	double I_CaNa_sl1 = (Fsl_CaL * ibarna_sl * y[4][id] * y[5][id] * (1.0 - y[7][id]) * pow(Q10CaL, Qpow)) * 0.45;
	I_Ca_junc = I_Ca_junc1;
	I_Ca_sl = I_Ca_sl1;
	I_CaNa_junc = I_CaNa_junc1;
	I_CaNa_sl = I_CaNa_sl1;
	I_CaK = I_CaK1; 
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::i_LCC_Markov_lat(int id, double dt, double &Nao, double &Ko, double &Cao, double &Fjunc_CaL, double &Fsl_CaL, double &I_Ca_junc, double &I_Ca_sl, double &I_CaNa_junc, double &I_CaNa_sl, double &I_CaK)
{
	double LCC_CKdyadp = y[87 + 6 + 45 + 2][id] / LCCtotDyad; 
	double LCCa_PKAp = y[87 + 6 + 45 + 6 + 28][id] / LCCtotBA;
	double LCCb_PKAp = y[87 + 6 + 45 + 6 + 29][id] / LCCtotBA;
	double K_Ica = 1;		 
	double pNa = K_Ica * 1.5e-8; 
	double pCa = K_Ica * 5.4e-4; 
	double pK = K_Ica * 2.7e-7;	 
	double Q10CaL = 1.8;
	double ibark;
	if (fabs(y[39][id] * FoRT) < 1e-4)
	{
		ibark = pK * Frdy * (0.75 * y[35][id] * exp(y[39][id] * FoRT) - 0.75 * Ko);
	}
	else
	{
		ibark = pK * (y[39][id] * Frdy * FoRT) * (0.75 * y[35][id] * exp(y[39][id] * FoRT) - 0.75 * Ko) / (exp(y[39][id] * FoRT) - 1.0);
	}
	double ibarna_j;
	if (fabs(y[39][id] * FoRT) < 1e-4)
	{
		ibarna_j = pNa * Frdy * (0.75 * y[32][id] * exp(y[39][id] * FoRT) - 0.75 * Nao);
	}
	else
	{
		ibarna_j = pNa * (y[39][id] * Frdy * FoRT) * (0.75 * y[32][id] * exp(y[39][id] * FoRT) - 0.75 * Nao) / (exp(y[39][id] * FoRT) - 1.0);
	}
	double ibarna_sl;
	if (fabs(y[39][id] * FoRT) < 1e-4)
	{
		ibarna_sl = pNa * Frdy * (0.75 * y[33][id] * exp(y[39][id] * FoRT) - 0.75 * Nao);
	}
	else
	{
		ibarna_sl = pNa * (y[39][id] * Frdy * FoRT) * (0.75 * y[33][id] * exp(y[39][id] * FoRT) - 0.75 * Nao) / (exp(y[39][id] * FoRT) - 1.0);
	}
	double cajLCC = y[36][id];
	double caslLCC = y[37][id];
	double taupo = 1.0; 
	double TBa = 450.0; 
	double s1o = 0.0221;
	double k1o = 0.03;
	double kop = 2.5e-3;   
	double cpbar = 8.0e-3; 
	double tca = 78.0312;
	double ICa_scale = 5.25;
	double recoveryReduc = 3.0;
	double fracLCCbp0 = 0.250657;	
	double fracLCCbpISO = 0.525870; 
	double a_favail = (1.56 - 1.0) / (fracLCCbpISO / fracLCCbp0 - 1.0);		
	double favail = (1.0 - a_favail) + a_favail * (LCCb_PKAp / fracLCCbp0); 
	ICa_scale = ICa_scale * favail;
	double SSAshift = 0.0;
	double SSIshift = 0.0;
	double poss = 1.0 / (1.0 + exp(-(y[39][id] + SSAshift) / 8.0));
	double fcaj = 1.0 / (1.0 + pow((kop / cajLCC), 3));
	double Rv = 10.0 + 4954.0 * exp(y[39][id] / 15.6);
	double PrLCC = 1.0 - 1.0 / (1.0 + exp(-(y[39][id] + 40.0) / 4.0));
	double PsLCC = 1.0 / (1.0 + exp(-(y[39][id] + 40.0 + SSIshift) / 11.32));
	double TCaj = (tca + 0.1 * (1.0 + (cajLCC / cpbar) * (cajLCC / cpbar))) / (1.0 + (cajLCC / cpbar) * (cajLCC / cpbar));
	double tauCaj = (Rv - TCaj) * PrLCC + TCaj;
	double tauBa = (Rv - TBa) * PrLCC + TBa;
	double alphaLCC = poss / taupo;
	double betaLCC = (1 - poss) / taupo;
	double r1 = 0.3; 
	double r2 = 3.0; 
	double s1 = s1o * fcaj;
	double s1p = 0.00195; 
	double k1 = k1o * fcaj;
	double k1p = 0.00413; 
	double k2 = 1.0e-4;	  
	double k2p = 0.00224; 
	double s2 = s1 * (k2 / k1) * (r1 / r2);
	double s2p = s1p * (k2p / k1p) * (r1 / r2);
	double k3 = exp(-(y[39][id] + 40.0) / 3.0) / (3.0 * (1 + exp(-(y[39][id] + 40.0) / 3.0)));
	double k3p = k3;
	double k5 = (1.0 - PsLCC) / tauCaj;
	double k6 = (fcaj * PsLCC) / tauCaj;
	double k5p = (1.0 - PsLCC) / tauBa;
	k5 = k5 / recoveryReduc;
	k5p = k5p / recoveryReduc;
	double k6p = PsLCC / tauBa;
	double k4 = k3 * (alphaLCC / betaLCC) * (k1 / k2) * (k5 / k6);
	double k4p = k3p * (alphaLCC / betaLCC) * (k1p / k2p) * (k5p / k6p);
	double Po_LCCj_m1 = 1.0 - y[60][id] - y[61][id] - y[62][id] - y[63][id] - y[64][id] - y[65][id];								  
	ydot[60][id] = betaLCC * y[61][id] + k5 * y[63][id] + k5p * y[65][id] - (k6 + k6p + alphaLCC) * y[60][id];						  
	ydot[61][id] = alphaLCC * y[60][id] + k2 * y[62][id] + k2p * y[64][id] + r2 * Po_LCCj_m1 - (r1 + betaLCC + k1 + k1p) * y[61][id]; 
	ydot[62][id] = k1 * y[61][id] + k4 * y[63][id] + s1 * Po_LCCj_m1 - (k2 + k3 + s2) * y[62][id];									  
	ydot[63][id] = k3 * y[62][id] + k6 * y[60][id] - (k4 + k5) * y[63][id];															  
	ydot[64][id] = k1p * y[61][id] + k4p * y[65][id] + s1p * Po_LCCj_m1 - (k2p + k3p + s2p) * y[64][id];							  
	ydot[65][id] = k3p * y[64][id] + k6p * y[60][id] - (k5p + k4p) * y[65][id];														  
	double ibarca_jm1;
	if (fabs(2 * y[39][id] * FoRT) < 1e-4)
	{
		ibarca_jm1 = (2.0 * pCa * Frdy) * (0.001 * exp(2.0 * y[39][id] * FoRT) - 0.341 * Cao);
	}
	else
	{
		ibarca_jm1 = (4.0 * pCa * y[39][id] * Frdy * FoRT) * (0.001 * exp(2.0 * y[39][id] * FoRT) - 0.341 * Cao) / (exp(2 * y[39][id] * FoRT) - 1.0);
	}
	double I_Ca_junc_m1 = (Fjunc_CaL * ibarca_jm1 * Po_LCCj_m1 * pow(Q10CaL, Qpow)) * ICa_scale;
	double s1om2 = 0.0221;
	double k1om2 = 0.03;
	double kopm2 = 2.5e-3;
	double cpbarm2 = 8.0e-3;
	double tcam2 = 78.0312;
	double possm2 = 1.0 / (1.0 + exp(-(y[39][id] + SSAshift) / 8.0));
	double fcajm2 = 1.0 / (1.0 + pow((kopm2 / cajLCC), 3)); 
	double Rvm2 = 10.0 + 4954.0 * exp(y[39][id] / 15.6);
	double PrLCCm2 = 1.0 - 1.0 / (1.0 + exp(-(y[39][id] + 40.0) / 4.0));
	double PsLCCm2 = 1.0 / (1.0 + exp(-(y[39][id] + 40.0 + SSIshift) / 11.32));
	double TCajm2 = (tcam2 + 0.1 * (1.0 + pow((cajLCC / cpbarm2), 2))) / (1.0 + pow((cajLCC / cpbarm2), 2)); 
	double tauCajm2 = (Rvm2 - TCajm2) * PrLCCm2 + TCajm2;													 
	double tauBam2 = (Rvm2 - TBa) * PrLCCm2 + TBa;
	double alphaLCCm2 = possm2 / taupo;
	double betaLCCm2 = (1 - possm2) / taupo;
	double r1m2 = 0.3;		 
	double r2m2 = 3.0 / 10.0; 
	double s1m2 = s1om2 * fcajm2;
	double s1pm2 = 0.00195; 
	double k1m2 = k1om2 * fcajm2;
	double k1pm2 = 0.00413; 
	double k2m2 = 1.0e-4;	
	double k2pm2 = 0.00224; 
	double s2m2 = s1m2 * (k2m2 / k1m2) * (r1m2 / r2m2);
	double s2pm2 = s1pm2 * (k2pm2 / k1pm2) * (r1m2 / r2m2);
	double k3m2 = exp(-(y[39][id] + 40.0) / 3.0) / (3.0 * (1 + exp(-(y[39][id] + 40.0) / 3.0)));
	double k3pm2 = k3m2;
	double k5m2 = (1.0 - PsLCCm2) / tauCajm2;
	double k6m2 = (fcajm2 * PsLCCm2) / tauCajm2;
	double k5pm2 = (1.0 - PsLCCm2) / tauBam2;
	k5m2 = k5m2 / recoveryReduc;   
	k5pm2 = k5pm2 / recoveryReduc; 
	double k6pm2 = PsLCCm2 / tauBam2;
	double k4m2 = k3m2 * (alphaLCCm2 / betaLCCm2) * (k1m2 / k2m2) * (k5m2 / k6m2);
	double k4pm2 = k3pm2 * (alphaLCCm2 / betaLCCm2) * (k1pm2 / k2pm2) * (k5pm2 / k6pm2);
	double Po_LCCj_m2 = 1.0 - y[66][id] - y[67][id] - y[68][id] - y[69][id] - y[70][id] - y[71][id];												  
	ydot[66][id] = betaLCCm2 * y[67][id] + k5m2 * y[69][id] + k5pm2 * y[71][id] - (k6m2 + k6pm2 + alphaLCCm2) * y[66][id];							  
	ydot[67][id] = alphaLCCm2 * y[66][id] + k2m2 * y[68][id] + k2pm2 * y[70][id] + r2m2 * Po_LCCj_m2 - (r1m2 + betaLCCm2 + k1m2 + k1pm2) * y[67][id]; 
	ydot[68][id] = k1m2 * y[67][id] + k4m2 * y[69][id] + s1m2 * Po_LCCj_m2 - (k2m2 + k3m2 + s2m2) * y[68][id];										  
	ydot[69][id] = k3m2 * y[68][id] + k6m2 * y[66][id] - (k4m2 + k5m2) * y[69][id];																	  
	ydot[70][id] = k1pm2 * y[67][id] + k4pm2 * y[71][id] + s1pm2 * Po_LCCj_m2 - (k2pm2 + k3pm2 + s2pm2) * y[70][id];								  
	ydot[71][id] = k3pm2 * y[70][id] + k6pm2 * y[66][id] - (k5pm2 + k4pm2) * y[71][id];																  
	double ibarca_jm2;
	if (fabs(2 * y[39][id] * FoRT) < 1e-4)
	{
		ibarca_jm2 = (2.0 * pCa * Frdy) * (0.001 * exp(2 * y[39][id] * FoRT) - 0.341 * Cao);
	}
	else
	{
		ibarca_jm2 = (4.0 * pCa * y[39][id] * Frdy * FoRT) * (0.001 * exp(2 * y[39][id] * FoRT) - 0.341 * Cao) / (exp(2 * y[39][id] * FoRT) - 1.0);
	}
	double I_Ca_junc_m2 = (Fjunc_CaL * ibarca_jm2 * (Po_LCCj_m2)*pow(Q10CaL, Qpow)) * ICa_scale;
	double fracLCCap0 = 0.219577; 
	double frac_fpkam2 = (0.15 * fracLCCap0) / (1.0 - fracLCCap0);
	double fpkam2 = (0.15 + frac_fpkam2) * LCCa_PKAp - frac_fpkam2; 
	double fckiim2 = LCC_CKdyadp * 0.1; 
	double junc_mode2 = fckiim2 + fpkam2;
	double I_Ca_junc2 = (1.0 - junc_mode2) * I_Ca_junc_m1 + junc_mode2 * I_Ca_junc_m2;
	double fcasl = 1.0 / (1.0 + pow((kop / caslLCC), 3)); 
	double TCasl = (tca + 0.1 * pow((1 + (caslLCC / cpbar)), 2)) / (1.0 + pow((caslLCC / cpbar), 2));
	double tauCasl = (Rv - TCasl) * PrLCC + TCasl;
	double s1sl = s1o * fcasl;
	double k1sl = k1o * fcasl;
	double s2sl = s1sl * (k2 / k1sl) * (r1 / r2);
	double s2psl = s1p * (k2p / k1p) * (r1 / r2);
	double k5sl = (1.0 - PsLCC) / tauCasl / recoveryReduc; 
	double k6sl = (fcasl * PsLCC) / tauCasl;
	double k4sl = k3 * (alphaLCC / betaLCC) * (k1sl / k2) * (k5sl / k6sl);
	double k4psl = k3p * (alphaLCC / betaLCC) * (k1p / k2p) * (k5p / k6p);
	double Po_LCCsl_m1 = 1.0 - y[72][id] - y[73][id] - y[74][id] - y[75][id] - y[76][id] - y[77][id];									 
	ydot[72][id] = betaLCC * y[73][id] + k5sl * y[75][id] + k5p * y[77][id] - (k6sl + k6p + alphaLCC) * y[72][id];						 
	ydot[73][id] = alphaLCC * y[72][id] + k2 * y[74][id] + k2p * y[76][id] + r2 * Po_LCCsl_m1 - (r1 + betaLCC + k1sl + k1p) * y[73][id]; 
	ydot[74][id] = k1sl * y[73][id] + k4sl * y[75][id] + s1sl * Po_LCCsl_m1 - (k2 + k3 + s2sl) * y[74][id];								 
	ydot[75][id] = k3 * y[74][id] + k6sl * y[72][id] - (k4sl + k5sl) * y[75][id];														 
	ydot[76][id] = k1p * y[73][id] + k4psl * y[77][id] + s1p * Po_LCCsl_m1 - (k2p + k3p + s2psl) * y[76][id];							 
	ydot[77][id] = k3p * y[76][id] + k6p * y[72][id] - (k5p + k4psl) * y[77][id];														 
	double ibarca_slm1;
	if (fabs(2.0 * y[39][id] * FoRT) < 1e-4)
	{
		ibarca_slm1 = (2.0 * pCa * Frdy) * (0.001 * exp(2.0 * y[39][id] * FoRT) - 0.341 * Cao);
	}
	else
	{
		ibarca_slm1 = (4.0 * pCa * y[39][id] * Frdy * FoRT) * (0.001 * exp(2.0 * y[39][id] * FoRT) - 0.341 * Cao) / (exp(2.0 * y[39][id] * FoRT) - 1.0);
	}
	double I_Casl_m1 = (Fsl_CaL * ibarca_slm1 * Po_LCCsl_m1 * pow(Q10CaL, Qpow)) * ICa_scale;
	double r2slm2 = r2m2;
	double s2slm2 = s1sl * (k2 / k1sl) * (r1 / r2slm2);
	double s2pslm2 = s1p * (k2p / k1p) * (r1 / r2slm2);
	double Po_LCCsl_m2 = 1.0 - y[78][id] - y[79][id] - y[80][id] - y[81][id] - y[82][id] - y[83][id]; 
	if (Po_LCCsl_m2 < 0)
		printf("%f\n", Po_LCCsl_m2);
	ydot[78][id] = betaLCC * y[79][id] + k5sl * y[81][id] + k5p * y[83][id] - (k6sl + k6p + alphaLCC) * y[78][id];							 
	ydot[79][id] = alphaLCC * y[78][id] + k2 * y[80][id] + k2p * y[82][id] + r2slm2 * Po_LCCsl_m2 - (r1 + betaLCC + k1sl + k1p) * y[79][id]; 
	ydot[80][id] = k1sl * y[79][id] + k4sl * y[81][id] + s1sl * Po_LCCsl_m2 - (k2 + k3 + s2slm2) * y[80][id];								 
	ydot[81][id] = k3 * y[80][id] + k6sl * y[78][id] - (k4sl + k5sl) * y[81][id];															 
	ydot[82][id] = k1p * y[79][id] + k4psl * y[83][id] + s1p * Po_LCCsl_m2 - (k2p + k3p + s2pslm2) * y[82][id];								 
	ydot[83][id] = k3p * y[82][id] + k6p * y[78][id] - (k5p + k4psl) * y[83][id];															 
	double ibarca_slm2;
	if (fabs(2.0 * y[39][id] * FoRT) < 1e-4)
	{
		ibarca_slm2 = (2.0 * pCa * Frdy) * (0.001 * exp(2.0 * y[39][id] * FoRT) - 0.341 * Cao);
	}
	else
	{
		ibarca_slm2 = (4.0 * pCa * y[39][id] * Frdy * FoRT) * (0.001 * exp(2.0 * y[39][id] * FoRT) - 0.341 * Cao) / (exp(2.0 * y[39][id] * FoRT) - 1.0);
	}
	double I_Casl_m2 = (Fsl_CaL * ibarca_slm2 * Po_LCCsl_m2 * pow(Q10CaL, Qpow)) * ICa_scale;
	double fckiim2_sl = 0.0; 
	double sl_mode2 = fckiim2_sl + fpkam2;
	double I_Ca_sl2 = (1.0 - sl_mode2) * I_Casl_m1 + sl_mode2 * I_Casl_m2;
	double I_CaKj2 = ibark * Fjunc_CaL * ((1.0 - junc_mode2) * Po_LCCj_m1 + junc_mode2 * Po_LCCj_m2) * pow(Q10CaL, Qpow) * ICa_scale;
	double I_CaKsl2 = ibark * Fsl_CaL * ((1.0 - sl_mode2) * Po_LCCsl_m1 + sl_mode2 * Po_LCCsl_m2) * pow(Q10CaL, Qpow) * ICa_scale;
	double I_CaK2 = I_CaKj2 + I_CaKsl2;
	double I_CaNa_junc2 = (Fjunc_CaL * ibarna_j * ((1.0 - junc_mode2) * Po_LCCj_m1 + junc_mode2 * Po_LCCj_m2) * pow(Q10CaL, Qpow)) * ICa_scale;
	double I_CaNa_sl2 = Fsl_CaL * ibarna_sl * ((1.0 - sl_mode2) * Po_LCCsl_m1 + sl_mode2 * Po_LCCsl_m2) * pow(Q10CaL, Qpow) * ICa_scale;
	I_Ca_junc = I_Ca_junc2;
	I_Ca_sl = I_Ca_sl2;
	I_CaNa_junc = I_CaNa_junc2;
	I_CaNa_sl = I_CaNa_sl2;
	I_CaK = I_CaK2; 
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_NCX_lat(int id, double &Nao, double &Cao, double &Fjunc_ncx, double &Fsl_ncx, double &I_ncx_junc, double &I_ncx_sl)
{
	double IbarNCX = 9.0;				 
	double KmCai = 3.59e-3;				 
	double KmCao = 1.3;					 
	double KmNai = 12.29;				 
	double KmNao = 87.5;				 
	double ksat = 0.27;					 
	double nu = 0.35;					 
	double Kdact =  0.256e-3; 
	double Q10NCX = 1.57; 
	double Ka_junc = 1.0 / (1.0 + pow((Kdact / y[36][id]), 3));
	double Ka_sl = 1.0 / (1.0 + pow((Kdact / y[37][id]), 3));
	double s1_junc = exp(nu * y[39][id] * FoRT) * pow(y[32][id], 3) * Cao;
	double s1_sl = exp(nu * y[39][id] * FoRT) * pow(y[33][id], 3) * Cao;
	double s2_junc = exp((nu - 1.0) * y[39][id] * FoRT) * pow(Nao, 3) * y[36][id];
	double s3_junc = (KmCai * pow(Nao, 3) * (1.0 + pow((y[32][id] / KmNai), 3)) + pow(KmNao, 3) * y[36][id] + pow(KmNai, 3) * Cao * (1.0 + y[36][id] / KmCai) + KmCao * pow(y[32][id], 3) + pow(y[32][id], 3) * Cao + pow(Nao, 3) * y[36][id]) * (1.0 + ksat * exp((nu - 1.0) * y[39][id] * FoRT));
	double s2_sl = exp((nu - 1.0) * y[39][id] * FoRT) * pow(Nao, 3) * y[37][id];
	double s3_sl = (KmCai * pow(Nao, 3) * (1.0 + pow((y[33][id] / KmNai), 3)) + pow(KmNao, 3) * y[37][id] + pow(KmNai, 3) * Cao * (1.0 + y[37][id] / KmCai) + KmCao * pow(y[33][id], 3) + pow(y[33][id], 3) * Cao + pow(Nao, 3) * y[37][id]) * (1.0 + ksat * exp((nu - 1) * y[39][id] * FoRT));
	I_ncx_junc = Fjunc_ncx * IbarNCX * pow(Q10NCX, Qpow) * Ka_junc * (s1_junc - s2_junc) / s3_junc;
	I_ncx_sl = Fsl_ncx * IbarNCX * pow(Q10NCX, Qpow) * Ka_sl * (s1_sl - s2_sl) / s3_sl;
	double I_ncx = I_ncx_junc + I_ncx_sl;
	return I_ncx;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_pca_lat(int id, double &Fjunc, double &Fsl, double &I_pca_junc, double &I_pca_sl)
{
	double IbarSLCaP = 0.0673; 
	double Q10SLCaP = 2.35;	   
	double KmPCa = 0.5e-3;	   
	I_pca_junc = Fjunc * pow(Q10SLCaP, Qpow) * IbarSLCaP * pow(y[36][id], 1.6) / (pow(KmPCa, 1.6) + pow(y[36][id], 1.6));
	I_pca_sl = Fsl * pow(Q10SLCaP, Qpow) * IbarSLCaP * pow(y[37][id], 1.6) / (pow(KmPCa, 1.6) + pow(y[37][id], 1.6));
	double I_pca = I_pca_junc + I_pca_sl;
	return I_pca;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_Cabk_lat(int id, double &Fjunc, double &Fsl, double &I_cabk_junc, double &I_cabk_sl, double &eca_junc, double &eca_sl)
{
	double GCaB =  2.513e-4; 
	I_cabk_junc = Fjunc * GCaB * (y[39][id] - eca_junc);
	I_cabk_sl = Fsl * GCaB * (y[39][id] - eca_sl);
	double I_cabk = I_cabk_junc + I_cabk_sl;
	return I_cabk;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::i_RyR_lat(int id, double &J_SRCarel, double &J_SRleak)
{
	double RyR_CKp = y[87 + 6 + 45 + 4][id] / RyRtot; 
	double RyR_PKAp = y[87 + 6 + 45 + 6 + 30][id] / RyRtotBA;
	double ec50SR = 0.45 ; 
	double ks = 25.0;	
	double koCa = 10.0; 
	double kom = 0.06;	
	double kiCa = 0.5;	
	double kim = 0.005; 
	double fCKII_ec50SR = 1.16 - 4.0 / 5.0 * RyR_CKp;
	ec50SR = fCKII_ec50SR * ec50SR; 
	double MaxSR = 15.0;
	double MinSR = 1.0;
	double kCaSR = MaxSR - (MaxSR - MinSR) / (1.0 + pow((ec50SR / y[31][id]), 2.5));
	double koSRCa = koCa / kCaSR;
	double kiSRCa = kiCa * kCaSR;
	double kleak =  5.348e-6; 
	double fCKII_RyR = (20 * RyR_CKp / 3 - 1 / 3); 
	double frac_RyRo = 0.204276;						  
	double a_RyR = (2.0 - 1.0) / (1.0 / frac_RyRo - 1.0); 
	double fPKA_RyR = 1.0 - a_RyR + a_RyR * (RyR_PKAp / frac_RyRo);
	koSRCa = (fCKII_RyR + fPKA_RyR - 1) * koSRCa;
	double RI = 1.0 - y[14][id] - y[15][id] - y[16][id];
	ydot[14][id] = (kim * RI - kiSRCa * y[36][id] * y[14][id]) - (koSRCa * y[36][id] * y[36][id] * y[14][id] - kom * y[15][id]);		
	ydot[15][id] = (koSRCa * y[36][id] * y[36][id] * y[14][id] - kom * y[15][id]) - (kiSRCa * y[36][id] * y[15][id] - kim * y[16][id]); 
	ydot[16][id] = (kiSRCa * y[36][id] * y[15][id] - kim * y[16][id]) - (kom * y[16][id] - koSRCa * y[36][id] * y[36][id] * RI);		
	J_SRCarel = ks * y[15][id] * (y[31][id] - y[36][id]);																				
	if (y[15][id] < 0)
		printf("%f\n", y[15][id]);
	kleak = (1 / 3 + 10 * RyR_CKp / 3)*kleak; 
	J_SRleak = kleak * (y[31][id] - y[36][id]);		   
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_SERCA_lat(int id)
{
	double PLB_CKp = y[87 + 6 + 45 + 5][id] / PLBtot; 
	double PLB_PKAn = (PLBtotBA - y[87 + 6 + 45 + 6 + 26][id]) / PLBtotBA;
	double Q10SRCaP = 2.6;					   
	double Vmax_SRCaP = 1.15 * 1.15 * 2.86e-4; 
	double Kmf = 0.246e-3;					   
	double Kmr = 1.7;						   
	double hillSRCaP = 1.787;				   
	double fCKII_PLB = (1.0 - 0.5 * PLB_CKp); 
	double fracPKA_PLBo = 1.0 - 0.079755;	  
	double fPKA_PLB = (PLB_PKAn / fracPKA_PLBo) * (100.0 - 55.31) / 100.0 + 55.31 / 100.0; 
	if (fCKII_PLB < fPKA_PLB)
	{
		Kmf = Kmf * fCKII_PLB; 
	}
	else if (fPKA_PLB < fCKII_PLB)
	{
		Kmf = Kmf * fPKA_PLB; 
	}
	double J_serca = pow(Q10SRCaP, Qpow) * Vmax_SRCaP * (pow((y[38][id] / Kmf), hillSRCaP) - pow((y[31][id] / Kmr), hillSRCaP)) / (1.0 + pow((y[38][id] / Kmf), hillSRCaP) + pow((y[31][id] / Kmr), hillSRCaP)); 
	return J_serca;
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::Myofilament(int id, int &mechFlag)
{
	double TnI_PKAp = y[87 + 6 + 45 + 6 + 31][id] / TnItotBA;
	double nc = 3.0; 
	double Myo_PKAp = TnI_PKAp;
	double Ap = 1008.0e+04; 
	double Aw = Ap / 5.0;	
	double alfa = 0.5;		
	double bet = 80.0;		
	double Bp = 0.5;		
	double Bw = 0.35;		
	double f = 0.0023;		
	double gama = 28000.0;	
	double hpr = 0.006;		
	double hwr = 0.0001;	
	double Ke = 105000.0;	
	double Lz = 0.97;		
	double La = 1.15;		
	double Lc = 1.05;		
	double Le = 10.0;		
	double RLa = 20.0;		
	double TSt = 0.07 / nc; 
	double Yb = 181.6e+06;
	double Yc = 4.0;		
	double Yd = 0.028;		
	double Yp = 0.1397;		
	double Yr = 0.1397;		
	double Yv = 0.9;		
	double Za = 0.0023;		
	double Zb = 0.1397;		
	double Zp = 0.2095;		
	double Zr = 7262.6e+06; 
	double Fh = 0.1;		
	double fracPKA_Myoo = 0.062698;	  
	double fracPKA_Myoiso = 0.868762; 
	double kPKA_Myo = (Myo_PKAp - fracPKA_Myoo) / (fracPKA_Myoiso - fracPKA_Myoo);
	double uMyo = 1.0;								 
	Ke = (1.0 + uMyo * (0.5 - 1.0) * kPKA_Myo) * Ke; 
	double uXBCa = 1.0;								  
	Zb = (1.0 + uXBCa * (4.2 - 1.0) * kPKA_Myo) * Zb; 
	Zr = (1.0 + uXBCa * (1.8 - 1.0) * kPKA_Myo) * Zr; 
	Yr = (1.0 + uXBCa * (2.2 - 1.0) * kPKA_Myo) * Yr;
	double uXBcy = 1.0; 
	Za = (1.0 + uXBcy * (1.24 - 1.0) * kPKA_Myo) * Za;
	f = (1.0 + uXBcy * (1.24 - 1.0) * kPKA_Myo) * f;
	RLa = (1.0 + uXBcy * (0.4 - 1.0) * kPKA_Myo) * RLa;
	Zp = (1.0 + uXBcy * (2.2 - 1.0) * kPKA_Myo) * Zp;
	Yp = (1.0 + uXBcy * (2.2 - 1.0) * kPKA_Myo) * Yp;
	Bp = (1.0 + uXBcy * (3.4 - 1.0) * kPKA_Myo) * Bp;
	Bw = (1.0 + uXBcy * (3.4 - 1.0) * kPKA_Myo) * Bw;
	Yc = (1.0 + uXBcy * (0.4 - 1.0) * kPKA_Myo) * Yc;
	Yd = (1.0 + uXBcy * (2.2 - 1.0) * kPKA_Myo) * Yd;
	Yv = (1.0 + uXBcy * (1.6 - 1.0) * kPKA_Myo) * Yv;
	double Liso = 1.05;
	double Lm = 1.05 * mechFlag + Liso * (1.0 - mechFlag);
	double Fm = 0.87 * mechFlag + 0.87 * (1.0 - mechFlag);
	double con = 0.0;
	double corL = 10.0;
	double L = 1.03;
	double FB, w, w1;
	if (abs(corL) > 0.00001)
	{
		if (mechFlag == 0)
		{
			Fm = alfa * (exp(bet * (Lm - L)) - 1.0);
		}
		FB = Ap * (y[85 + 4][id] + y[87 + 4][id]) * (L - y[88 + 4][id]) + Aw * y[86 + 4][id] * (L - y[89 + 4][id]);
		w = FB + Ke * pow((L - Lz), 5) + Le * (L - Lz) - Fm;
		w1 = Ap * (y[85 + 4][id] + y[87 + 4][id]) + Aw * y[86 + 4][id] + 5.0 * Ke * pow((L - Lz), 4) + Le + bet * (Fm + alfa);
		corL = -w / w1;
		L = L + 0.1 * corL;
		con = con + 1.0;
	}
	double TS = TSt - y[84 + 4][id] - y[85 + 4][id] - y[86 + 4][id] - y[87 + 4][id];
	double ER = exp(-RLa * (L - La) * (L - La));
	double Yh = Yv * (1 - exp(-gama * (L - y[89 + 4][id] - hwr) * (L - y[89 + 4][id] - hwr)));
	if ((L - y[89 + 4][id]) > hwr)
	{
		Yh = Fh * Yh;
	}
	double fa = f * ER;
	double ga = Za + Yh;
	double gd = Yd * exp(-Yc * (L - Lc));
	ydot[84 + 4][id] = ga * y[86 + 4][id] - fa * y[84 + 4][id] + Yb * TS * pow(y[38][id], nc) - Zb * y[84 + 4][id];			   
	ydot[85 + 4][id] = Zr * y[87 + 4][id] * pow(y[38][id], nc) - Yr * y[85 + 4][id] + Yp * y[86 + 4][id] - Zp * y[85 + 4][id]; 
	ydot[86 + 4][id] = Zp * y[85 + 4][id] - Yp * y[86 + 4][id] + fa * y[84 + 4][id] - ga * y[86 + 4][id];					   
	ydot[87 + 4][id] = -gd * y[87 + 4][id] + Yr * y[85 + 4][id] - Zr * y[87 + 4][id] * pow(y[38][id], nc);					   
	ydot[88 + 4][id] = Bp * (L - y[88 + 4][id] - hpr);																		   
	ydot[89 + 4][id] = Bw * (L - y[89 + 4][id] - hwr);																		   
	if (mechFlag == 1)
	{
		Lm = L + log((Fm + alfa) / alfa) / bet;
	}
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ void myocyte<NMYOS, NCable>::Na_Ca_buffer(int id, double &Mgi, double &Vmyo, double &Vjunc, double &Vsl, double &J_CaB_cytosol, double &J_CaB_junction, double &J_CaB_sl)
{
	double nc = 3.0; 
	double Bmax_Naj = 7.561;							
	double Bmax_Nasl = 1.65;							
	double koff_na = 1.0e-3;							
	double kon_na = 0.1e-3;								
	double Bmax_TnChigh = 140.0e-3;						
	double koff_tnchca = 0.032e-3;						
	double kon_tnchca = 2.37;							
	double koff_tnchmg = 3.33e-3;						
	double kon_tnchmg = 3.0e-3;							
	double Bmax_NMYOSin = 140.0e-3;						
	double koff_myoca = 0.46e-3;						
	double kon_myoca = 13.8;							
	double koff_myomg = 0.057e-3;						
	double kon_myomg = 0.0157;							
	double Bmax_SR = 19.0 * 0.9e-3;						
	double koff_sr = 60.0e-3;							
	double kon_sr = 100.0;								
	double Bmax_SLlowsl = 37.38e-3 * Vmyo / Vsl;		
	double Bmax_SLlowj = 4.62e-3 * Vmyo / Vjunc * 0.1;	
	double koff_sll = 1300.0e-3;						
	double kon_sll = 100.0;								
	double Bmax_SLhighsl = 13.35e-3 * Vmyo / Vsl;		
	double Bmax_SLhighj = 1.65e-3 * Vmyo / Vjunc * 0.1; 
	double koff_slh = 30.0e-3;							
	double kon_slh = 100.0;								
	double PLB_CKp = y[87 + 6 + 45 + 5][id] / PLBtot; 
	double PLB_PKAn = (PLBtotBA - y[87 + 6 + 45 + 6 + 26][id]) / PLBtotBA;
	double fCKII_PLB = (1.0 - 0.5 * PLB_CKp);											   
	double fracPKA_PLBo = 1.0 - 0.079755;												   
	double fPKA_PLB = (PLB_PKAn / fracPKA_PLBo) * (100.0 - 55.31) / 100.0 + 55.31 / 100.0; 
	if (fCKII_PLB < fPKA_PLB)
	{
		koff_sr = koff_sr * fCKII_PLB;
	}
	else if (fPKA_PLB < fCKII_PLB)
	{
		koff_sr = koff_sr * fPKA_PLB;
	}
	ydot[17][id] = kon_na * y[32][id] * (Bmax_Naj - y[17][id]) - koff_na * y[17][id];  
	ydot[18][id] = kon_na * y[33][id] * (Bmax_Nasl - y[18][id]) - koff_na * y[18][id]; 
	ydot[19][id] = nc * (ydot[84 + 4][id] + ydot[85 + 4][id] + ydot[86 + 4][id]);							  
	ydot[20][id] = kon_tnchca * y[38][id] * (Bmax_TnChigh - y[20][id] - y[21][id]) - koff_tnchca * y[20][id]; 
	ydot[21][id] = kon_tnchmg * Mgi * (Bmax_TnChigh - y[20][id] - y[21][id]) - koff_tnchmg * y[21][id];		  
	ydot[22][id] = 0.0;																						  
	ydot[23][id] = kon_myoca * y[38][id] * (Bmax_NMYOSin - y[23][id] - y[24][id]) - koff_myoca * y[23][id]; 
	ydot[24][id] = kon_myomg * Mgi * (Bmax_NMYOSin - y[23][id] - y[24][id]) - koff_myomg * y[24][id];		
	ydot[25][id] = kon_sr * y[38][id] * (Bmax_SR - y[25][id]) - koff_sr * y[25][id];						
	J_CaB_cytosol = ydot[19][id] + ydot[20][id] + ydot[22][id] + ydot[23][id] + ydot[25][id];
	ydot[26][id] = kon_sll * y[36][id] * (Bmax_SLlowj - y[26][id]) - koff_sll * y[26][id];	 
	ydot[27][id] = kon_sll * y[37][id] * (Bmax_SLlowsl - y[27][id]) - koff_sll * y[27][id];	 
	ydot[28][id] = kon_slh * y[36][id] * (Bmax_SLhighj - y[28][id]) - koff_slh * y[28][id];	 
	ydot[29][id] = kon_slh * y[37][id] * (Bmax_SLhighsl - y[29][id]) - koff_slh * y[29][id]; 
	J_CaB_junction = ydot[26][id] + ydot[28][id];
	J_CaB_sl = ydot[27][id] + ydot[29][id];
}
template <unsigned int NMYOS, unsigned int NCable>
__host__ __device__ double myocyte<NMYOS, NCable>::i_ChR2(int id, double &Ilight)
{
	double E_ChR2 = 0.0;
	double g_ChR2 = 0.75; 
	double gama0 = 0.1;
	double Gd2 = 0.05;
	double epsilon1 = 0.8535;
	double epsilon2 = 0.14;
	double tau_ChR2 = 1.3;
	double e12d = 0.011;
	double e21d = 0.008;
	double lambda = 470;
	double w_loss = 0.77;
	double G = (10.6408 - 14.6408 * exp(-y[39][id] / 42.7671)) / y[39][id];
	double Gd1 = (0.075 + 0.043 * tanh((y[39][id] + 20) / -20));
	double Gr = 4.34587 * 1.0e-5 * exp(-0.0211539274 * y[39][id]);
	double e12 = e12d + 0.005 * log(1 + Ilight / 0.024);
	double e21 = e21d + 0.004 * log(1 + Ilight / 0.024);
	double sita = 100 * Ilight;
	double S0_sita = 0.5 * (1 + tanh(120 * (sita - 0.1)));
	double F = 0.0006 * Ilight * lambda / w_loss;
	double k1_ChR2 = epsilon1 * F * p_ChR2_myo[id];
	double k2_ChR2 = epsilon2 * F * p_ChR2_myo[id];
	d_O1_ChR2_myo[id] = (-(Gd1 + e12) * O1_ChR2_myo[id] + e21 * O2_ChR2_myo[id] + k1_ChR2 * C1_ChR2_myo[id]);
	d_O2_ChR2_myo[id] = (e12 * O1_ChR2_myo[id] - (Gd2 + e21) * O2_ChR2_myo[id] + k2_ChR2 * C2_ChR2_myo[id]);
	d_C1_ChR2_myo[id] = (Gd1 * O1_ChR2_myo[id] - k1_ChR2 * C1_ChR2_myo[id] + Gr * C2_ChR2_myo[id]);
	d_C2_ChR2_myo[id] = (Gd2 * O2_ChR2_myo[id] - (k2_ChR2 + Gr) * C2_ChR2_myo[id]);
	d_p_ChR2_myo[id] = (S0_sita - p_ChR2_myo[id]) / tau_ChR2;
	double I_ChR2_myo = g_ChR2 * G * (O1_ChR2_myo[id] + gama0 * O2_ChR2_myo[id]) * (y[39][id] - E_ChR2); 
	return I_ChR2_myo;
}
