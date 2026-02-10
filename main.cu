#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <fstream>
#include <vector>
#include <sstream> 
#include <stdlib.h>
#include <cublas_v2.h>          
#include <sys/time.h>
#include <time.h>
using namespace std;
#include "SpatialTissue.cuh"
#define N_PARAMS 10 
#define INPUTFROMFILE
const unsigned int NBeat = 2;
const double StimStrength = 50;
const double StimDur = 1;
double dt = 0.001;
const unsigned int StimSteps = static_cast<int>(StimDur / dt);
const double MAXDVDT = 25;
const unsigned int TS_ADAPT_MULT = 100; 
const unsigned int NUMofMYOS = number_of_myocyte;       
const unsigned int NUMofCable = number_of_1D_cables;     
const int current_flag = 0;     
const int parameter_flag = 0;     
int main(int argc, char *argv[])
{
	cout << "number_of_myocyte: " << number_of_myocyte << endl;
	cout << "number_of_1D_cables: " << number_of_1D_cables << endl;
	cout << endl;
	cout << endl;
	std::string stim_location = argv[1];
	std::string spcl = "1000";
	std::string sg_cleft_control = "3";
	std::string swidth_of_cleft_myo = "18";
	std::string sG_gap_myo = "0";
	std::string srate_aj_af = "0.1";	                           
	std::string srate_af_am = "0.2";	                           
	std::string sNUMofMYOS = std::to_string(NUMofMYOS);            
	std::string sNUMofCable = std::to_string(NUMofCable);          
	std::vector<int> cell_type_and_location(NUMofMYOS * NUMofCable, 1);
	int id, id_cell, id_cleft;
	for (int idy = 0; idy < NUMofCable; idy++){
		cout<<"Cell location for Cable "<< idy + 1 <<" :"<<endl;
		for(int idx = 0; idx < NUMofMYOS; idx++){
			id = idx + idy * NUMofMYOS;
			cout << cell_type_and_location[id] << " ";
		}
		cout << endl;
		cout << endl;
	}
	cout << endl;
	int current_stim_location = std::stoi(stim_location) - 1;
	class SpatialTissue<NUMofMYOS, NUMofCable> *wholeTissue = 
		new SpatialTissue<NUMofMYOS, NUMofCable>(current_stim_location);
	double pcl = atof(spcl.c_str());
	double rate_af_am = atof(srate_af_am.c_str());
	double rate_aj_af = atof(srate_aj_af.c_str());
	double g_cleft_control[NUMofCable];
	double width_of_cleft_myo[NUMofCable];
	double G_gap_myo[NUMofCable];
	FILE *outputFile[NUMofCable];
		int idy = 0;
		std::string widthinput = argv[2];
		std::string gapinput = argv[3];
		g_cleft_control[idy] = 3.0f;
		width_of_cleft_myo[idy] = std::stof(widthinput);
		G_gap_myo[idy] = std::stof(gapinput) * 1e-6;
		std::string filename = "s"+stim_location+"_NofM_"+ sNUMofMYOS +
							"_PCL_" + spcl +
							"_width_of_cleft_myo_" + widthinput + 
							"_G_gap_myo_" + gapinput +
							".txt";
		outputFile[idy] = fopen(filename.c_str(), "w");
	wholeTissue -> arrange_cells_surface_and_Capacitance(rate_af_am, rate_aj_af);
	wholeTissue -> arrange_cleft_parameter(g_cleft_control, width_of_cleft_myo, G_gap_myo);
	wholeTissue -> calculate_invert_cublas();              
	wholeTissue -> cpu_copy_gpu();	
	int myo_number_counter;
	for (int idy = 0; idy < NUMofCable; idy++){
		myo_number_counter = 0;
		fprintf(outputFile[idy], "t ");
		for (int idx = 0; idx < NUMofMYOS - 1; idx++){
			fprintf(outputFile[idy], "myo_%d cleft_%d ", myo_number_counter + 1, idx + 1);
			myo_number_counter++;
		}
		fprintf(outputFile[idy], "myo_%d", myo_number_counter + 1);
		if(parameter_flag){
			for (int idx = 0; idx < NUMofMYOS - 2; idx++){
				fprintf(outputFile[idy], " Na_cleft_%d K_cleft_%d", idx + 1, idx + 1);
			}
			fprintf(outputFile[idy], " Na_cleft_%d K_cleft_%d", NUMofMYOS - 2 + 1, NUMofMYOS - 2 + 1);
		}
		if(current_flag){
			myo_number_counter = 0;
			for (int idx = 0; idx < NUMofMYOS - 1; idx++){
				fprintf(outputFile[idy], " myo_%d_lat_current myo_%d_lat_Na_current myo_%d_lat_K1_current myo_%d_left_Na_current myo_%d_right_Na_current myo_%d_left_K1_current myo_%d_right_K1_current I_Na_cleft_%d I_K1_cleft_%d",
										myo_number_counter + 1,myo_number_counter + 1,myo_number_counter + 1, myo_number_counter + 1, myo_number_counter + 1, myo_number_counter + 1, myo_number_counter + 1,
										idx + 1, idx + 1);
				myo_number_counter++;
			}
			fprintf(outputFile[idy], " myo_%d_lat_current myo_%d_lat_Na_current myo_%d_lat_K1_current myo_%d_left_Na_current myo_%d_right_Na_current myo_%d_left_K1_current myo_%d_right_K1_current", 
									myo_number_counter + 1,myo_number_counter + 1,myo_number_counter + 1, myo_number_counter + 1, myo_number_counter + 1, myo_number_counter + 1, myo_number_counter + 1);
		}
		fprintf(outputFile[idy], "\n");
	}
	double time = 0.0;
	struct timeval startTime, endTime;
	gettimeofday(&startTime, NULL);
	for (int numbeat = 0; numbeat < NBeat; numbeat++)
	{
		double t0 = pcl * numbeat; 
		int Nsteps = static_cast<int>(pcl / dt);
		for (int i = 0; i < Nsteps; i++)
		{
			time = t0 + i * dt;
			bool output_time_check;
			if (i > 0 && i % 100 == 0)
			{
				output_time_check = true;
				for (int idy = 0; idy < NUMofCable; idy++){
					fprintf(outputFile[idy], "%f", time);
					for (int idx = 0; idx < NUMofMYOS - 1; idx++)
					{
						id_cell = idx + idy * NUMofMYOS;
						id_cleft = idx + idy * (NUMofMYOS - 1);
						fprintf(outputFile[idy], " %f %f", wholeTissue -> cable_grid -> V_Cell[id_cell], 
														wholeTissue -> cable_grid -> V_cleft[id_cleft]);
					}
					fprintf(outputFile[idy], " %f", wholeTissue -> cable_grid -> V_Cell[NUMofMYOS - 1 + NUMofMYOS * idy]);
					if(parameter_flag){
						for (int idx = 0; idx < NUMofMYOS - 2; idx++)
						{
							id_cleft = idx + idy * (NUMofMYOS - 1);
							fprintf(outputFile[idy], " %f %f", wholeTissue -> cable_grid -> Na_cleft[id_cleft],
																wholeTissue -> cable_grid -> K_cleft[id_cleft]);
						}
						id_cleft = NUMofMYOS - 2 + idy * (NUMofMYOS - 1);
						fprintf(outputFile[idy], " %f %f", wholeTissue -> cable_grid -> Na_cleft[id_cleft],
															wholeTissue -> cable_grid -> K_cleft[id_cleft]);
					}
					if(current_flag){
						for (int idx = 0; idx < NUMofMYOS - 1; idx++){
							id_cell = idx + idy * NUMofMYOS;
							id_cleft = idx + idy * (NUMofMYOS - 1);
							fprintf(outputFile[idy], " %e %e %e %e %e %e %e %e %e",
													wholeTissue -> cable_grid -> Cell_lat_current[id_cell],
													wholeTissue -> cable_grid -> Cell_lat_Ina_current_foroutput[id_cell],
													wholeTissue -> cable_grid -> Cell_lat_IK1_current_foroutput[id_cell],
													wholeTissue -> cable_grid -> Cell_Na_junc_current[0][id_cell],
													wholeTissue -> cable_grid -> Cell_Na_junc_current[1][id_cell],
													wholeTissue -> cable_grid -> Cell_K1_junc_current[0][id_cell],
													wholeTissue -> cable_grid -> Cell_K1_junc_current[1][id_cell], 
													wholeTissue -> cable_grid -> I_Na_cleft[id_cleft],
													wholeTissue -> cable_grid -> I_K_cleft[id_cleft]);
						}
						fprintf(outputFile[idy], " %e %e %e %e %e %e %e",
												wholeTissue -> cable_grid -> Cell_lat_current[NUMofMYOS - 1 + NUMofMYOS * idy],
												wholeTissue -> cable_grid -> Cell_lat_Ina_current_foroutput[NUMofMYOS - 1 + NUMofMYOS * idy],
												wholeTissue -> cable_grid -> Cell_lat_IK1_current_foroutput[NUMofMYOS - 1 + NUMofMYOS * idy],
												wholeTissue -> cable_grid -> Cell_Na_junc_current[0][NUMofMYOS - 1 + NUMofMYOS * idy],
												wholeTissue -> cable_grid -> Cell_Na_junc_current[1][NUMofMYOS - 1 + NUMofMYOS * idy],
												wholeTissue -> cable_grid -> Cell_K1_junc_current[0][NUMofMYOS - 1 + NUMofMYOS * idy],
												wholeTissue -> cable_grid -> Cell_K1_junc_current[1][NUMofMYOS - 1 + NUMofMYOS * idy]);
					}
					fprintf(outputFile[idy], "\n");
				}
			}
			else if(i > 0 && i % 1 == 0 && i < StimSteps * 30){
				output_time_check = true;
				for (int idy = 0; idy < NUMofCable; idy++){
					fprintf(outputFile[idy], "%f", time);
					for (int idx = 0; idx < NUMofMYOS - 1; idx++)
					{
						id_cell = idx + idy * NUMofMYOS;
						id_cleft = idx + idy * (NUMofMYOS - 1);
						fprintf(outputFile[idy], " %f %f", wholeTissue -> cable_grid -> V_Cell[id_cell], 
														wholeTissue -> cable_grid -> V_cleft[id_cleft]);
					}
					fprintf(outputFile[idy], " %f", wholeTissue -> cable_grid -> V_Cell[NUMofMYOS - 1 + NUMofMYOS * idy]);
					if(parameter_flag){
						for (int idx = 0; idx < NUMofMYOS - 2; idx++)
						{
							id_cleft = idx + idy * (NUMofMYOS - 1);
							fprintf(outputFile[idy], " %f %f", wholeTissue -> cable_grid -> Na_cleft[id_cleft],
																wholeTissue -> cable_grid -> K_cleft[id_cleft]);
						}
						id_cleft = NUMofMYOS - 2 + idy * (NUMofMYOS - 1);
						fprintf(outputFile[idy], " %f %f", wholeTissue -> cable_grid -> Na_cleft[id_cleft],
															wholeTissue -> cable_grid -> K_cleft[id_cleft]);
					}
					if(current_flag){
						for (int idx = 0; idx < NUMofMYOS - 1; idx++){
							id_cell = idx + idy * NUMofMYOS;
							id_cleft = idx + idy * (NUMofMYOS - 1);
							fprintf(outputFile[idy], " %e %e %e %e %e %e %e %e %e",
													wholeTissue -> cable_grid -> Cell_lat_current[id_cell],
													wholeTissue -> cable_grid -> Cell_lat_Ina_current_foroutput[id_cell],
													wholeTissue -> cable_grid -> Cell_lat_IK1_current_foroutput[id_cell],
													wholeTissue -> cable_grid -> Cell_Na_junc_current[0][id_cell],
													wholeTissue -> cable_grid -> Cell_Na_junc_current[1][id_cell],
													wholeTissue -> cable_grid -> Cell_K1_junc_current[0][id_cell],
													wholeTissue -> cable_grid -> Cell_K1_junc_current[1][id_cell], 
													wholeTissue -> cable_grid -> I_Na_cleft[id_cleft],
													wholeTissue -> cable_grid -> I_K_cleft[id_cleft]);
						}
					fprintf(outputFile[idy], " %e %e %e %e %e %e %e",
												wholeTissue -> cable_grid -> Cell_lat_current[NUMofMYOS - 1 + NUMofMYOS * idy],
												wholeTissue -> cable_grid -> Cell_lat_Ina_current_foroutput[NUMofMYOS - 1 + NUMofMYOS * idy],
												wholeTissue -> cable_grid -> Cell_lat_IK1_current_foroutput[NUMofMYOS - 1 + NUMofMYOS * idy],
												wholeTissue -> cable_grid -> Cell_Na_junc_current[0][NUMofMYOS - 1 + NUMofMYOS * idy],
												wholeTissue -> cable_grid -> Cell_Na_junc_current[1][NUMofMYOS - 1 + NUMofMYOS * idy],
												wholeTissue -> cable_grid -> Cell_K1_junc_current[0][NUMofMYOS - 1 + NUMofMYOS * idy],
												wholeTissue -> cable_grid -> Cell_K1_junc_current[1][NUMofMYOS - 1 + NUMofMYOS * idy]);
					}
					fprintf(outputFile[idy], "\n");
				}
			}
			else{
				output_time_check = false;
			}
			if(i == 0) output_time_check = true;
			double Istim = (i < StimSteps && numbeat < NBeat) ? StimStrength : 0.0;
			double in_time = i*dt; 
			bool success = wholeTissue -> step(dt, Istim, MAXDVDT, output_time_check,in_time);
			if (!success)
			{
				for (int n = 0; n < TS_ADAPT_MULT; n++)
				{
					wholeTissue -> step(dt / TS_ADAPT_MULT, Istim, MAXDVDT, output_time_check,in_time);
				}
			}
		}
	}
	gettimeofday(&endTime, NULL);
	printf("Gpu use time: %ld s\n",
		(endTime.tv_sec - startTime.tv_sec));
	delete wholeTissue;
	return 0;
 }
