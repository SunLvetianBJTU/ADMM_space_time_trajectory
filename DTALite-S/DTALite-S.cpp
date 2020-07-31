//  Portions Copyright 2010 Xuesong Zhou

//   If you help write or modify the code, please also list your names here.
//   The reason of having copyright info here is to ensure all the modified version, as a whole, under the GPL 
//   and further prevent a violation of the GPL.

// More about "How to use GNU licenses for your own software"
// http://www.gnu.org/licenses/gpl-howto.html

//    This file is part of DTALite.

//    DTALite is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    DTALite is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with DTALite.  If not, see <http://www.gnu.org/licenses/>.

// DTALite.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <list> 
#include <omp.h>
#include <algorithm>
#include <time.h>
#include "CSVParser.h"
#include "DTALite-S.h"
#include <functional>
#include<stdio.h>   
#include<tchar.h>
#define _MAX_LABEL_COST 999999999
#define _MAX_NUMBER_OF_PAX 100
#define _MAX_NUMBER_OF_VEHICLES 100
#define _MAX_NUMBER_OF_TIME_INTERVALS 150
#define _MAX_NUMBER_OF_DEMAND_TYPES 20
#define _MAX_NUMBER_OF_PHYSICAL_NODES 10000

#define _MAX_STATES 2
// Linear congruential generator 
#define LCG_a 17364
#define LCG_c 0
#define LCG_M 65521  // it should be 2^32, but we use a small 16-bit number to save memory

// The one and only application object

CWinApp theApp;
using namespace std;
TCHAR g_SettingFileName[_MAX_PATH] = _T("./Settings.txt");

FILE* g_pFileDebugLog = NULL;

FILE* g_pFileOutputLog = NULL;
FILE* g_pFileDebugLog_LR = NULL;
FILE* g_pFileDebugLog_ADMM = NULL;
FILE* g_pTSViewOutput = NULL;
FILE* g_pNGSIMOuputLog = NULL;
FILE* g_ptrainDelayLog = NULL;

int enum_waiting_link_type = 5;  //To do: to be changed. 
int enum_road_capacity_link_type = 0;  //To do: to be changed. 
int enum_request_link_type = 100;  //To do: to be changed. 

int g_number_of_threads = 4;
int g_shortest_path_debugging_flag = 0;
int g_number_of_agents;
int g_number_of_demand_types = 1;

int time_index = 0;

double g_number_of_seconds_per_interval = 6;  // 0.2 seconds for 300 intervals per min
int g_number_of_simulation_intervals = 600 * 60 / g_number_of_seconds_per_interval;    // 60min
int g_number_of_optimization_time_intervals = 60;

int g_Simulation_StartTimeInMin = 9999;
int g_Simulation_EndTimeInMin = 0;
int g_start_simu_interval_no, g_end_simu_interval_no;

int g_Post_Simulation_DurationInMin = 120;
int g_dp_algorithm_debug_flag = 0;
float g_penalty_RHO = 0.5;
int g_optimization_method = 2;
float LR_multiplier;
//g_convert_abs_simu_interval_to_relative_simu_interval
int g_A2R_simu_interval(int abs_simu_interval)
{
	return abs_simu_interval - g_start_simu_interval_no;

}

std::map<int, int> g_link_key_to_seq_no_map;  // hush table, map key to internal link sequence no. 

int g_TAU;


std::map<int, int> g_internal_node_seq_no_map;  // hush table, map external node number to internal node sequence no. 
std::map<int, int> g_internal_node_seq_no_to_node_id_map;  // hush table, map external node number to internal node sequence no. 


long g_GetLinkSeqNo(int from_node_id, int to_node_id)
{
	if (g_internal_node_seq_no_map.find(from_node_id) == g_internal_node_seq_no_map.end())
	{
		return -1; //have not been defined
	}

	if (g_internal_node_seq_no_map.find(from_node_id) == g_internal_node_seq_no_map.end())
	{
		return -1; //have not been defined
	}

	int from_node_seq_no = g_internal_node_seq_no_map[from_node_id];
	int to_node_seq_no = g_internal_node_seq_no_map[to_node_id];

	long link_key = from_node_seq_no * _MAX_NUMBER_OF_PHYSICAL_NODES + to_node_seq_no;

	if (g_link_key_to_seq_no_map.find(link_key) != g_link_key_to_seq_no_map.end())
		return g_link_key_to_seq_no_map[link_key];
	else
		return -1;
}

class CLink
{
public:
	CLink()  // construction 
	{
		//m_End_of_RedTime_in_simu_interval = 0;
		link_flag = 0;
		external_link_id = 0;
		service_type = 0;
		service_price = 0;
		VRP_load_id = -1;
		base_price = 20;
		VRP_load_difference = 0;
		LR_multiplier = 0;
		cost = 0;
		BRP_alpha = 0.15f;
		BRP_beta = 4.0f;
		link_capacity = 1000;
		free_flow_travel_time_in_min = 1;
		flow_volume = 0;
		number_of_lanes = 1;
		// mfd
		mfd_zone_id = 0;

		max_allowed_waiting_time = 0;

		VRP_time_window_begin = -1;
		VRP_time_window_end = 10000;
	}

	~CLink()
	{
		DeallocateMemory();
	}

	std::list<int>  m_waiting_traveler_queue;

	// all alocated as relative time
	float* m_LinkOutFlowCapacity;
	float* m_LinkInFlowCapacity;
	int m_End_of_RedTime_in_simu_interval;

	int* m_LinkCumulativeArrival;
	int* m_LinkCumulativeDeparture;
	int* m_LinkCumulativeVirtualDelay;
	float* m_LinkTravelTime;

	int m_CumulativeArrivalCount;
	int m_CumulativeDepartureCount;
	int m_CumulativeVirtualDelayCount;

	float LR_multiplier;
	float link_cost;

	int sub_flag;

	float free_speed;
	int VRP_time_window_begin, VRP_time_window_end;
	float base_price;
	std::vector<float> travel_time_vector;
	std::vector<float> time_dependent_LR_multiplier_vector, time_dependent_external_cost_vector, time_dependent_ADMM_multiplier_vector, time_dependent_discharge_rate, time_dependent_inflow_rate;

	std::vector<int> time_dependent_visit_counts, time_dependent_ADMM_visit_counts,
		time_depedent_capacity_vector;
	std::vector<float> time_dependent_link_cost;
	std::vector<float> time_dependent_travel_time_vector;

	//ADMM_multiplier_matrix is used in the searching process and updated by LR_multiplier_matrix
	int max_allowed_waiting_time;

	void Setup_State_Dependent_Data_Matrix(int number_of_optimization_time_intervals)
	{

		for (int t = 0; t < number_of_optimization_time_intervals; t++)
		{
			time_dependent_visit_counts.push_back(0);
			time_dependent_LR_multiplier_vector.push_back(0);
			time_depedent_capacity_vector.push_back(link_capacity);

			time_dependent_external_cost_vector.push_back(0);
			time_dependent_ADMM_multiplier_vector.push_back(0);
			travel_time_vector.push_back((int)(free_flow_travel_time_in_min));  //assume simulation time interval as free-flow travel time per cell 

		}

		VRP_time_window_begin = max(0, VRP_time_window_begin);
		VRP_time_window_end = min(number_of_optimization_time_intervals - 1, VRP_time_window_end);

	}

	float GetCapacityPerSimuInterval(float link_capacity_per_hour)
	{
		return link_capacity_per_hour / 3600.0 *g_number_of_seconds_per_interval;
	}

	void AllocateMemory()
	{
		m_LinkOutFlowCapacity = new float[g_number_of_simulation_intervals];
		m_LinkInFlowCapacity = new float[g_number_of_simulation_intervals];
		m_LinkCumulativeArrival = new int[g_number_of_simulation_intervals];
		m_LinkCumulativeDeparture = new int[g_number_of_simulation_intervals];
		m_LinkCumulativeVirtualDelay = new int[g_number_of_simulation_intervals];
		m_LinkTravelTime = new float[g_number_of_simulation_intervals];


		for (int t = 0; t < g_number_of_simulation_intervals; t++)
		{
			m_LinkOutFlowCapacity[t] = GetCapacityPerSimuInterval(link_capacity);
			m_LinkCumulativeArrival[t] = 0;
			m_LinkCumulativeDeparture[t] = 0;
			m_LinkCumulativeVirtualDelay[t] = 0;

		}

		free_flow_travel_time_in_simu_interval = int(free_flow_travel_time_in_min*60.0 / g_number_of_seconds_per_interval + 0.5);
	}

	void ResetMOE()
	{
		m_CumulativeArrivalCount = 0;
		m_CumulativeDepartureCount = 0;
		m_CumulativeVirtualDelayCount = 0;


	}

	void DeallocateMemory()
	{
		//if(m_LinkOutFlowCapacity != NULL) delete m_LinkOutFlowCapacity;
		//if (m_LinkInFlowCapacity != NULL) delete m_LinkInFlowCapacity;
		//if (m_LinkCumulativeArrival != NULL) delete m_LinkCumulativeArrival;
		//if (m_LinkCumulativeDeparture != NULL) delete m_LinkCumulativeDeparture;
		//if (m_LinkTravelTime != NULL) delete m_LinkTravelTime;

	}
	int link_flag;
	int external_link_id;
	int link_seq_no;  // internal seq no
	int from_node_seq_no;
	int to_node_seq_no;
	float cost;
	float free_flow_travel_time_in_min;
	int free_flow_travel_time_in_simu_interval;
	int number_of_lanes;
	bool demand_type_code[_MAX_NUMBER_OF_DEMAND_TYPES];
	float demand_type_TTcost[_MAX_NUMBER_OF_DEMAND_TYPES];

	int type;
	int service_type; // 0: moving, -1: drop off, +1, pick up

	float service_price; // for pick up or drop off
	int VRP_load_id;
	int VRP_group_id;

	int VRP_load_difference; // we use a single point time window now

	int link_capacity;
	float flow_volume;
	float travel_time;
	float BRP_alpha;
	float BRP_beta;
	float length;
	// mfd
	int mfd_zone_id;

	void CalculateBPRFunctionAndCost()
	{
		travel_time = free_flow_travel_time_in_min*(1 + BRP_alpha*pow(flow_volume / max(0.00001, link_capacity), BRP_beta));
		cost = travel_time;
	}

	float get_VOC_ratio()
	{
		return flow_volume / max(0.00001, link_capacity);

	}

	float get_speed()
	{
		return length / max(travel_time, 0.0001) * 60;  // per hour
	}
	// mfd 

	float get_link_in_flow_per_min(int time_in_min)
	{
		if (time_in_min < g_Simulation_EndTimeInMin + g_Post_Simulation_DurationInMin - 1)
		{
			int next_time_in_internval = (time_in_min + 1) * 60 / g_number_of_seconds_per_interval;
			int time_in_interval = time_in_min * 60 / g_number_of_seconds_per_interval;
			return m_LinkCumulativeArrival[g_A2R_simu_interval(next_time_in_internval)] - m_LinkCumulativeArrival[g_A2R_simu_interval(time_in_interval)];
		}
		else
			return 0;

	}

	float get_link_out_flow_per_min(int time_in_min)
	{
		if (time_in_min < g_Simulation_EndTimeInMin + g_Post_Simulation_DurationInMin - 1)
		{
			int next_time_in_internval = (time_in_min + 1) * 60 / g_number_of_seconds_per_interval;
			int time_in_interval = time_in_min * 60 / g_number_of_seconds_per_interval;

			return m_LinkCumulativeDeparture[g_A2R_simu_interval(next_time_in_internval)] - m_LinkCumulativeDeparture[g_A2R_simu_interval(time_in_interval)];
		}
		else
			return 0;

	}

	float get_number_of_vehicles(int time_in_min)
	{
		if (time_in_min < g_Simulation_EndTimeInMin + g_Post_Simulation_DurationInMin)
		{
			int time_in_interval = g_A2R_simu_interval((time_in_min) * 60 / g_number_of_seconds_per_interval);


			return m_LinkCumulativeArrival[time_in_interval] - m_LinkCumulativeDeparture[time_in_interval];
		}
		else
			return 0;

	}


	float get_avg_delay_in_min(int time_in_min, int time_duration)
	{

		if (m_LinkCumulativeVirtualDelay != NULL && time_in_min + time_duration < g_Simulation_EndTimeInMin + g_Post_Simulation_DurationInMin)
		{
			int time_next_in_interval = g_A2R_simu_interval((time_in_min + time_duration) * 60 / g_number_of_seconds_per_interval);
			int time_in_interval = g_A2R_simu_interval(time_in_min * 60 / g_number_of_seconds_per_interval);
			float total_delay = (m_LinkCumulativeVirtualDelay[time_next_in_interval] - m_LinkCumulativeVirtualDelay[time_in_interval]);
			return total_delay / max(1, get_number_of_vehicles(time_in_min));
		}
		else
		{
			return 0;
		}
	}

};


class CNode
{
public:
	CNode()
	{
		zone_id = 0;
		bOriginNode_ForAgents = false;
		m_OriginNodeSeqNo = -1;
	}

	int node_seq_no;  // sequence number 
	int external_node_id;      //external node number 
	int zone_id;
	double x;
	double y;
	bool bOriginNode_ForAgents;
	int m_OriginNodeSeqNo;
	int sub_flag;
	std::vector<CLink> m_outgoing_node_vector;
	std::vector<CLink> m_sub_outgoing_node_vector;
};
class Cdemand
{
	public:
	Cdemand()
	{
	}
	int o;
	int d;
	float ODvolume;
	
};
std::vector<Cdemand> g_demand_vector;
std::vector<CNode> g_node_vector;
std::vector<CLink> g_link_vector;

class CAgent
{
public:
	unsigned int m_RandomSeed;
	bool m_bGenereated;
	CAgent()
	{
	}
	int demand_type;
	int agent_id;
	int agent_vector_seq_no;
	int agent_service_type;
	int origin_node_id;
	int destination_node_id;
	int time_index;

	float departure_time_in_min;
	float arrival_time_in_min;
	std::vector<int> path_link_seq_no_vector;
	std::vector<float> time_seq_no_vector;

	std::vector<int> path_node_id_vector;

	std::vector<int> likely_node_id_vector;
	std::vector<int> likely_link_id_vector;

	// STS
	std::vector<float> time_seq_vector;

};

vector<CAgent> g_agent_vector;
std::map<int, int> g_map_agent_id_to_agent_vector_seq_no;
class Subnetwork
{
public: Subnetwork() 
	{
	}
		int O_id;
		int D_id;
		vector<int> node_sequence;
		vector<int> link_sequence;

};
int g_number_of_links = 0;
int g_number_of_nodes = 0;
int g_number_of_zones = 0;


int time_max = 0;
int total_arrival = 0;
vector<Subnetwork>g_subnetwork_vector;
void g_ReadInputData()
{
	g_number_of_nodes = 0;
	g_number_of_links = 0;  // initialize  the counter to 0

	int internal_node_seq_no = 0;
	double x, y;
	// step 1: read node file 
	CCSVParser parser;
	if (parser.OpenCSVFile("node.csv", true))
	{
		std::map<int, int> node_id_map;

		while (parser.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{

			string name;

			int node_type;
			int node_id;

			if (parser.GetValueByFieldName("node_id", node_id) == false)
				continue;

			if (g_internal_node_seq_no_map.find(node_id) != g_internal_node_seq_no_map.end())
			{
				continue; //has been defined
			}
			g_internal_node_seq_no_map[node_id] = internal_node_seq_no;
			g_internal_node_seq_no_to_node_id_map[internal_node_seq_no] = node_id;

			parser.GetValueByFieldName("x_coord", x, false);
			parser.GetValueByFieldName("y_coord", y, false);


			CNode node;  // create a node object 

			node.external_node_id = node_id;
			node.node_seq_no = internal_node_seq_no;
			parser.GetValueByFieldName("zone_id", node.zone_id);
			node.x = x;
			node.y = y;
			internal_node_seq_no++;

			g_node_vector.push_back(node);  // push it to the global node vector

			g_number_of_nodes++;
			if (g_number_of_nodes % 1000 == 0)
				cout << "reading " << g_number_of_nodes << " nodes.. " << endl;
		}

		cout << "number of nodes = " << g_number_of_nodes << endl;

		fprintf(g_pFileOutputLog, "number of nodes =,%d\n", g_number_of_nodes);
		parser.CloseCSVFile();
	}
	else
	{
		cout << "node.csv is not opened." << endl;
		g_ProgramStop();
	}

	// step 2: read link file 

	CCSVParser parser_link;

	if (parser_link.OpenCSVFile("road_link.csv", true))
	{
		while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			int from_node_id = 0;
			int to_node_id = 0;
			if (parser_link.GetValueByFieldName("from_node_id", from_node_id) == false)
				continue;
			if (parser_link.GetValueByFieldName("to_node_id", to_node_id) == false)
				continue;

			// add the to node id into the outbound (adjacent) node list

			int internal_from_node_seq_no = g_internal_node_seq_no_map[from_node_id];  // map external node number to internal node seq no. 
			int internal_to_node_seq_no = g_internal_node_seq_no_map[to_node_id];

			CLink link;  // create a link object 

			parser_link.GetValueByFieldName("road_link_id", link.external_link_id);
			parser_link.GetValueByFieldName("link_cost", link.link_cost);

			link.from_node_seq_no = internal_from_node_seq_no;
			link.to_node_seq_no = internal_to_node_seq_no;
			link.link_seq_no = g_number_of_links;
			parser_link.GetValueByFieldName("facility_type", link.type);

			parser_link.GetValueByFieldName("lanes", link.number_of_lanes);

			parser_link.GetValueByFieldName("capacity", link.link_capacity);


			parser_link.GetValueByFieldName("free_speed", link.free_speed);
			parser_link.GetValueByFieldName("length", link.length);

			int travel_time = 60 * link.length / link.free_speed;
			link.travel_time = max(1,travel_time);

			for (int t = 0;t < time_index;t++)
			{
				link.time_dependent_link_cost.push_back(travel_time);
				link.time_dependent_ADMM_visit_counts.push_back(0);
				link.time_dependent_LR_multiplier_vector.push_back(0);
			}
			

			g_node_vector[internal_from_node_seq_no].m_outgoing_node_vector.push_back(link);  // add this link to the corresponding node as part of outgoing node/link
			g_link_vector.push_back(link);
			g_number_of_links++;

			if (g_number_of_links % 1000 == 0)
				cout << "reading " << g_number_of_links << " links.. " << endl;
		}
	}
	else
	{
		cout << "road_link.csv is not opened." << endl;
		g_ProgramStop();
	}


	cout << "number of links = " << g_number_of_links << endl;

	fprintf(g_pFileOutputLog, "number of links =,%d\n", g_number_of_links);

	parser_link.CloseCSVFile();

	/*
	CCSVParser parser_demand;
	if (parser_demand.OpenCSVFile("demand.csv", true))   
	{
		while (parser_demand.ReadRecord())
		{
			Cdemand demand;
			parser_demand.GetValueByFieldName("o", demand.o);
			parser_demand.GetValueByFieldName("d", demand.d);
			parser_demand.GetValueByFieldName("value", demand.ODvolume);
			g_demand_vector.push_back(demand);
		}
		int no_of_agent = 0;
		for (int i = 0;i < g_demand_vector.size();i++)
		{
			for (int a = 0;a < g_demand_vector[i].ODvolume;a++)
			{
				CAgent agent;
				agent.agent_id = no_of_agent;
				agent.origin_node_id = g_demand_vector[i].o;
				agent.destination_node_id = g_demand_vector[i].d;
				agent.departure_time_in_min = 0;
				agent.arrival_time_in_min = time_index;
				g_agent_vector.push_back(agent);
			}
		}
	}
	parser_demand.CloseCSVFile();
	*/
	/*
	g_number_of_agents = 0;
	CCSVParser parser_agent;
	std::vector<int> path_node_sequence;
	string path_node_sequence_str;

	std::vector<int> path_schedule_time_sequence;
	string path_schedule_time_sequence_str;
	
	if (parser_agent.OpenCSVFile("input_agent.csv", true))   // read agent as demand input 
	{
		while (parser_agent.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			CAgent agent;  // create an agent object 
			if (parser_agent.GetValueByFieldName("agent_id", agent.agent_id) == false)
				continue;

			parser_agent.GetValueByFieldName("agent_service_type", agent.agent_service_type);

			int origin_node_id = 0;
			int destination_node_id = 0;
			parser_agent.GetValueByFieldName("from_origin_node_id", origin_node_id, false);

			agent.origin_node_id = origin_node_id;
			parser_agent.GetValueByFieldName("to_destination_node_id", destination_node_id, false);
			agent.destination_node_id = destination_node_id;

			parser_agent.GetValueByFieldName("demand_type", agent.demand_type, false);

			if (agent.demand_type > g_number_of_demand_types)
				g_number_of_demand_types = agent.demand_type;

			ASSERT(g_number_of_demand_types + 1 < _MAX_NUMBER_OF_DEMAND_TYPES);

			parser_agent.GetValueByFieldName("departure_time_in_min", agent.departure_time_in_min);
			parser_agent.GetValueByFieldName("arrival_time_in_min", agent.arrival_time_in_min);

			agent.path_node_id_vector = ParseLineToIntegers(path_node_sequence_str);

			g_agent_vector.push_back(agent);
			g_number_of_agents++;
			if (g_number_of_agents % 1000 == 0)
				cout << "reading = " << g_number_of_agents / 1000 << " k agents..." << endl;

		}
	}

	cout << "number of agents = " << g_agent_vector.size() << endl;

	cout << " Sort agents... " << endl;
	std::sort(g_agent_vector.begin(), g_agent_vector.end());

	cout << " Sorting ends. ..." << endl;

	parser_agent.CloseCSVFile();
	*/
}
std::vector<float> departure_time_vector;
void g_ReadSubInputData()
{
	int no_of_sub = 0;
	CCSVParser parser_subnetwork;
	if (parser_subnetwork.OpenCSVFile("agent.csv", true))
	{
		while (parser_subnetwork.ReadRecord())
		{
			Subnetwork sub;
			parser_subnetwork.GetValueByFieldName("o_zone_id", sub.O_id);
			parser_subnetwork.GetValueByFieldName("d_zone_id", sub.D_id);

			string node_sequence;
			parser_subnetwork.GetValueByFieldName("node_sequence", node_sequence);
			char *node_input = (char *)node_sequence.c_str();
			const char * split_node = ";";
			char *p_node = strtok(node_input, split_node);
			int a;
			while (p_node != NULL)
			{
				sscanf(p_node, "%d", &a);
				sub.node_sequence.push_back(a);
				p_node = strtok(NULL, split_node);
			}

			string link_sequence;
			parser_subnetwork.GetValueByFieldName("link_sequence", link_sequence);
			char *link_input = (char *)link_sequence.c_str();
			const char * split_link = ";";
			char *p_link = strtok(link_input, split_link);
			int b;
			while (p_link != NULL)
			{
				sscanf(p_link, "%d", &b);
				sub.link_sequence.push_back(b);
				p_link = strtok(NULL, split_link);
			}

			g_subnetwork_vector.push_back(sub);
			no_of_sub++;
		}
		cout << "no of sub network = " << no_of_sub << endl;
	}

	else
	{
		cout << "agent.csv is not opened." << endl;
		g_ProgramStop();
	}

}

std::vector<CLink> g_sub_link_vector;
std::vector<CNode> g_sub_node_vector;
void g_read_input_agent()
{
	int no_of_agent = 0;
	int travel_time;

	CCSVParser parser_agent;
	if (parser_agent.OpenCSVFile("agent.csv", true))
	{
		while (parser_agent.ReadRecord())
		{
			int current_no_of_agent;
			parser_agent.GetValueByFieldName("volume", current_no_of_agent);

			std::vector<int> likely_node_sequence;
			std::vector<int> likely_link_sequence;

			string node_sequence;
			parser_agent.GetValueByFieldName("node_sequence", node_sequence);
			char *node_input = (char *)node_sequence.c_str();
			const char * split_node = ";";
			char *p_node = strtok(node_input, split_node);
			int a;
			while (p_node != NULL)
			{
				sscanf(p_node, "%d", &a);
				likely_node_sequence.push_back(a);
				p_node = strtok(NULL, split_node);
			}

			string link_sequence;
			parser_agent.GetValueByFieldName("link_sequence", link_sequence);
			char *link_input = (char *)link_sequence.c_str();
			const char * split_link = ";";
			char *p_link = strtok(link_input, split_link);
			int b;
			while (p_link != NULL)
			{
				sscanf(p_link, "%d", &b);
				likely_link_sequence.push_back(b);
				p_link = strtok(NULL, split_link);
			}
			int departure_time, arrival_time;

			parser_agent.GetValueByFieldName("departure_time", departure_time);
			parser_agent.GetValueByFieldName("arrival_time", arrival_time);

			// space time network time index setting:
			time_index = max(time_index, arrival_time - departure_time);

			for (int i = 0;i < current_no_of_agent;i++)
			{
				CAgent agent;
				parser_agent.GetValueByFieldName("o_zone_id", agent.origin_node_id);
				parser_agent.GetValueByFieldName("d_zone_id", agent.destination_node_id);
				parser_agent.GetValueByFieldName("departure_time", agent.departure_time_in_min);
				parser_agent.GetValueByFieldName("arrival_time", agent.arrival_time_in_min);
				agent.time_index = agent.arrival_time_in_min - agent.departure_time_in_min;
				agent.likely_node_id_vector = likely_node_sequence;
				agent.likely_link_id_vector = likely_link_sequence;
				g_agent_vector.push_back(agent);

				no_of_agent++;
			}


		}
		cout << "no of agent = " << no_of_agent << endl;
		time_index++;
		cout << "time index of this DP:" << time_index << endl;
	}

	else
	{
		cout << "agent.csv is not opened." << endl;
		g_ProgramStop();
	}

}

class NetworkForSP  // mainly for shortest path calculation
{
public:
	int m_threadNo;  // internal thread number 
	std::list<int>  m_SENodeList;  //scan eligible list as part of label correcting algorithm 

	float** m_node_label_cost; // label cost
	int** m_node_time_departure_flag; 
	int** m_vertex_node_predecessor;
	int** m_vertex_time_predecessor;
	int** m_vertex_link_predecessor;

	int* m_node_predecessor;  // predecessor for nodes
	int* m_time_of_node;
	int* m_time_predecessor;
	int* m_node_status_array; // update status 
	float** new_to_node_cost;
	int* m_link_predecessor;  // predecessor for this node points to the previous link that updates its label cost (as part of optimality condition) (for easy referencing)

	FILE* pFileAgentPathLog;  // file output

	float** m_link_volume_array; // link volume for all agents assigned in this network (thread)
	float** m_link_cost_array; // link cost 


	int m_private_origin_seq_no;
	std::vector<int>  m_agent_vector; // assigned agents for computing 
	std::vector<int>  m_node_vector; // assigned nodes for computing 

	NetworkForSP()
	{
		pFileAgentPathLog = NULL;
		m_private_origin_seq_no = -1;

	}

	void AllocateMemory(int number_of_nodes, int number_of_links, int time_index)
	{
		m_node_label_cost = new float*[number_of_nodes];
		m_link_cost_array = new float*[number_of_links];
		new_to_node_cost = new float*[number_of_nodes];
		m_node_time_departure_flag = new int*[number_of_nodes];
		m_vertex_node_predecessor = new int*[number_of_nodes];
		m_vertex_time_predecessor = new int*[number_of_nodes];
		m_vertex_link_predecessor = new int*[number_of_nodes];

		for (int j = 0;j < number_of_nodes;j++)
		{
			m_node_label_cost[j] = new float[time_index];
			new_to_node_cost[j] = new float[time_index];
			m_node_time_departure_flag[j] = new int[time_index];
			m_vertex_node_predecessor[j] = new int[time_index];
			m_vertex_time_predecessor[j] = new int[time_index];
			m_vertex_link_predecessor[j] = new int[time_index];

		}
	}

	~NetworkForSP()
	{

		if (m_node_label_cost != NULL)
			delete m_node_label_cost;

		if (m_node_predecessor != NULL)
			delete m_node_predecessor;

		if (m_node_status_array != NULL)
			delete m_node_status_array;

		if (m_link_predecessor != NULL)
			delete m_link_predecessor;

		if (m_link_volume_array != NULL)
			delete m_link_volume_array;

		if (m_link_cost_array != NULL)
			delete m_link_cost_array;

		if (pFileAgentPathLog != NULL)
			fclose(pFileAgentPathLog);
	}

	int optimal_label_correcting(int origin_node, int destination_node, int departure_time, int arrival_time, int agent_id)
		// time-dependent label correcting algorithm with double queue implementation
	{
		int internal_debug_flag = 0;

		for (int i = 0; i < g_number_of_nodes; i++) //Initialization for all nodes
		{
			//by slt
			for (int t = 0;t < time_index;t++)
			{
				m_node_label_cost[i][t] = _MAX_LABEL_COST;
				m_node_time_departure_flag[i][t] = -1;
				m_vertex_time_predecessor[i][t] = -1;
				m_vertex_node_predecessor[i][t] = -1;
				m_vertex_link_predecessor[i][t] = -1;
			}
		}

		//Initialization for origin node at the preferred departure time, at departure time, cost = 0, otherwise, the g_A2R_simu_interval at origin node

		m_node_label_cost[origin_node][departure_time] = 0;

		/*
		vector<int> vertex; //(i,t)
		vector<vector<int>> SElist;//((i,t),(i,t)...)
		vertex.push_back(origin_node);
		vertex.push_back(departure_time);
		SElist.push_back(vertex);
		while (SElist.size() != 0)
		{
			vector<int> from_vertex = SElist.front();
			int from_node = from_vertex.front();
			int departure_time = from_vertex.back();
			SElist.erase(SElist.begin());
			//SElist.pop_front;
			for (int i = 0;i < g_node_vector[from_node].m_outgoing_node_vector.size();i++)
			{
				int to_node = g_node_vector[from_node].m_outgoing_node_vector[i].to_node_seq_no;
				int travel_time = g_node_vector[from_node].m_outgoing_node_vector[i].travel_time;
				int next_time = departure_time + travel_time;
				int current_link = g_node_vector[from_node].m_outgoing_node_vector[i].link_seq_no;

				if (next_time <= arrival_time)// cannot exceed final ST vertex
				{
					float current_cost = g_link_vector[current_link].time_dependent_link_cost[departure_time];
					float new_cost = m_node_label_cost[from_node][departure_time] + current_cost;
					if (new_cost < m_node_label_cost[to_node][next_time])
					{
						m_node_label_cost[to_node][next_time] = new_cost;
						m_vertex_node_predecessor[to_node][next_time] = from_node;
						m_vertex_time_predecessor[to_node][next_time] = departure_time;
						m_vertex_link_predecessor[to_node][next_time] = current_link;

						vector<int> new_vertex;
						new_vertex.push_back(to_node);
						new_vertex.push_back(next_time);
						SElist.push_back(new_vertex);
					}
				}

			}
		}
		*/
		

		for (int i = 0; i < g_number_of_nodes; i++) //Initialization for all nodes
		{
										 //by slt
			for (int t = 0;t < time_index;t++)
			{
				m_node_label_cost[i][t] = _MAX_LABEL_COST;
				m_node_time_departure_flag[i][t] = -1;
				m_vertex_time_predecessor[i][t] = -1;
				m_vertex_node_predecessor[i][t] = -1;
				m_vertex_link_predecessor[i][t] = -1;
			}
		}


		//Initialization for origin node at the preferred departure time, at departure time, cost = 0, otherwise, the g_A2R_simu_interval at origin node
		
		m_node_label_cost[origin_node][departure_time] = 0;
		// initialization
		m_node_time_departure_flag[origin_node][departure_time] = 1;//initialization
		for (int t = departure_time;t < arrival_time;t++)
		{
			for (int i = 0;i < g_sub_node_vector.size();i++)
			{
				int current_node = g_sub_node_vector[i].node_seq_no;
				if (m_node_time_departure_flag[current_node][t] == 1)// for each available space-time vertex
				{
					for (int j = 0; j < g_sub_node_vector[i].m_sub_outgoing_node_vector.size();j++)
					{
						int current_travel_time = g_sub_node_vector[i].m_sub_outgoing_node_vector[j].travel_time;
						int next_time = t + current_travel_time;
						if (next_time <= arrival_time)// cannot exceed final ST vertex
						{
							int current_link = g_sub_node_vector[i].m_sub_outgoing_node_vector[j].link_seq_no;
							int current_link_no;
							for (int u = 0;u < g_sub_link_vector.size();u++)
							{
								if (g_sub_link_vector[u].link_seq_no == current_link)
								{
									current_link_no = u;
								}
							}
							float current_cost = g_sub_link_vector[current_link_no].time_dependent_link_cost[t];
							int to_node = g_sub_node_vector[i].m_sub_outgoing_node_vector[j].to_node_seq_no;
							float new_cost = m_node_label_cost[current_node][t] + current_cost;
							if (new_cost < m_node_label_cost[to_node][next_time])
							{
								m_node_label_cost[to_node][next_time] = new_cost;
								m_vertex_node_predecessor[to_node][next_time] = current_node;
								m_vertex_time_predecessor[to_node][next_time] = t;
								m_vertex_link_predecessor[to_node][next_time] = current_link_no;
								m_node_time_departure_flag[to_node][next_time] = 1;
								//cout << "next_time:" << next_time << "to_node:" << to_node << "agent:" << agent_id << endl;
								//cout << "time_pre:" << t << endl;
								//cout << "link_pre:" << current_link_no << endl;
							}
						}
					}
				}
			}
		}
		

		CAgent* p_agent = &(g_agent_vector[agent_id]);
		p_agent->path_link_seq_no_vector.clear();  // reset;
		p_agent->path_node_id_vector.clear();  // reset;
		p_agent->time_seq_no_vector.clear();

		int t = arrival_time;
		int i = destination_node;
		p_agent->time_seq_no_vector.push_back(t);
		p_agent->path_node_id_vector.push_back(i);
		while (m_vertex_time_predecessor[i][t] != departure_time)
		{
			//cout << m_vertex_time_predecessor[i][t] <<  m_vertex_node_predecessor[i][t]  << endl;
			int l = m_vertex_link_predecessor[i][t];
			p_agent->path_link_seq_no_vector.push_back(l);

			int node = m_vertex_node_predecessor[i][t];
			p_agent->path_node_id_vector.push_back(node);

			int time = m_vertex_time_predecessor[i][t];
			p_agent->time_seq_no_vector.push_back(time);

			t = time;
			i = node;
		}

		p_agent->path_node_id_vector.push_back(origin_node);
		p_agent->time_seq_no_vector.push_back(departure_time);
		p_agent->path_link_seq_no_vector.push_back(m_vertex_link_predecessor[i][t]);
				
				
		std::reverse(std::begin(p_agent->path_node_id_vector),
			std::end(p_agent->path_node_id_vector));
		std::reverse(std::begin(p_agent->path_link_seq_no_vector),
			std::end(p_agent->path_link_seq_no_vector));
		std::reverse(std::begin(p_agent->time_seq_no_vector),
			std::end(p_agent->time_seq_no_vector));
		//update mu:
		for (int i = 0;i < p_agent->time_seq_no_vector.size()-1;i++)
		{
			g_link_vector[p_agent->path_link_seq_no_vector[i]].time_dependent_ADMM_visit_counts[p_agent->time_seq_no_vector[i]]++;
		}


		if (destination_node == -1)
			return 1;  // one to all shortest path
		else
			return -1;


	}

};



int g_number_of_CPU_threads()
{
	int number_of_threads = omp_get_max_threads();

	int max_number_of_threads = 8;

	if (number_of_threads > max_number_of_threads)
		number_of_threads = max_number_of_threads;

	return number_of_threads;

}

NetworkForSP* pNetworkForSP = NULL;

int g_state_to_load_mapping[_MAX_STATES];
int g_initial_state_no = 0;   // customized 

int g_number_of_LR_iterations = 20;
int g_number_of_ADMM_iterations = 100;
int g_CurrentLRIterationNumber = 0;
int g_Number_Of_Iterations_With_Memory = 5;

float g_best_upper_bound = 99999;
float g_best_lower_bound = -99999;
float g_stepSize = 0;
float g_penalty_PHO = 0.1;
float g_minimum_subgradient_step_size = 0.1;

NetworkForSP* g_pNetworkForSP = NULL;
int LR_iteration_no = 10;

void reset_subnetwork(int origin_node,int destination_node)
{
	for (int i = 0;i < g_link_vector.size();i++)
	{
		g_link_vector[i].sub_flag = 0;
	}
	g_sub_link_vector.clear();
	for (int i = 0;i < g_node_vector.size();i++)
	{
		g_node_vector[i].sub_flag = 0;
		g_node_vector[i].m_sub_outgoing_node_vector.clear();
	}
	for (int i = 0;i < g_sub_node_vector.size();i++)
	{
		g_sub_node_vector[i].m_sub_outgoing_node_vector.clear();
	}
	g_sub_node_vector.clear();
	//sub_network setting
	//locate node and link sequence of current OD; set sub_flag = 1.
	for (int i = 0;i < g_subnetwork_vector.size();i++)
	{
		if (g_subnetwork_vector[i].O_id == origin_node && g_subnetwork_vector[i].D_id == destination_node)
		{
			for (int j = 0;j < g_subnetwork_vector[i].node_sequence.size();j++)
			{
				int node = g_subnetwork_vector[i].node_sequence[j] - 1;
				g_node_vector[node].sub_flag = 1;
			}
			for (int j = 0;j < g_subnetwork_vector[i].link_sequence.size();j++)
			{
				int link = g_subnetwork_vector[i].link_sequence[j] - 1;
				g_link_vector[link].sub_flag = 1;
			}
		}
	}
	//put sub_link into sub_link_vector including waiting link(at origin and destination).
	for (int i = 0;i < g_link_vector.size();i++)
	{
		if (g_link_vector[i].sub_flag == 1)
		{
			g_sub_link_vector.push_back(g_link_vector[i]);
		}
		if (g_link_vector[i].from_node_seq_no == origin_node - 1 && g_link_vector[i].to_node_seq_no == origin_node - 1)
		{
			g_link_vector[i].sub_flag = 1;
			g_sub_link_vector.push_back(g_link_vector[i]);
		}
		if (g_link_vector[i].from_node_seq_no == destination_node - 1 && g_link_vector[i].to_node_seq_no == destination_node - 1)
		{
			g_link_vector[i].sub_flag = 1;
			g_sub_link_vector.push_back(g_link_vector[i]);
		}
	}
	// set nodes' outgoing link and put them into sub_outgoing_vector.
	for (int i = 0;i < g_node_vector.size();i++)
	{
		if (g_node_vector[i].sub_flag == 1)
		{
			int node_seq = g_node_vector[i].node_seq_no;
			for (int j = 0;j < g_sub_link_vector.size();j++)
			{
				if (node_seq == g_sub_link_vector[j].from_node_seq_no && g_sub_link_vector[j].sub_flag == 1)
				{
					g_node_vector[i].m_sub_outgoing_node_vector.push_back(g_sub_link_vector[j]);
				}
			}
			g_sub_node_vector.push_back(g_node_vector[i]);
		}
	}

}

void find_sp_for_each_agent_use_LR()
{
	g_pNetworkForSP = new NetworkForSP; // create n copies of network, each for a subset of agents to use	
	g_pNetworkForSP->AllocateMemory(g_number_of_nodes, g_number_of_links,time_index);
	for (int a = 0; a < g_agent_vector.size(); a++)
	{
		for (int l = 0; l < g_number_of_links; l++)
		{
			for (int t = 0;t < time_index;t++)
			{
				g_link_vector[l].time_dependent_link_cost[t] +=  g_link_vector[l].time_dependent_LR_multiplier_vector[t];
			}
		}
		if (a == 0)
		{
			reset_subnetwork(g_agent_vector[a].origin_node_id,g_agent_vector[a].destination_node_id);
			g_pNetworkForSP->optimal_label_correcting(g_agent_vector[a].origin_node_id-1 , g_agent_vector[a].destination_node_id-1, g_agent_vector[a].departure_time_in_min, g_agent_vector[a].arrival_time_in_min, a);
		}
		else if (a >= 1 && g_agent_vector[a].origin_node_id == g_agent_vector[a - 1].origin_node_id&& g_agent_vector[a].destination_node_id == g_agent_vector[a - 1].destination_node_id)
		{
			g_pNetworkForSP->optimal_label_correcting(g_agent_vector[a].origin_node_id-1, g_agent_vector[a].destination_node_id-1, g_agent_vector[a].departure_time_in_min, g_agent_vector[a].arrival_time_in_min, a);
		}
		else
		{
			reset_subnetwork(g_agent_vector[a].origin_node_id, g_agent_vector[a].destination_node_id);
			//cout << a << " agents..." << "arrival time:" << g_agent_vector[a].arrival_time_in_min << endl;
			g_pNetworkForSP->optimal_label_correcting(g_agent_vector[a].origin_node_id-1, g_agent_vector[a].destination_node_id-1, g_agent_vector[a].departure_time_in_min, g_agent_vector[a].arrival_time_in_min, a);
		}
		if (a % 10000 == 0)
		{
			cout << a << " agents..." << endl;
		}
	}
}
clock_t start_t, end_t, total_t;

void LR_process()
{
	float stepsize;

	for (int k = 1; k < LR_iteration_no; k++)
	{
		//LR process
		//initialization
		cout << "iteration" << k << "..." << endl;		
		if (k == 1)
		{
			find_sp_for_each_agent_use_LR();
		}
		stepsize = 1.0 / (k+1);
		//update lambda
		cout << "update lambda..." << endl;
		for (int l = 0;l < g_link_vector.size();l++)
		{
			for (int t = 0;t < time_index;t++)
			{
				int current_v = 0;
				/*for (int a = 0;a < g_agent_vector.size();a++)
				{
					if (t < g_agent_vector[a].path_link_seq_no_vector.size())
					{
						if (g_agent_vector[a].path_link_seq_no_vector[t] == l)
						{
							current_v++;
						}
					}
				}*/
				current_v = g_link_vector[l].time_dependent_ADMM_visit_counts[t];
				g_link_vector[l].time_dependent_LR_multiplier_vector[t] = max(0,g_link_vector[l].time_dependent_LR_multiplier_vector[t] + (current_v - g_link_vector[l].link_capacity)*g_penalty_RHO);
			}
		}
		//ADMM process
		//set mu
		g_pNetworkForSP = new NetworkForSP; // create n copies of network, each for a subset of agents to use	
		g_pNetworkForSP->AllocateMemory(g_number_of_nodes, g_number_of_links, time_index);
		cout << "update mu..." << endl;
		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			//set mu; all the other = all - self;
			/*
			for (int i = 0;i < g_agent_vector[a].time_seq_no_vector.size() - 1;i++)
			{
				g_link_vector[g_agent_vector[a].path_link_seq_no_vector[i]].time_dependent_ADMM_visit_counts[g_agent_vector[a].time_seq_no_vector[i]]--;
			}
			for (int l = 0;l < g_link_vector.size()-g_node_vector.size();l++)//: not updating waiting cost
			{
				for (int t = 0;t <time_index;t++)
				{
					int mu;
					for (int a1 = 0;a1 < g_agent_vector.size();a1++)
					{
						if (a != a1)
						{
							if (t < g_agent_vector[a1].path_link_seq_no_vector.size())
							{
								if (g_agent_vector[a1].path_link_seq_no_vector[t] == l)
									mu++;
							}
						}
					}
					
					mu = g_link_vector[l].time_dependent_ADMM_visit_counts[t];
					 g_link_vector[l].time_dependent_link_cost[t] = g_link_vector[l].time_dependent_link_cost[t]+ g_link_vector[l].time_dependent_LR_multiplier_vector[t] + g_penalty_RHO / 2 + g_penalty_RHO*(mu - g_link_vector[l].link_capacity);

				}
			}*/
			if (a == 0)
			{
				for (int i = 0;i < g_agent_vector[a].time_seq_no_vector.size() - 1;i++)
				{
					g_link_vector[g_agent_vector[a].path_link_seq_no_vector[i]].time_dependent_ADMM_visit_counts[g_agent_vector[a].time_seq_no_vector[i]]--;
				}


				reset_subnetwork(g_agent_vector[a].origin_node_id, g_agent_vector[a].destination_node_id);
				for (int l = 0;l < g_sub_link_vector.size();l++)//: not updating waiting cost
				{
					for (int t = 0;t < g_agent_vector[a].time_index;t++)
					{
						int mu;
						mu = g_sub_link_vector[l].time_dependent_ADMM_visit_counts[t];
						g_sub_link_vector[l].time_dependent_link_cost[t] = g_sub_link_vector[l].time_dependent_link_cost[t] + g_sub_link_vector[l].time_dependent_LR_multiplier_vector[t] + g_penalty_RHO / 2 + g_penalty_RHO * (mu - g_sub_link_vector[l].link_capacity);

					}
				}
				g_pNetworkForSP->optimal_label_correcting(g_agent_vector[a].origin_node_id-1, g_agent_vector[a].destination_node_id-1, g_agent_vector[a].departure_time_in_min, g_agent_vector[a].arrival_time_in_min, a);
			}
			else if (a >= 1 && g_agent_vector[a].origin_node_id == g_agent_vector[a - 1].origin_node_id&& g_agent_vector[a].destination_node_id == g_agent_vector[a - 1].destination_node_id)
			{
				for (int i = 0;i < g_agent_vector[a].time_seq_no_vector.size() - 1;i++)
				{
					g_link_vector[g_agent_vector[a].path_link_seq_no_vector[i]].time_dependent_ADMM_visit_counts[g_agent_vector[a].time_seq_no_vector[i]]--;
				}
				for (int l = 0;l < g_sub_link_vector.size();l++)//: update link cost
				{
					for (int t = 0;t < g_agent_vector[a].time_index;t++)
					{
						int mu;
						mu = g_sub_link_vector[l].time_dependent_ADMM_visit_counts[t];
						g_sub_link_vector[l].time_dependent_link_cost[t] = g_sub_link_vector[l].time_dependent_link_cost[t] + g_sub_link_vector[l].time_dependent_LR_multiplier_vector[t] + g_penalty_RHO / 2 + g_penalty_RHO * (mu - g_sub_link_vector[l].link_capacity);

					}
				}
				g_pNetworkForSP->optimal_label_correcting(g_agent_vector[a].origin_node_id-1, g_agent_vector[a].destination_node_id-1, g_agent_vector[a].departure_time_in_min, g_agent_vector[a].arrival_time_in_min, a);
			}
			else
			{
				for (int i = 0;i < g_agent_vector[a].time_seq_no_vector.size() - 1;i++)
				{
					g_link_vector[g_agent_vector[a].path_link_seq_no_vector[i]].time_dependent_ADMM_visit_counts[g_agent_vector[a].time_seq_no_vector[i]]--;
				}
				reset_subnetwork(g_agent_vector[a].origin_node_id, g_agent_vector[a].destination_node_id);
				for (int l = 0;l < g_sub_link_vector.size();l++)//: not updating waiting cost
				{
					for (int t = 0;t < g_agent_vector[a].time_index;t++)
					{
						int mu;
						mu = g_sub_link_vector[l].time_dependent_ADMM_visit_counts[t];
						g_sub_link_vector[l].time_dependent_link_cost[t] = g_sub_link_vector[l].time_dependent_link_cost[t] + g_sub_link_vector[l].time_dependent_LR_multiplier_vector[t] + g_penalty_RHO / 2 + g_penalty_RHO * (mu - g_sub_link_vector[l].link_capacity);

					}
				}
				g_pNetworkForSP->optimal_label_correcting(g_agent_vector[a].origin_node_id-1, g_agent_vector[a].destination_node_id-1, g_agent_vector[a].departure_time_in_min, g_agent_vector[a].arrival_time_in_min, a);
			}
			if (a % 10000 == 0)
			{
				cout << a << " agents..." << endl;
			}

		}
		end_t = clock();

		total_t = (end_t - start_t);

		cout << "CPU Running Time = " << total_t << " milliseconds for iteration" << k << endl;

	}


}

void g_OutputFiles_1()
{
	cout << "outputing files... " << endl;

	//output_agent.csv
	FILE* g_pFileAgent = NULL;
	g_pFileAgent = fopen("output_agent.csv", "w");
	if (g_pFileAgent == NULL)
	{
		cout << "File output_agent.csv cannot be opened." << endl;
		g_ProgramStop();
	}
	else
	{

		fprintf(g_pFileAgent, "agent_id,agent_type,o_node_id,d_node_id,node_sequence,link_sequence,time_sequence\n");

		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			CAgent* p_agent = &(g_agent_vector[a]);

			fprintf(g_pFileAgent, "%d,%d,%d,%d,",
				p_agent->agent_id,
				p_agent->agent_service_type,
				p_agent->origin_node_id,
				p_agent->destination_node_id
			);

			// path node id sequence
			for (int i = 0; i < p_agent->path_node_id_vector.size(); i++)
			{
				int internal_node_id = p_agent->path_node_id_vector[i];
				int external_node_id = internal_node_id + 1;
				if (i == p_agent->path_node_id_vector.size()-1)
				{
					fprintf(g_pFileAgent, "%d", external_node_id);
				}
				else
				{
					fprintf(g_pFileAgent, "%d;", external_node_id);
				}
			}

			fprintf(g_pFileAgent, ",");
			// path link id sequence 
			for (int i = 0; i < p_agent->path_link_seq_no_vector.size(); i++)
			{
				int internal_link_id = p_agent->path_link_seq_no_vector[i];
				int external_link_id = g_link_vector[internal_link_id].external_link_id;
				if (i == p_agent->path_link_seq_no_vector.size() - 1)
				{
					fprintf(g_pFileAgent, "%d", external_link_id);
				}
				else
				{
					fprintf(g_pFileAgent, "%d;", external_link_id);
				}
			}
			fprintf(g_pFileAgent, ",");
			// time sequence 
			for (int i = 0; i < p_agent->time_seq_no_vector.size(); i++)
			{
				int internal_time = p_agent->time_seq_no_vector[i];
				//time transform
				internal_time = internal_time + 1200;
				if (i == p_agent->time_seq_no_vector.size()-1)
				{
					fprintf(g_pFileAgent, "%d", internal_time);
				}
				else
				{
					fprintf(g_pFileAgent, "%d;", internal_time);
				}
			}
			fprintf(g_pFileAgent, ",");

			fprintf(g_pFileAgent, "\n");

		}
		//
		fclose(g_pFileAgent);
	}

}

int main(int argc, TCHAR* argv[], TCHAR* envp[])
{
	g_pFileDebugLog = fopen("Debug.txt", "w");
	if (g_pFileDebugLog == NULL)
	{
		cout << "File Debug.txt cannot be opened." << endl;
		g_ProgramStop();
	}
	g_pFileOutputLog = fopen("output_solution.csv", "w");
	if (g_pFileOutputLog == NULL)
	{
		cout << "File output_solution.csv cannot be opened." << endl;
		g_ProgramStop();
	}
	g_read_input_agent();

	g_ReadSubInputData();
	g_ReadInputData();

	LR_process();
	g_OutputFiles_1();

	cout << "it is done!" << endl;

	return 1;
	
}

