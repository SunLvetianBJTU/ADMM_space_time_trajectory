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
#define _MAX_LABEL_COST 99999
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

int time_index = 11;

double g_number_of_seconds_per_interval = 6;  // 0.2 seconds for 300 intervals per min
int g_number_of_simulation_intervals = 600 * 60 / g_number_of_seconds_per_interval;    // 60min
int g_number_of_optimization_time_intervals = 60;

int g_Simulation_StartTimeInMin = 9999;
int g_Simulation_EndTimeInMin = 0;
int g_start_simu_interval_no, g_end_simu_interval_no;

int g_Post_Simulation_DurationInMin = 120;
int g_dp_algorithm_debug_flag = 0;
float g_penalty_RHO = 2;
int g_optimization_method = 2;
float LR_multiplier;
std::map<int, int> g_link_key_to_seq_no_map;  // hush table, map key to internal link sequence no. 


//mfd
int g_TAU;
std::map<int, int> g_internal_node_seq_no_map;  // hush table, map external node number to internal node sequence no. 
std::map<int, int> g_internal_node_seq_no_to_node_id_map;  // hush table, map external node number to internal node sequence no. 

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

		m_LinkOutFlowCapacity = NULL;
		m_LinkInFlowCapacity = NULL;
		m_LinkCumulativeArrival = NULL;
		m_LinkCumulativeDeparture = NULL;
		m_LinkCumulativeVirtualDelay = NULL;

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
	int VRP_time_window_begin, VRP_time_window_end;
	float base_price;
	std::vector<float> travel_time_vector;
	std::vector<float> time_dependent_LR_multiplier_vector, time_dependent_external_cost_vector, time_dependent_ADMM_multiplier_vector, time_dependent_discharge_rate, time_dependent_inflow_rate;
	std::vector<float> time_dependent_LR_multiplier_vector_for_ADMM;

	std::vector<int> time_dependent_visit_counts, time_dependent_ADMM_visit_counts,
		time_depedent_capacity_vector;
	std::vector<float> time_dependent_link_cost;
	std::vector<float> time_dependent_link_cost_for_LR, time_dependent_link_cost_for_ADMM;

	std::vector<float> time_dependent_travel_time_vector;

	//ADMM_multiplier_matrix is used in the searching process and updated by LR_multiplier_matrix
	int max_allowed_waiting_time;

	void Setup_State_Dependent_Data_Matrix(int number_of_optimization_time_intervals)
	{

		for (int t = 0; t < number_of_optimization_time_intervals; t++)
		{
			time_dependent_visit_counts.push_back(0);
			time_dependent_ADMM_visit_counts.push_back(0);
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
};

class CNode
{
public:
	CNode()
	{
		zone_id = 0;
		accessible_node_count = 0;
		bOriginNode_ForAgents = false;
		m_OriginNodeSeqNo = -1;
	}

	int accessible_node_count;

	int node_seq_no;  // sequence number 
	int external_node_id;      //external node number 
	int zone_id;
	float departure_time;
	double x;
	double y;

	int waiting_flag;
	int external_travel_time;
	float waiting_cost;

	bool bOriginNode_ForAgents;
	int m_OriginNodeSeqNo;

	std::vector<CLink> m_outgoing_node_vector;

};

std::vector<CNode> g_node_vector;
std::vector<CLink> g_link_vector;
class CCost
{
public:
	std::vector<float> time_dependent_link_cost_for_ADMM;
	std::vector<float> time_dependent_link_cost_for_LR;
};

class CAgent
{
public:
	unsigned int m_RandomSeed;
	bool m_bGenereated;
	CAgent()
	{
		agent_vector_seq_no = -1;
		agent_service_type = 0;  //0: pax vehicle 1: travler 2: scheduled transportation vehicle
		m_bMoveable = 1;
		fixed_path_flag = 0;
		vehicle_seat_capacity = 1;
		m_bGenereated = false;
		PCE_factor = 1.0;
		path_cost = 0;
		m_Veh_LinkArrivalTime_in_simu_interval = NULL;
		m_Veh_LinkDepartureTime_in_simu_interval = NULL;
		m_bCompleteTrip = false;
		departure_time_in_min = 0;
		// vrp

		transportation_time_cost = 0;
		schedule_early_cost = 0;
		schedule_delay_cost = 0;

		earliest_departure_time = 0;
		departure_time_window = 1;

		latest_arrival_time = 0;
		arrival_time_window = 1;

	}
	int fixed_path_flag;
	int demand_type;
	int agent_id;
	int agent_vector_seq_no;
	int agent_service_type;
	int m_bMoveable;
	int origin_node_id;
	int destination_node_id;

	int origin_zone_seq_no;
	int destination_zone_seq_no;

	float departure_time_in_min;
	int departure_time_in_simu_interval;
	float arrival_time_in_min;
	float PCE_factor;  // passenger car equivalent : bus = 3
	float path_cost;
	std::vector<int> path_link_seq_no_vector;
	std::vector<int> path_timestamp_vector;

	std::vector<int> path_link_seq_no_vector_for_LR;
	std::vector<float> time_seq_no_vector_for_LR;

	std::vector<int> path_link_seq_no_vector_for_ADMM;
	std::vector<float> time_seq_no_vector_for_ADMM;

	int m_path_link_seq_no_vector_size;

	std::vector<int> path_node_id_vector;
	std::vector<int> path_schedule_time_vector;

	std::vector<int> path_node_id_vector_for_LR;
	std::vector<int> path_node_id_vector_for_ADMM;

	std::vector<CCost> time_dependent_link_cost_for_LR;
	std::vector<CCost> time_dependent_link_cost_for_ADMM;

	int m_current_link_seq_no;
	int* m_Veh_LinkArrivalTime_in_simu_interval;
	int* m_Veh_LinkDepartureTime_in_simu_interval;

	int vehicle_seat_capacity;
	int VRP_group_id;

	std::list<int>  m_PassengerList;

	vector<int> set_of_allowed_links;
	vector<int> m_set_of_allowed_links_flag;
	vector<int> set_of_allowed_nodes;
	vector<int> set_of_allowed_links_LR;
	vector<int> m_set_of_allowed_links_flag_LR;

	// STS
	float transportation_time_cost;
	float schedule_early_cost;
	float schedule_delay_cost;

	float earliest_departure_time;
	int departure_time_in_simulation_interval;
	float departure_time_window;

	float latest_arrival_time;
	float arrival_time_window;
	std::vector<float> time_seq_vector;

	// STS
	void Pickup(int p)
	{
		if (m_PassengerList.size() < vehicle_seat_capacity)
		{
			m_PassengerList.push_back(p);
		}

	}


	int GetRemainingCapacity()
	{

		return vehicle_seat_capacity - m_PassengerList.size();
	}

	//above are simulated 

	bool m_bCompleteTrip;
	bool operator<(const CAgent &other) const
	{
		return departure_time_in_min < other.departure_time_in_min;
	}

	std::map<int, int> m_VRP_ADMM_link_time_map;

};

class CFlow
{
public:
	CFlow()
	{};
	int time_index;
	int link_no;
	float cumulative_arrival;
	float outflow_rate;
	float cumulative_departure;
	float time_dependent_inflow_rate;
};
vector<CFlow> g_flow_vector;
vector<CAgent> g_agent_vector;

std::map<int, int> g_map_agent_id_to_agent_vector_seq_no;

int g_number_of_links = 0;
int g_number_of_nodes = 0;
int g_number_of_zones = 0;
int time_max = 0;
int total_arrival = 0;
void g_ReadInputData()
{
	g_number_of_nodes = 0;
	g_number_of_links = 0;  // initialize  the counter to 0

	int internal_node_seq_no = 0;
	double x, y;
	// step 1: read node file 
	CCSVParser parser;
	if (parser.OpenCSVFile("input_node.csv", true))
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

			parser.GetValueByFieldName("x", x, false);
			parser.GetValueByFieldName("y", y, false);


			CNode node;  // create a node object 

			node.external_node_id = node_id;
			node.node_seq_no = internal_node_seq_no;
			parser.GetValueByFieldName("zone_id", node.zone_id);
			node.departure_time = 0;
			node.x = x;
			node.y = y;
			internal_node_seq_no++;
			parser.GetValueByFieldName("waiting_flag", node.waiting_flag);
			if (node.waiting_flag == 1)
			{
				parser.GetValueByFieldName("external_travel_time", node.external_travel_time);
				parser.GetValueByFieldName("waiting_cost", node.waiting_cost);
			}

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
		cout << "input_node.csv is not opened." << endl;
		g_ProgramStop();
	}

	// step 2: read link file 

	CCSVParser parser_link;

	if (parser_link.OpenCSVFile("input_link.csv", true))
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

			parser_link.GetValueByFieldName("link_id", link.external_link_id);
			parser_link.GetValueByFieldName("link_cost", link.link_cost);

			link.from_node_seq_no = internal_from_node_seq_no;
			link.to_node_seq_no = internal_to_node_seq_no;
			link.link_seq_no = g_number_of_links;
			//link.to_node_seq_no = internal_to_node_seq_no;
			link.LR_multiplier = 0;
			float m_End_of_RedTime_in_min = 0;
			parser_link.GetValueByFieldName("m_End_of_RedTime_in_min", m_End_of_RedTime_in_min);
			link.m_End_of_RedTime_in_simu_interval = m_End_of_RedTime_in_min * 60 / g_number_of_seconds_per_interval;

			parser_link.GetValueByFieldName("link_flag", link.link_flag);
			parser_link.GetValueByFieldName("link_type", link.type);
			parser_link.GetValueByFieldName("service_type", link.service_type, false);

			if (link.service_type != 0)
			{
				parser_link.GetValueByFieldName("VRP_load_id", link.VRP_load_id, false);
				parser_link.GetValueByFieldName("VRP_group_id", link.VRP_group_id, false);
				parser_link.GetValueByFieldName("VRP_time_window_begin", link.VRP_time_window_begin, false);
				parser_link.GetValueByFieldName("VRP_time_window_end", link.VRP_time_window_end, false);
				parser_link.GetValueByFieldName("VRP_load_difference", link.VRP_load_difference, false);
			}
			float length = 1; // km or mile
			float speed_limit = 1;
			parser_link.GetValueByFieldName("length", length);
			parser_link.GetValueByFieldName("speed_limit", speed_limit);
			parser_link.GetValueByFieldName("base_price", link.base_price);
			parser_link.GetValueByFieldName("BPR_alpha_term", link.BRP_alpha);
			parser_link.GetValueByFieldName("BPR_beta_term", link.BRP_beta);

			int number_of_lanes = 1;
			float lane_cap = 1000;
			parser_link.GetValueByFieldName("number_of_lanes", link.number_of_lanes);
			parser_link.GetValueByFieldName("lane_cap", lane_cap);

			parser_link.GetValueByFieldName("link_cap", link.link_capacity);

			//link.link_capacity = lane_cap* number_of_lanes;

			string demand_type_code;
			for (int d = 0; d < _MAX_NUMBER_OF_DEMAND_TYPES; d++)
			{
				link.demand_type_code[d] = true;
				link.demand_type_TTcost[d] = 0;
			}

			parser_link.GetValueByFieldName("demand_type_code", demand_type_code);

			if (demand_type_code.size() > 0)  //demand type code has a string
			{
				for (int d = 0; d < _MAX_NUMBER_OF_DEMAND_TYPES; d++)
				{
					link.demand_type_code[d] = false;
					CString demand_type_number;
					demand_type_number.Format(_T("%d"), d);

					std::string str_number = CString2StdString(demand_type_number);

					if (demand_type_code.find(str_number) != std::string::npos)   // find this number
					{
						link.demand_type_code[d] = true;  // allow this demand type
					}
				}
			}


			for (int d = 0; d < _MAX_NUMBER_OF_DEMAND_TYPES; d++)
			{
				CString demand_type_number;
				demand_type_number.Format(_T("demand_type_%d_TTcost"), d);

				std::string str_number = CString2StdString(demand_type_number);
				parser_link.GetValueByFieldName(str_number, link.demand_type_TTcost[d]);

			}
			for (int t = 0;t < time_index;t++)
			{
				float link_cost;
				parser_link.GetValueByFieldName("link_cost", link_cost);
				link.time_dependent_link_cost_for_LR.push_back(link_cost);
				link.time_dependent_link_cost_for_ADMM.push_back(link_cost);

				int link_travel_time;
				parser_link.GetValueByFieldName("travel_time", link_travel_time);
				link.time_dependent_travel_time_vector.push_back(link_travel_time);

				float LR_multiplier = 0;
				link.time_dependent_LR_multiplier_vector.push_back(LR_multiplier);
				link.time_dependent_LR_multiplier_vector_for_ADMM.push_back(LR_multiplier);
				link.time_dependent_discharge_rate.push_back(link.link_capacity);
			}

			link.free_flow_travel_time_in_min = length / speed_limit * 60;

			float external_travel_time = -1;
			parser_link.GetValueByFieldName("external_travel_time", external_travel_time, false);

			if (external_travel_time >= 0.1)
			{  // reset 
				link.free_flow_travel_time_in_min = external_travel_time;
			}

			link.length = length;
			link.cost = length / speed_limit * 60; // min // calculate link cost based length and speed limit // later we should also read link_capacity, calculate g_A2R_simu_interval 

			g_node_vector[internal_from_node_seq_no].m_outgoing_node_vector.push_back(link);  // add this link to the corresponding node as part of outgoing node/link

			long link_key = internal_from_node_seq_no * _MAX_NUMBER_OF_PHYSICAL_NODES + internal_to_node_seq_no;

			g_link_key_to_seq_no_map[link_key] = link.link_seq_no;

			link.CalculateBPRFunctionAndCost(); // initial link travel time value
			g_link_vector.push_back(link);
			g_number_of_links++;

			if (g_number_of_links % 1000 == 0)
				cout << "reading " << g_number_of_links << " links.. " << endl;
		}
	}
	else
	{
		cout << "input_link.csv is not opened." << endl;
		g_ProgramStop();
	}


	cout << "number of links = " << g_number_of_links << endl;

	fprintf(g_pFileOutputLog, "number of links =,%d\n", g_number_of_links);

	parser_link.CloseCSVFile();

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
			int origin_node_id = 0;
			int destination_node_id = 0;
			int departure_time_in_min = 0;
			int arrival_time_in_min = 0;

			parser_agent.GetValueByFieldName("from_origin_node_id", origin_node_id, false);
			agent.origin_node_id = origin_node_id;
			parser_agent.GetValueByFieldName("to_destination_node_id", destination_node_id, false);
			agent.destination_node_id = destination_node_id;
			parser_agent.GetValueByFieldName("departure_time_in_min", departure_time_in_min,false);
			agent.departure_time_in_min = departure_time_in_min;
			parser_agent.GetValueByFieldName("arrival_time_in_min", arrival_time_in_min,false);
			agent.arrival_time_in_min = arrival_time_in_min;


			parser_agent.GetValueByFieldName("PCE", agent.PCE_factor);
			parser_agent.GetValueByFieldName("path_node_sequence", path_node_sequence_str);

			agent.path_node_id_vector = ParseLineToIntegers(path_node_sequence_str);

			parser_agent.GetValueByFieldName("path_schdule_time_sequence", path_schedule_time_sequence_str);

			agent.path_schedule_time_vector = ParseLineToIntegers(path_schedule_time_sequence_str);

			for (int l = 0;l < g_link_vector.size();l++)
			{
				CCost cost;
				for (int t = 0;t < agent.arrival_time_in_min-agent.departure_time_in_min;t++)
				{
					float current_cost = g_link_vector[l].time_dependent_link_cost_for_ADMM[t];
					cost.time_dependent_link_cost_for_ADMM.push_back(current_cost);
					cost.time_dependent_link_cost_for_LR.push_back(current_cost);
				}
				agent.time_dependent_link_cost_for_ADMM.push_back(cost);
				agent.time_dependent_link_cost_for_LR.push_back(cost);
			}
			g_agent_vector.push_back(agent);
		}
	}

	cout << "number of agents = " << g_agent_vector.size() << endl;

	cout << " Sort agents... " << endl;
	std::sort(g_agent_vector.begin(), g_agent_vector.end());

	// simulation
	for (int a = 0; a < g_agent_vector.size(); a++)
	{
		g_agent_vector[a].agent_vector_seq_no = a;

		g_map_agent_id_to_agent_vector_seq_no[g_agent_vector[a].agent_id] = a;// based on agent_id to find agent_vector
	}
	cout << " Sorting ends. ..." << endl;

	// use absoluate time scale
	g_start_simu_interval_no = g_Simulation_StartTimeInMin * 60 / g_number_of_seconds_per_interval;
	g_end_simu_interval_no = g_start_simu_interval_no + g_number_of_simulation_intervals;

	parser_agent.CloseCSVFile();

	CCSVParser parser_cumulative;
	if (parser_cumulative.OpenCSVFile("input_file.csv", true))   // read agent as demand input 
	{
		while (parser_cumulative.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			CFlow flow;
			parser_cumulative.GetValueByFieldName("link_no", flow.link_no);
			parser_cumulative.GetValueByFieldName("time_index", flow.time_index);
			parser_cumulative.GetValueByFieldName("inflow_rate", flow.time_dependent_inflow_rate);
			g_flow_vector.push_back(flow);
			time_max++;
			total_arrival = total_arrival + flow.time_dependent_inflow_rate;
		}

	}
	parser_cumulative.CloseCSVFile();

}
std::vector<float> departure_time_vector;
class NetworkForSP  // mainly for shortest path calculation
{
public:
	int m_threadNo;  // internal thread number 

	std::list<int>  m_SENodeList;  //scan eligible list as part of label correcting algorithm 

	float** m_node_label_cost;  // label cost 
	int* m_node_predecessor;  // predecessor for nodes
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
		m_node_predecessor = new int[number_of_nodes];
		m_time_predecessor = new int[number_of_nodes];
		m_node_status_array = new int[number_of_nodes];
		m_node_label_cost = new float*[number_of_nodes];
		m_link_predecessor = new int[number_of_links];   // note that, the size is still the number of nodes, as each node has only one link predecessor
		m_link_cost_array = new float*[number_of_links];
		new_to_node_cost = new float*[number_of_nodes];
		for (int i = 0;i < number_of_links;i++)
		{
			m_link_cost_array[i] = new float[time_index];
		}
		for (int j = 0;j < number_of_nodes;j++)
		{
			m_node_label_cost[j] = new float[time_index];
			new_to_node_cost[j] = new float[time_index];
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

	// SEList: scan eligible List implementation: the reason for not using STL-like template is to avoid overhead associated pointer allocation/deallocation
	void SEList_clear()
	{
		m_SENodeList.clear();
	}

	void SEList_push_front(int node)
	{
		m_SENodeList.push_front(node);
	}

	void SEList_push_back(int node)
	{
		m_SENodeList.push_back(node);
	}

	bool SEList_empty()
	{
		return m_SENodeList.empty();
	}

	int SEList_front()
	{
		return m_SENodeList.front();
	}

	void SEList_pop_front()
	{
		m_SENodeList.pop_front();
	}

	int optimal_label_correcting_for_LR(int origin_node, int destination_node, int departure_time, int arrival_time, int agent_id)
		// time-dependent label correcting algorithm with double queue implementation
	{
		int internal_debug_flag = 0;

		if (g_node_vector[origin_node].m_outgoing_node_vector.size() == 0)
		{
			return 0;
		}

		for (int i = 0; i < g_number_of_nodes; i++) //Initialization for all nodes
		{
			m_node_status_array[i] = 0;  // not scanned
			m_node_predecessor[i] = -1;  // pointer to previous NODE INDEX from the current label at current node and time
			m_link_predecessor[i] = -1;  // pointer to previous NODE INDEX from the current label at current node and time
			//by slt
			for (int t = 0;t < time_index;t++)
			{
				m_node_label_cost[i][t] = _MAX_LABEL_COST;
			}
		}
		for (int t = 0;t < time_index;t++)
		{
			m_time_predecessor[t] = -1;
		}

		//Initialization for origin node at the preferred departure time, at departure time, cost = 0, otherwise, the g_A2R_simu_interval at origin node
		
		m_node_label_cost[origin_node][departure_time] = departure_time;
		// initialization
		for (int t = departure_time ;t < arrival_time;t++)
		{
			if (t >= 1)
			{
				for (int i = 0; i < g_node_vector.size(); i++)// for each node
				{
					for (int j = 0;j < g_node_vector[i].m_outgoing_node_vector.size();j++)// for each node's to_node 
					{
						int to_node = g_node_vector[i].m_outgoing_node_vector[j].to_node_seq_no;
						int link_seq_no = g_node_vector[i].m_outgoing_node_vector[j].link_seq_no;
						for (int current_departure_time = 0;current_departure_time < t - departure_time;current_departure_time++)// for each travel time
						{
							int current_travel_time = t - current_departure_time;
							if (current_travel_time == g_link_vector[link_seq_no].time_dependent_travel_time_vector[current_departure_time])
							{
								float current_travel_cost = g_agent_vector[agent_id].time_dependent_link_cost_for_LR[link_seq_no].time_dependent_link_cost_for_LR[current_departure_time];
								new_to_node_cost[to_node][t] = m_node_label_cost[i][current_departure_time] + current_travel_cost;
								if (m_node_label_cost[to_node][t] > new_to_node_cost[to_node][t])// update condition
								{
									m_node_label_cost[to_node][t] = new_to_node_cost[to_node][t];
									m_node_predecessor[to_node] = i;  
									m_link_predecessor[to_node] = g_node_vector[i].m_outgoing_node_vector[j].link_seq_no; 
									m_time_predecessor[t] = current_departure_time;
									cout << m_node_label_cost[to_node][t] << endl;

								}
							}
						}
					}
				}
			}
		}
		

		CAgent* p_agent = &(g_agent_vector[agent_id]);
		int current_node_seq_no;
		int current_link_seq_no=-1;
		int current_time_seq_no;
		p_agent->path_link_seq_no_vector_for_LR.clear();  // reset;
		p_agent->path_node_id_vector_for_LR.clear();  // reset;
		p_agent->time_seq_no_vector_for_LR.clear();
		current_node_seq_no = g_internal_node_seq_no_map[p_agent->destination_node_id];
		int current_travel_time = 0;
		int t = arrival_time-1;
		p_agent->time_seq_no_vector_for_LR.push_back(t);
		while (m_time_predecessor[t] != departure_time)
		{
			//cout << m_time_predecessor[t] << endl;
			t = m_time_predecessor[t];
			p_agent->time_seq_no_vector_for_LR.push_back(t);
		}
		p_agent->time_seq_no_vector_for_LR.push_back(departure_time);
				while (current_node_seq_no >= 0 && p_agent->path_link_seq_no_vector_for_LR.size() != p_agent->time_seq_no_vector_for_LR.size())  // this is valid node 
				{
					if (m_link_predecessor[current_node_seq_no] >= 0)
					{
						current_link_seq_no = m_link_predecessor[current_node_seq_no];
					}
					if (current_link_seq_no >= 0)
					{
						p_agent->path_link_seq_no_vector_for_LR.push_back(current_link_seq_no);
					}

					p_agent->path_node_id_vector_for_LR.push_back(current_node_seq_no);
					// by slt
				
					current_node_seq_no = m_node_predecessor[current_node_seq_no];
				}
				if (p_agent->fixed_path_flag != 1)
				{
					std::reverse(std::begin(p_agent->path_node_id_vector_for_LR),
						std::end(p_agent->path_node_id_vector_for_LR));
					std::reverse(std::begin(p_agent->path_link_seq_no_vector_for_LR),
						std::end(p_agent->path_link_seq_no_vector_for_LR));
					std::reverse(std::begin(p_agent->time_seq_no_vector_for_LR),
						std::end(p_agent->time_seq_no_vector_for_LR));
				}
				//if (destination_node >= 0 && m_node_label_cost[destination_node] < _MAX_LABEL_COST)
					//return 1;
				 if (destination_node == -1)
					return 1;  // one to all shortest path
				else
					return -1;


	}

	int optimal_label_correcting_for_ADMM(int origin_node, int destination_node, int departure_time, int arrival_time, int agent_id)
		// time-dependent label correcting algorithm with double queue implementation
	{
		int internal_debug_flag = 0;
		for (int i = 0;i < g_node_vector.size();i++)
		{
			if (g_node_vector[i].external_node_id == origin_node+1)
			{
				origin_node = g_node_vector[i].node_seq_no;
			}
			if (g_node_vector[i].external_node_id == destination_node + 1)
			{
				destination_node = g_node_vector[i].node_seq_no;
			}
		}
		if (g_node_vector[origin_node].m_outgoing_node_vector.size() == 0)
		{
			return 0;
		}

		for (int i = 0; i < g_number_of_nodes; i++) //Initialization for all nodes
		{
			m_node_predecessor[i] = -1;  // pointer to previous NODE INDEX from the current label at current node and time
			m_link_predecessor[i] = -1;  // pointer to previous NODE INDEX from the current label at current node and time
										 //by slt
			for (int t = 0;t < time_index;t++)
			{
				m_node_label_cost[i][t] = _MAX_LABEL_COST;
			}
		}
		for (int t = 0;t < time_index;t++)
		{
			m_time_predecessor[t] = -1;
		}

		//Initialization for origin node at the preferred departure time, at departure time, cost = 0, otherwise, the g_A2R_simu_interval at origin node

		m_node_label_cost[origin_node][departure_time] = departure_time;
		// initialization
		for (int t = departure_time;t < arrival_time+1;t++)
		{
			if (t >= 1)
			{
				for (int i = 0; i < g_node_vector.size(); i++)// for each node
				{
					for (int j = 0;j < g_node_vector[i].m_outgoing_node_vector.size();j++)// for each node's to_node 
					{
						int to_node = g_node_vector[i].m_outgoing_node_vector[j].to_node_seq_no;
						int link_seq_no = g_node_vector[i].m_outgoing_node_vector[j].link_seq_no;

						int current_t=0;
						int test_value=100;
						for (int t_test = t;t_test > 0;t_test--)
						{
							if (m_node_label_cost[i][t_test] <= test_value)
							{
								current_t = t_test;
								test_value = m_node_label_cost[i][t_test];
							}
						}
						for (int current_departure_time = current_t;current_departure_time < t - departure_time;current_departure_time++)// for each travel time
						{
							int current_travel_time = t - current_departure_time;
							if (current_travel_time == g_link_vector[link_seq_no].time_dependent_travel_time_vector[current_departure_time])
							{
								float current_travel_cost = g_agent_vector[agent_id].time_dependent_link_cost_for_ADMM[link_seq_no].time_dependent_link_cost_for_ADMM[current_departure_time];
								new_to_node_cost[to_node][t] = m_node_label_cost[i][current_departure_time] + current_travel_cost;
								if (m_node_label_cost[to_node][t] > new_to_node_cost[to_node][t])// update condition
								{
									m_node_label_cost[to_node][t] = new_to_node_cost[to_node][t];
									m_node_predecessor[to_node] = i;
									m_link_predecessor[to_node] = g_node_vector[i].m_outgoing_node_vector[j].link_seq_no;
									m_time_predecessor[t] = current_departure_time;
									//cout << m_node_label_cost[to_node][t] << endl;

									int j1 = m_node_label_cost[to_node][t];
									int j2 = m_node_label_cost[to_node][t];

								}
							}
						}
					}
				}
			}
		}

		
		CAgent* p_agent = &(g_agent_vector[agent_id]);
		int current_node_seq_no;
		int current_link_seq_no = -1;
		int current_time_seq_no;
		p_agent->path_link_seq_no_vector_for_ADMM.clear();  // reset;
		p_agent->path_node_id_vector_for_ADMM.clear();  // reset;
		p_agent->time_seq_no_vector_for_ADMM.clear();
		current_node_seq_no = g_internal_node_seq_no_map[p_agent->destination_node_id];
		int current_travel_time = 0;
		int t = arrival_time;
		p_agent->time_seq_no_vector_for_ADMM.push_back(t);
		m_node_predecessor[0] = 0;

		while (m_time_predecessor[t] != departure_time)
		{
			//cout << m_time_predecessor[t] << endl;
			t = m_time_predecessor[t];
			p_agent->time_seq_no_vector_for_ADMM.push_back(t);
		}
		p_agent->time_seq_no_vector_for_ADMM.push_back(departure_time);
		while (current_node_seq_no >= 0 && p_agent->path_node_id_vector_for_ADMM.size() != p_agent->time_seq_no_vector_for_ADMM.size())  // this is valid node 
		{
			p_agent->path_node_id_vector_for_ADMM.push_back(current_node_seq_no);
			// by slt
			current_node_seq_no = m_node_predecessor[current_node_seq_no];
		}
		current_node_seq_no = g_internal_node_seq_no_map[p_agent->destination_node_id];
		while (current_node_seq_no >= 0 && p_agent->path_link_seq_no_vector_for_ADMM.size() != p_agent->time_seq_no_vector_for_ADMM.size() - 1) 
		{
			if (m_link_predecessor[current_node_seq_no] >= 0)
			{
				current_link_seq_no = m_link_predecessor[current_node_seq_no];
			}
			if (current_link_seq_no >= 0)
			{
				p_agent->path_link_seq_no_vector_for_ADMM.push_back(current_link_seq_no);
			}
			// by slt
			current_node_seq_no = m_node_predecessor[current_node_seq_no];
		}

		if (p_agent->fixed_path_flag != 1)
		{
			std::reverse(std::begin(p_agent->path_node_id_vector_for_ADMM),
				std::end(p_agent->path_node_id_vector_for_ADMM));
			std::reverse(std::begin(p_agent->path_link_seq_no_vector_for_ADMM),
				std::end(p_agent->path_link_seq_no_vector_for_ADMM));
			std::reverse(std::begin(p_agent->time_seq_no_vector_for_ADMM),
				std::end(p_agent->time_seq_no_vector_for_ADMM));
		}
		//if (destination_node >= 0 && m_node_label_cost[destination_node] < _MAX_LABEL_COST)
		//return 1;
		if (destination_node == -1)
			return 1;  // one to all shortest path
		else
			return -1;


	}

};

class FLow_based_NetworkForSP  // mainly for shortest path calculation
{
public:
	int m_max_cumulative_arrival, m_max_cumulative_departure, m_time_index;
	int*** m_waiting_label_cost;  // label cost 
	int*** new_waiting_cost;
	int*** pred_A;
	int*** pred_D;
	int*** outflow_rate;
	int*** queue_length;
	int** cumulative_departure;

	void AllocateMemory(int max_cumulative_arrival, int max_cumulative_departure, int time_index)
	{
		m_max_cumulative_arrival = max_cumulative_arrival;
		m_max_cumulative_departure = max_cumulative_departure;
		m_time_index = time_index;

		m_waiting_label_cost = Allocate3DDynamicArray<int>(m_max_cumulative_arrival, m_max_cumulative_departure, m_time_index);
		new_waiting_cost = Allocate3DDynamicArray<int>(m_max_cumulative_arrival, m_max_cumulative_departure, m_time_index);
		pred_A = Allocate3DDynamicArray<int>(m_max_cumulative_arrival, m_max_cumulative_departure, m_time_index);
		pred_D = Allocate3DDynamicArray<int>(m_max_cumulative_arrival, m_max_cumulative_departure, m_time_index);
		outflow_rate = Allocate3DDynamicArray<int>(m_max_cumulative_arrival, m_max_cumulative_departure, m_time_index);
		queue_length = Allocate3DDynamicArray<int>(m_max_cumulative_arrival, m_max_cumulative_departure, m_time_index);
		cumulative_departure = new int*[max_cumulative_arrival];
		for (int t = 0;t < time_index;t++)
		{
			cumulative_departure[t] = new int[time_index];
		}


		for (int ca = 0; ca < total_arrival;ca++)
		{
			for (int cd = 0;cd < total_arrival;cd++)
			{
				for (int t = 0;t < time_max;t++)
				{
					if (ca == 0 && cd == 0 && t == 0)
					{
						m_waiting_label_cost[ca][cd][t] = 0;
						queue_length[ca][cd][t] = 0;
						cumulative_departure[ca][t] = 0;
					}
					m_waiting_label_cost[ca][cd][t] = _MAX_LABEL_COST;
					pred_A[ca][cd][t] = -1;
					pred_D[ca][cd][t] = -1;
					outflow_rate[ca][cd][t] = 0;
				}
			}
		}
	}
	~FLow_based_NetworkForSP()
	{
		Deallocate3DDynamicArray<int>(m_waiting_label_cost, m_max_cumulative_arrival, m_max_cumulative_departure);
		Deallocate3DDynamicArray<int>(pred_A, m_max_cumulative_arrival, m_max_cumulative_departure);
		Deallocate3DDynamicArray<int>(pred_D, m_max_cumulative_arrival, m_max_cumulative_departure);
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
//parameters of LR
int g_number_of_LR_iterations = 20;
int g_number_of_ADMM_iterations = 100;
int g_CurrentLRIterationNumber = 0;
int g_Number_Of_Iterations_With_Memory = 5;

float g_best_upper_bound = 99999;
float g_best_lower_bound = -99999;
float g_stepSize = 0;
float g_penalty_PHO = 1;
float g_minimum_subgradient_step_size = 0.1;

NetworkForSP* g_pNetworkForSP_LR = NULL;
NetworkForSP* g_pNetworkForSP_ADMM = NULL;

int LR_iteration_no = 50;
void output_ADMM_files()
{
	FILE* g_pFilecumulative_ADMM = NULL;
	g_pFilecumulative_ADMM = fopen("output_ADMM_trajectory.csv", "w");
	if (g_pFilecumulative_ADMM == NULL)
	{
		cout << "File output_cumulative.csv cannot be opened." << endl;
		g_ProgramStop();
	}
	else
	{
		fprintf(g_pFilecumulative_ADMM, "agent,time_seq,node_seq,link_seq\n");
		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			fprintf(g_pFilecumulative_ADMM, "%d,",a);

			for (int i = 0;i < g_agent_vector[a].time_seq_no_vector_for_ADMM.size();i++)
			{
				int external_time = g_agent_vector[a].time_seq_no_vector_for_ADMM[i];
				fprintf(g_pFilecumulative_ADMM, "%d;", external_time);
			}
			fprintf(g_pFilecumulative_ADMM, ",");

			for (int i = 0;i < g_agent_vector[a].path_node_id_vector_for_ADMM.size();i++)
			{
				int external_node = g_agent_vector[a].path_node_id_vector_for_ADMM[i];
				fprintf(g_pFilecumulative_ADMM, "%d;", external_node);
			}
			fprintf(g_pFilecumulative_ADMM, ",");

			for (int i = 0;i < g_agent_vector[a].path_link_seq_no_vector_for_ADMM.size();i++)
			{
				int external_link = g_agent_vector[a].path_link_seq_no_vector_for_ADMM[i];
				fprintf(g_pFilecumulative_ADMM, "%d;", external_link);
			}
			fprintf(g_pFilecumulative_ADMM, "\n");
		}
		fclose(g_pFilecumulative_ADMM);

	}

}

void find_sp_for_each_agent_use_LR()
{
	g_pNetworkForSP_LR = new NetworkForSP; // create n copies of network, each for a subset of agents to use	
	g_pNetworkForSP_LR->AllocateMemory(g_number_of_nodes, g_number_of_links,time_index);
	for (int a = 0; a < g_agent_vector.size(); a++)
	{
		for (int l = 0; l < g_number_of_links; l++)
		{
			for (int t = 0;t < time_index;t++)
			{
				g_agent_vector[a].time_dependent_link_cost_for_LR[l].time_dependent_link_cost_for_ADMM[t] = g_agent_vector[a].time_dependent_link_cost_for_LR[l].time_dependent_link_cost_for_ADMM[t] + g_link_vector[l].time_dependent_LR_multiplier_vector[t];
			}
		}
		g_pNetworkForSP_LR->optimal_label_correcting_for_LR(g_agent_vector[a].origin_node_id - 1, g_agent_vector[a].destination_node_id - 1, g_agent_vector[a].departure_time_in_min,g_agent_vector[a].arrival_time_in_min,a);
	}
}

void LR_process()
{
	float stepsize;
	g_pNetworkForSP_ADMM = new NetworkForSP; // create n copies of network, each for a subset of agents to use	
	g_pNetworkForSP_ADMM->AllocateMemory(g_number_of_nodes, g_number_of_links, time_index);
	for (int k = 1; k < 5; k++)
	{
		//LR process
		cout << "iteration " << k << endl;
		stepsize = 1.0 / k;
		//update LR_multilier
		for (int l = 0;l < g_link_vector.size();l++)
		{
			for (int t = 0;t < time_index-1;t++)
			{
				int current_v = 0;
				for (int a = 0;a < g_agent_vector.size();a++)
				{
					for (int no = 0;no < g_agent_vector[a].path_link_seq_no_vector_for_LR.size();no++)
					{
						if (g_agent_vector[a].path_link_seq_no_vector_for_LR[no] == g_link_vector[l].link_seq_no && g_agent_vector[a].time_seq_no_vector_for_LR[no] == t)
						{
							current_v++;
						}
					}
				}
				g_link_vector[l].time_dependent_LR_multiplier_vector[t] = max(0, g_link_vector[l].time_dependent_LR_multiplier_vector[t] + (current_v - g_link_vector[l].link_capacity)*stepsize);
			}
		}
		//find_sp_for_each_agent_use_LR();
		//ADMM process
		//set mu
		for (int l = 0;l < g_link_vector.size();l++)
		{
			for (int t = 0;t < time_index - 1;t++)
			{
				int current_v = 0;
				for (int a = 0;a < g_agent_vector.size();a++)
				{
					for (int no = 0;no < g_agent_vector[a].path_link_seq_no_vector_for_ADMM.size();no++)
					{
						if (g_agent_vector[a].path_link_seq_no_vector_for_ADMM[no] == g_link_vector[l].link_seq_no && g_agent_vector[a].time_seq_no_vector_for_ADMM[no] == t)
						{
							current_v++;
						}
					}
				}
				g_link_vector[l].time_dependent_LR_multiplier_vector_for_ADMM[t] = max(0, g_link_vector[l].time_dependent_LR_multiplier_vector_for_ADMM[t] + (current_v - g_link_vector[l].link_capacity)*stepsize);
			}
		}

		for (int a = 0; a < g_agent_vector.size(); a++)
		{
			for (int l = 0;l < g_link_vector.size();l++)
			{
				for (int t = g_agent_vector[a].departure_time_in_min;t <g_agent_vector[a].arrival_time_in_min;t++)
				{
					int mu = 0;
					for (int a1 = 0;a1 < g_agent_vector.size();a1++)
					{
						if (a1 != a)
						{
							for (int no = 0;no < g_agent_vector[a1].path_link_seq_no_vector_for_ADMM.size();no++)
							{
								if (g_agent_vector[a1].path_link_seq_no_vector_for_ADMM[no] == g_link_vector[l].link_seq_no&&g_agent_vector[a1].time_seq_no_vector_for_ADMM[no] == t)
								{
									mu++;
								}
							}
						}
					}
					g_agent_vector[a].time_dependent_link_cost_for_ADMM[l].time_dependent_link_cost_for_ADMM[t] =  g_link_vector[l].link_cost+ g_link_vector[l].time_dependent_LR_multiplier_vector_for_ADMM[t] + max(0, g_penalty_RHO / 2 + g_penalty_PHO*(mu - g_link_vector[l].link_capacity));
					//cout << "agent "<<a<<" at link " << l << "at time " << t << "cost is " << g_agent_vector[a].time_dependent_link_cost_for_ADMM[l].time_dependent_link_cost_for_ADMM[t] << endl;
				}
			}
			g_pNetworkForSP_ADMM->optimal_label_correcting_for_ADMM(g_agent_vector[a].origin_node_id - 1, g_agent_vector[a].destination_node_id - 1, g_agent_vector[a].departure_time_in_min, g_agent_vector[a].arrival_time_in_min, a);
		}
	}

	output_ADMM_files();
}


float k_critical = 4;
float k_jam = 8;

class information
{
public:
	int pred_CA_positon;
	int current_no;
	int value;
};
class CA
{
public:
	std::vector<information> m_CA_vector;
	std::vector<information> m_CD_vector;
	std::vector<information> m_outflow_rate_vector;
	std::vector<information> queue_length;
	std::vector<information> current_cost;
};
class CA_new
{
public:
	int cumulative_arrival;
	int cumulative_departure;
	int outflow_rate;
	int queue_length;
	int current_cost;
	int vacant_no;
};
std::vector<int>CA_vector;
std::vector<int>CD_vector;
void output_cumulative_files(int link_no)
{
	FILE* g_pFilecumulative_CA = NULL;
	g_pFilecumulative_CA = fopen("output_CAfile.csv", "w");
	if (g_pFilecumulative_CA == NULL)
	{
		cout << "File output_cumulative.csv cannot be opened." << endl;
		g_ProgramStop();
	}
	else
	{
		fprintf(g_pFilecumulative_CA, "t,CA,link_id\n");
		for (int i = 0; i < CA_vector.size(); i++)
		{
			fprintf(g_pFilecumulative_CA, "%d,%d,%d\n",i,CA_vector[i],link_no);
		}
		fclose(g_pFilecumulative_CA);

	}
	FILE* g_pFilecumulative_CD = NULL;
	g_pFilecumulative_CD = fopen("output_CDfile.csv", "w");
	if (g_pFilecumulative_CD == NULL)
	{
		cout << "File output_cumulative.csv cannot be opened." << endl;
		g_ProgramStop();
	}
	else
	{
		fprintf(g_pFilecumulative_CD, "t,CD,link_id\n");
		for (int i = 0; i < CD_vector.size(); i++)
		{
			fprintf(g_pFilecumulative_CD, "%d,%d,%d\n", i, CD_vector[i],link_no);
		}
		fclose(g_pFilecumulative_CD);
	}

}

FLow_based_NetworkForSP* g_pFLow_based_NetworkForSP = NULL;
void Cumulative_Curve_Based_DP()
{
	//initialization
	std::vector<int>link_no_sequence;
	link_no_sequence.push_back(1);
	int max_inflow = 5;
	int max_cap = 4;
	for (int link_no = 0;link_no < link_no_sequence.size();link_no++)
	{
		std::vector<CA> cumulative_vector;
		std::vector<int> inflow_rate;
		for (int i = 0;i <= max_inflow;i++)
		{
			inflow_rate.push_back(i);
		}

		for (int t = 0;t < g_flow_vector.size() + 1;t++)
		{
			CA CA;
			if (t == 0)
			{
				information information;
				information.value = 0;
				CA.m_CA_vector.push_back(information);
				CA.m_CD_vector.push_back(information);
				CA.m_outflow_rate_vector.push_back(information);
				CA.queue_length.push_back(information);
				CA.current_cost.push_back(information);
			}
			cumulative_vector.push_back(CA);
		}
		for (int i = 0;i <= 5;i++)
		{
			inflow_rate[i] = i;
		}


		for (int t = 0;t < g_flow_vector.size();t++)// for each t
		{

			for (int no_of_CA = 0;no_of_CA < cumulative_vector[t].m_CA_vector.size();no_of_CA++)
			{
				//for each CA/CD
				int current_A = cumulative_vector[t].m_CA_vector[no_of_CA].value;
				int current_D = cumulative_vector[t].m_CD_vector[no_of_CA].value;
				int current_outflow_rate = cumulative_vector[t].m_outflow_rate_vector[no_of_CA].value;

				//for each inflow_rate, from 2 to 5;
				for (int no_of_current_inflow = 2; no_of_current_inflow < 5;no_of_current_inflow++)
				{
					int current_no = no_of_CA * 4 + no_of_current_inflow - 2;
					int new_A = current_A + inflow_rate[no_of_current_inflow];
					int new_D = current_D + current_outflow_rate;

					information information_ca;
					information_ca.value = new_A;
					information_ca.current_no = current_no;
					information_ca.pred_CA_positon = int(current_no / 4);
					if (new_A <= 58)
					{
						cumulative_vector[t + 1].m_CA_vector.push_back(information_ca);
					}
					else
					{
						information_ca.value = 58;
						cumulative_vector[t + 1].m_CA_vector.push_back(information_ca);
					}

					information information_cd;
					information_cd.value = new_D;
					information_cd.current_no = current_no;
					information_cd.pred_CA_positon = int(current_no / 4);
					if (new_D <= 58)
					{
						cumulative_vector[t + 1].m_CD_vector.push_back(information_cd);
					}
					else
					{
						information_cd.value = 58;
						cumulative_vector[t + 1].m_CD_vector.push_back(information_cd);
					}

					int new_queue = new_A - new_D;
					information information_queue;
					information_queue.value = new_queue;
					information_queue.current_no = current_no;
					information_queue.pred_CA_positon = int(current_no / 4);
					cumulative_vector[t + 1].queue_length.push_back(information_queue);

					int new_outflow_rate;
					if (new_queue <= k_critical)
					{
						new_outflow_rate = min(max_cap, new_queue);
					}
					if (new_queue > k_critical)
					{
						new_outflow_rate = max_cap*(k_jam - new_queue) / (k_jam - k_critical);
					}
					information information_outflowrate;
					information_outflowrate.value = new_outflow_rate;
					information_outflowrate.current_no = current_no;
					information_outflowrate.pred_CA_positon = int(current_no / 4);
					cumulative_vector[t + 1].m_outflow_rate_vector.push_back(information_outflowrate);

					information information_cost;
					int current_cost = new_queue * 1;
					information_cost.value = current_cost + cumulative_vector[t].current_cost[int(current_no / 4)].value;
					information_cost.current_no = current_no;
					information_cost.pred_CA_positon = int(current_no / 4);
					cumulative_vector[t + 1].current_cost.push_back(information_cost);
				}
			}
		}
		//generate a CA_vector

		int minimum_cost = 1000;
		int minimun_no;
		int pred_ca;
		for (int t = g_flow_vector.size();t >= 0;t--)
		{
			if (t == g_flow_vector.size())
			{
				for (int no = 0; no < cumulative_vector[t].current_cost.size();no++)
				{
					if (cumulative_vector[t].current_cost[no].value < minimum_cost)
					{
						minimum_cost = cumulative_vector[t].current_cost[no].value;
						minimun_no = cumulative_vector[t].current_cost[no].current_no;
						pred_ca = cumulative_vector[t].current_cost[no].pred_CA_positon;
					}
				}
				CA_vector.push_back(cumulative_vector[t].m_CA_vector[minimun_no].value);
				CD_vector.push_back(cumulative_vector[t].m_CD_vector[minimun_no].value);
			}
			else
			{
				CA_vector.push_back(cumulative_vector[t].m_CA_vector[pred_ca].value);
				CD_vector.push_back(cumulative_vector[t].m_CD_vector[pred_ca].value);

				pred_ca = cumulative_vector[t].m_CA_vector[pred_ca].pred_CA_positon;
			}
		}
		std::reverse(std::begin(CA_vector),
			std::end(CA_vector));
		std::reverse(std::begin(CD_vector),
			std::end(CD_vector));	
		output_cumulative_files(link_no_sequence[link_no]);
		CA_vector.clear();
		CD_vector.clear();
	}
}


//simulation and adjustment
void Cumulative_Curve_Based_DP_1()
{
	//initialization
	int max_inflow = 4;
	int max_cap=2;
	std::vector<CA_new> cumulative_vector;
	for (int t = 0;t < g_flow_vector.size() + 1;t++)
	{
		CA_new CA;
			CA.cumulative_arrival = 0;
			CA.cumulative_departure = 0;
			CA.current_cost = 0;
			CA.outflow_rate = 0;
			CA.queue_length = 0;
		cumulative_vector.push_back(CA);
	}
	
	//simulation process
	for (int t = 0;t < g_flow_vector.size();t++)// for each t
	{
		int new_A = cumulative_vector[t].cumulative_arrival + g_flow_vector[t].time_dependent_inflow_rate;
		int new_D = cumulative_vector[t].cumulative_departure + cumulative_vector[t].outflow_rate;
		int new_queue = new_A - new_D;
		cumulative_vector[t + 1].cumulative_arrival = new_A;
		cumulative_vector[t + 1].cumulative_departure = new_D;
		cumulative_vector[t + 1].queue_length = new_queue;
		if (new_queue <= k_critical)
		{
			cumulative_vector[t+1].outflow_rate = min(max_cap, new_queue);
		}
		if (new_queue > k_critical)
		{
			cumulative_vector[t + 1].outflow_rate = max_cap*(k_jam - new_queue) / (k_jam - k_critical);
		}
		cumulative_vector[t + 1].vacant_no = max_cap - g_flow_vector[t].time_dependent_inflow_rate;
		cumulative_vector[t + 1].current_cost = cumulative_vector[t].current_cost + cumulative_vector[t + 1].queue_length * 1;

	}
	
	//adjust process
	for (int t = 1;t < g_flow_vector.size();t++)
	{
		if (cumulative_vector[t].outflow_rate < max_cap) // look for the queue
		{

			int vacant_in_need = g_flow_vector[t-1].time_dependent_inflow_rate - max_cap;

			for (int t1 = t+1;t1 < g_flow_vector.size();t1++)
			{
				if (cumulative_vector[t1].vacant_no > 0)
				{
					if (vacant_in_need <= cumulative_vector[t1].vacant_no&&vacant_in_need > 0)
					{
						g_flow_vector[t1].time_dependent_inflow_rate = g_flow_vector[t1].time_dependent_inflow_rate + vacant_in_need;
						g_flow_vector[t].time_dependent_inflow_rate = g_flow_vector[t].time_dependent_inflow_rate - vacant_in_need;
						vacant_in_need = 0;
					}
					else
					{
						g_flow_vector[t1].time_dependent_inflow_rate = g_flow_vector[t1].time_dependent_inflow_rate + cumulative_vector[t1].vacant_no;
						vacant_in_need = vacant_in_need - cumulative_vector[t1].vacant_no;
					}
				}
			}
		}
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

	g_ReadInputData();  // step 1: read input data of network and demand agent 
	//g_TrafficAssignment();
	//g_TrafficSimulation();
	////LR
	////// step 2: initialize the resource matrix in upper bound
	////// step 3: Lagrangian optimization for Lower and Upper bound
	//g_LR_Optimization(0);
	//g_LR_Optimization_based_on_dynamic_state_array();
	//g_OutputFiles();
	LR_process();
	//Cumulative_Curve_Based_DP();
	//cout << "End of Optimization " << endl;
	//cout << "free memory.." << endl;
	//cout << "done." << endl;
	//delete[] pNetworkForSP;

	//g_node_vector.clear();
	//g_link_vector.clear();
	//g_agent_vector.clear();

	cout << "it is done!" << endl;

	return 1;
	
}

