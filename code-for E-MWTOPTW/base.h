#include<bits/stdc++.h>
#include<ilcplex/ilocplex.h>
using namespace std;
#define eps 1e-1
#define sqr(x) ((x)*(x))
typedef IloArray <IloNumVarArray> IloNumVarArray2;
typedef IloArray <IloNumVarArray2> IloNumVarArray3;
int **A;
int *Count;
int bad_times;
int MAXBAD;
const int INF = 0x3f3f3f3f;
int **Distances;// Distance_from_to
int **TimeCost;
int *CostOfWorkers;
int *TabuTable;
int N;// Number of customer
int M; // Maximum NumberOfWorks
double UnitChargeCost; // 单位时间充电成本
int MAXIteration;
int Capicity;
int K;
int speed;
int H;
double ChargeSpeed;
const int TabuLength = 10;
struct Vertex{
	int id;
	int Coord_x, Coord_y;
	int TimeWindowL, TimeWindowR;
	int Demand, Revenue;
	int serviceTime; //这个与需求和雇佣人数相关

	void ReadInfo(FILE *fp){
		fscanf(fp, "%d %d %d %d %d %d %d",&id, &Coord_x, &Coord_y, &TimeWindowL, &TimeWindowR, &Demand, &Revenue);
	}

	int CalcServiceTime(int NumberOfWorks){
		double decent = pow((double)NumberOfWorks, 0.9);
		return serviceTime = (int)((double)Demand / decent + 0.5);
	}// serviceTime = Demand / (N ^ 0.9)

	//friend int CalcDistance(const Vertex &a, const Vertex &b);

};

Vertex *Customers;

struct ServeNode{
	Vertex Customer;
	int EF, LS; // 最早结束时间和最晚开始时间
	int NumberOfWorks; // 确定雇佣的工人数量
	double TimeOfCharge; // 充电时间
};
typedef vector<ServeNode> ServeList;
ServeList AllCustomers;

inline int CalcDistance(const Vertex &a, const Vertex &b){
	return (int)(sqrt(sqr(a.Coord_x - b.Coord_x) + sqr(a.Coord_y - b.Coord_y)));
}

struct Vehicle{
    ServeList path;
	int Revenue; // 总收入
	double Porfit; // 总利润
	int Distance; // 总行驶距离
	double ChargeCost; // 总充电成本
	int EmpolyCost; // 总雇佣成本
	int CurrentCapicity; //　当前容量

    void GenerateNewVeh(){
        path.push_back(AllCustomers[0]);
        path.push_back(AllCustomers[0]);
        Revenue = 0;
        Porfit = 0;
        Distance = 0;
        ChargeCost = 0;
        EmpolyCost = 0;
		CurrentCapicity = 0;
        return;
    }

	int CalcRevenue(){
		Revenue = 0;
		for (size_t i = 0; i < path.size(); i ++){
			Revenue += path[i].Customer.Revenue;
		}
		return Revenue;
	}

	int CalcDistance(){
		Distance = 0;
		for (size_t i = 1; i < path.size(); i ++){
			Distance += Distances[path[i-1].Customer.id][path[i].Customer.id];
		}
		return Distance;
	}

	double CalcChargeCost(){
		ChargeCost = 0;
		for (size_t i = 0; i < path.size(); i ++){
			ChargeCost += path[i].TimeOfCharge * UnitChargeCost;
		}
		return ChargeCost;
	}

	int CalcEmpolyCost(){
		EmpolyCost = 0;
		CurrentCapicity = 0;
		for (size_t i = 0; i < path.size(); i ++){
			EmpolyCost += CostOfWorkers[path[i].NumberOfWorks];
			CurrentCapicity += path.at(i).Customer.Demand;
		}

		return EmpolyCost;
	}

	double CalcPorfit(){
		CalcEmpolyCost();
		CalcRevenue();
		CalcChargeCost();
		CalcDistance();
		Porfit = Revenue - ChargeCost - EmpolyCost;
		return Porfit;
	}

    void CalcEFandLS(){
		//最早结束时间
		path[0].EF = path[0].Customer.TimeWindowL;
        for (size_t i = 1; i < path.size(); i ++){
            path[i].EF = max(path[i].Customer.TimeWindowL
							+ path[i].Customer.serviceTime,
							path[i-1].EF + path[i].Customer.serviceTime
							+ TimeCost[path[i-1].Customer.id][path[i].Customer.id]);
        }

		//最晚开始时间
		int n = path.size() - 1;
		path[n].LS = path[n].Customer.TimeWindowR;
		for (int i = n - 1; i >= 0; i --){
			path[i].LS = min(path[i].Customer.TimeWindowR - path[i].Customer.serviceTime,
							 path[i + 1].LS - TimeCost[path[i].Customer.id][path[i+1].Customer.id]
						 	 - path[i].Customer.serviceTime);
		}
    }


	bool InsertInto(ServeNode &L){
		bool SuccessfulInsert = 0;
		CalcEFandLS();
		if (L.Customer.Demand + CurrentCapicity > Capicity) return SuccessfulInsert;
		for (int l = 1; l <= M; l ++){
            if (CostOfWorkers[l] > L.Customer.Revenue) break;
			if (SuccessfulInsert) break;
			int s = L.Customer.CalcServiceTime(l);
			for (size_t i = 0; i < path.size() - 1; i ++){
				L.EF = max(L.Customer.TimeWindowL + s, path[i].EF + TimeCost[path[i].Customer.id][L.Customer.id] + s);
				if (L.EF > L.Customer.TimeWindowR) continue;
				int TD = path[i+1].LS - path[i].EF;
				int TC = TimeCost[path[i].Customer.id][L.Customer.id] + s + TimeCost[L.Customer.id][path[i + 1].Customer.id];
				int TW = max(0, L.Customer.TimeWindowL - TimeCost[path[i].Customer.id][L.Customer.id] - path[i].EF);
				if (TD < TC + TW) continue;
				L.NumberOfWorks = l;
				L.Customer.CalcServiceTime(l);
				CurrentCapicity += L.Customer.Demand;
				SuccessfulInsert = 1;
				path.insert(path.begin() + i + 1, L);
				break;
			}
		}
		return SuccessfulInsert;
	}

    int GreedyInsert(ServeNode &L, int &LocationID, double &Gprofit, bool RealInsert = 0){
        bool SuccessfulInsert = 0;
		CalcEFandLS();
		if (L.Customer.Demand + CurrentCapicity > Capicity) return INF;
        int final_L, final_i;
        int MinDelta = INF;
        Gprofit = 0;
        double tmpprofit = 0;
        for (int l = 1; l <= M; l ++){
            if (CostOfWorkers[l] > L.Customer.Revenue)  return INF;
			int s = L.Customer.CalcServiceTime(l);
			for (size_t i = 0; i < path.size() - 1; i ++){
				L.EF = max(L.Customer.TimeWindowL + s, path[i].EF + TimeCost[path[i].Customer.id][L.Customer.id] + s);
				if (L.EF > L.Customer.TimeWindowR) continue;
				int TD = path[i+1].LS - path[i].EF;
				int TC = TimeCost[path[i].Customer.id][L.Customer.id] + s + TimeCost[L.Customer.id][path[i + 1].Customer.id];
				int TW = max(0, L.Customer.TimeWindowL - TimeCost[path[i].Customer.id][L.Customer.id] - path[i].EF);
				if (TD < TC + TW) continue;
                int DeltaDistance = - Distances[path.at(i).Customer.id][path.at(i+1).Customer.id]
                                    + Distances[path.at(i).Customer.id][L.Customer.id]
                                    + Distances[L.Customer.id][path.at(i+1).Customer.id];
                tmpprofit = L.Customer.Revenue - CostOfWorkers[l];
                if (DeltaDistance > MinDelta && Gprofit >= tmpprofit) continue;
                else{
                    MinDelta = DeltaDistance;
                    Gprofit = tmpprofit;
                    final_i = i;
                    final_L = l;
                    SuccessfulInsert = 1;
                }
			}
			if (SuccessfulInsert) break;
		}

        if (SuccessfulInsert && RealInsert){
            L.Customer.CalcServiceTime(final_L);
            L.NumberOfWorks = final_L;
            CurrentCapicity += L.Customer.Demand;
            path.insert(path.begin() + final_i + 1, L);
            LocationID = final_i;
        }

		return MinDelta;
    }

    void UpdateInsert(ServeNode &L, int final_i){
        Revenue += L.Customer.Revenue;
        Distance += Distances[path.at(final_i).Customer.id][L.Customer.id]
                    + Distances[L.Customer.id][path.at(final_i + 2).Customer.id]
                    - Distances[path.at(final_i).Customer.id][path.at(final_i + 2).Customer.id];
        EmpolyCost += CostOfWorkers[L.NumberOfWorks];
        Porfit = Revenue - EmpolyCost - ChargeCost;
        return;
    }

    void print(FILE *fp){
        fprintf(fp, "[");
        for (int i = 0; i < path.size()-1; i ++){
            fprintf(fp, "%d,", path[i].Customer.id);
        }
        fprintf(fp, "%d]\n", path.back().Customer.id);
    }

    bool FeasibleCheck(bool debug = 0, bool charge = 0){
        int tmpCapicity = 0;
        int tmpTime = 0;
        int tmpCharge = H;
        for (int i = 1; i < path.size(); i ++){
            tmpTime += TimeCost[path.at(i-1).Customer.id][path.at(i).Customer.id];
            tmpCharge -= Distances[path.at(i-1).Customer.id][path.at(i).Customer.id];
            tmpCharge += path.at(i).TimeOfCharge * ChargeSpeed;
            tmpCapicity += path.at(i).Customer.Demand;
            if (tmpCharge < 0 && charge == 1) {if (debug) printf("power\n"); return 0;}
            if (tmpCapicity > Capicity) {if (debug) {printf("Capicity\n"); printf("%d %d\n", tmpCapicity, CurrentCapicity);} return 0;}
            if (tmpTime < path.at(i).Customer.TimeWindowL) tmpTime = path.at(i).Customer.TimeWindowL;
            tmpTime += path.at(i).Customer.serviceTime;
            if (tmpTime > path.at(i).Customer.TimeWindowR) {if (debug) printf("Timewindow\n"); return 0;}
        }
        return 1;
    }
};
vector <Vehicle> SolutionsPool;

struct Solution{
	vector <Vehicle> Vehicles;
	int TotalRevenue; // 总收入
	double TotalPorfit; // 总利润
	int TotalDistance; // 总行驶距离
	double TotalChargeCost; // 总充电成本
    deque <ServeNode> NotVisited;
	int TotalEmpolyCost; // 总雇佣成本y
	double CalcTotalPorfit(){
        TotalPorfit = 0;
        TotalRevenue = 0;
        TotalEmpolyCost = 0;
        TotalDistance = 0;
        TotalChargeCost = 0;
		for (size_t i = 0; i < Vehicles.size(); i ++){
			Vehicles[i].CalcPorfit();
			TotalPorfit += Vehicles[i].Porfit;
			TotalRevenue += Vehicles[i].Revenue;
			TotalDistance += Vehicles[i].Distance;
			TotalChargeCost += Vehicles[i].ChargeCost;
			TotalEmpolyCost += Vehicles[i].EmpolyCost;
		}
		return TotalPorfit;
	}

    void buildUnvisit(){
        NotVisited.clear();
        while (NotVisited.size()) NotVisited.pop_back();
        vector<bool> vis;
        for (int i = 0; i <= N; i ++) vis.push_back(0);
        for (int i = 0; i < Vehicles.size(); i ++){
            Vehicle V = Vehicles[i];
            for (int j = 0; j < V.path.size(); j ++){
                vis[V.path.at(j).Customer.id] = 1;
            }
        }
        for (int i = 1; i <= N; i ++){
            if (!vis[i]){
               NotVisited.push_back(AllCustomers[i]);
            }
        }
    }

	void GenerateInitSolution(){
        Vehicles.clear();
        Vehicles.shrink_to_fit();
        while (NotVisited.size()) NotVisited.pop_back();
		Vehicle emptyV;
		emptyV.GenerateNewVeh();
		for (int k = 0; k < K; k ++)
			Vehicles.push_back(emptyV);
		bool InsertSucess = 0;
		bool RandInsert = rand() % 2;
		for (int i = 1; i <= N; i ++){
            InsertSucess = 0;
			for (size_t t = 0; t < Vehicles.size(); t ++){
                RandInsert = rand() % 2;
                if (!RandInsert){
                    continue;
                }
				InsertSucess = Vehicles[t].InsertInto(AllCustomers[i]);
				if (InsertSucess) break;
			}
			if (!InsertSucess) {
                NotVisited.push_back(AllCustomers[i]);
            }
		}
	}

	void UpdateDec(int VehicleID){
        TotalChargeCost -= Vehicles[VehicleID].ChargeCost;
        TotalDistance -= Vehicles[VehicleID].Distance;
        TotalEmpolyCost -= Vehicles[VehicleID].EmpolyCost;
        TotalPorfit -= Vehicles[VehicleID].Porfit;
        TotalRevenue -= Vehicles[VehicleID].Revenue;
        return;
	}

	void UpdateAdd(int VehicleID){
        TotalChargeCost += Vehicles[VehicleID].ChargeCost;
        TotalDistance += Vehicles[VehicleID].Distance;
        TotalEmpolyCost += Vehicles[VehicleID].EmpolyCost;
        TotalPorfit += Vehicles[VehicleID].Porfit;
        TotalRevenue += Vehicles[VehicleID].Revenue;
        return;
	}

    bool GreedyInsert(ServeNode &L, int &VehicleID, int &LocationID){
        int MinDelta = INF;
        int FinalV = -1;
        double Gporfit;
        double MaxPorfit = 0;
        for (int i = 0; i < Vehicles.size(); i ++){
            int Delta = Vehicles.at(i).GreedyInsert(L, LocationID, Gporfit);
            if (Delta < MinDelta || Gporfit > MaxPorfit){
                FinalV = i;
                MinDelta = Delta;
                MaxPorfit = Gporfit;
            }
        }
        if (MinDelta == INF) return 0;
        Vehicles.at(FinalV).GreedyInsert(L, LocationID, Gporfit, 1);
        VehicleID = FinalV;
        return 1;
    }

    void Eliminator(){
        const double Plow = 0.2;
        const double Phigh = 0.8;
        double MeanProfit = TotalPorfit / (N - NotVisited.size() * 1.0);
        for (int i = 0; i < Vehicles.size(); i ++){
            Vehicle &L = Vehicles[i];
            for (int j = 1; j < L.path.size() - 1; j ++){
                double profit = L.path.at(j).Customer.Revenue;
                profit -= CostOfWorkers[L.path.at(j).NumberOfWorks];
                profit -= L.path.at(j).TimeOfCharge * UnitChargeCost;
                double prob = (rand() % 10000) / (1.0 * 10000);
                if (profit >= MeanProfit){
                    if (prob < Plow){
                        NotVisited.push_back(L.path[j]);
                        L.CurrentCapicity -= L.path.at(j).Customer.Demand;
                        L.path.erase(j + L.path.begin());
                        j --;
                    }
                    continue;
                }

                if (profit < MeanProfit){
                    if (prob < Phigh){
                        NotVisited.push_back(L.path[j]);
                        L.CurrentCapicity -= L.path.at(j).Customer.Demand;
                        L.path.erase(j + L.path.begin());
                        j --;
                    }
                    continue;
                }
            }
        }

        int NowSize = NotVisited.size();
        bool SuccessInsert = 0;
        int VehicleID;
        int LocationID;
        while (NowSize --){
            ServeNode L = NotVisited.front();
            NotVisited.pop_front();
            SuccessInsert = GreedyInsert(L, VehicleID, LocationID);
            if (SuccessInsert){
                continue;
            }
            else {
                NotVisited.push_back(L);
            }
        }
        CalcTotalPorfit();
    }


    bool ZeroRelocate(){
        int NowSize = NotVisited.size();
        bool SuccessInsert = 0;
        int VehicleID;
        int LocationID;
        while (NowSize --){
            ServeNode L = NotVisited.front();
            NotVisited.pop_front();
            SuccessInsert = GreedyInsert(L, VehicleID, LocationID);
            if (SuccessInsert){
                continue;
            }
            else NotVisited.push_back(L);

        }
        CalcTotalPorfit();
    }

    void OneRelocate(bool SA = 0, bool Tabu = 0, int iterT = 0){
        for (int i = 0; i < Vehicles.size(); i ++){
            Vehicle &L = Vehicles[i];
            for (int j = 1; j < L.path.size()-1; j ++){
                ServeNode S = L.path.at(j);
                if (Tabu && TabuTable[L.path.at(j).Customer.id] > iterT && L.path.at(j).Customer.Revenue - CostOfWorkers[L.path.at(j).NumberOfWorks] > TotalPorfit / (N - NotVisited.size())){
                    continue;
                }
                if (Tabu){
                    TabuTable[L.path.at(j).Customer.id] = iterT + TabuLength;
                }
                L.CurrentCapicity -= L.path.at(j).Customer.Demand;
                L.path.erase(j + L.path.begin());
                j --;
                int LocationID;
                double tmp;
                if (!SA && !Tabu)
                    L.GreedyInsert(S, LocationID, tmp, 1);
                else{
                        L.InsertInto(S);
                }
                j++;
            }
        }
        CalcTotalPorfit();
    }

    void TwoRelocate(bool SA = 0, bool Tabu = 0, int iterT = 0){
        for (int i = 0; i < Vehicles.size(); i ++){
            Vehicle &L = Vehicles[i];
            for (int j = 1; j < L.path.size()-1; j ++){
                ServeNode S = L.path.at(j);
                if (Tabu && TabuTable[L.path.at(j).Customer.id] > iterT){
                    continue;
                }
                if (Tabu){
                    TabuTable[L.path.at(j).Customer.id] = iterT + TabuLength;
                }
                L.CurrentCapicity -= L.path.at(j).Customer.Demand;
                L.path.erase(j + L.path.begin());
                j --;
                int LocationID, VehicleID;
                if (!SA && !Tabu)
                    GreedyInsert(S, VehicleID, LocationID);
                else{
                    bool SuccessInsert = 0;
                    for (VehicleID = 0; VehicleID < Vehicles.size(); VehicleID ++){
                        SuccessInsert = Vehicles[VehicleID].InsertInto(S);
                        if (SuccessInsert) break;
                    }
                }
                if (VehicleID == i) j ++;
            }
        }
        CalcTotalPorfit();
    }

    void ZeroExchange(){
        int NowSize = NotVisited.size();
        while (NowSize --){
            ServeNode L = NotVisited.front();
            NotVisited.pop_front();
            bool SuccessSwap = 0;
            for (int i = 0; i < Vehicles.size(); i ++){
                if (SuccessSwap) break;
                Vehicle &V = Vehicles[i];
                V.CalcEFandLS();
                for (int j = 1; j < V.path.size() - 1; j ++){
                    if (SuccessSwap) break;
                    ServeNode &S = V.path[j];

                    for (int l = 1; l <= M; l ++){
                        if (S.Customer.Revenue - CostOfWorkers[S.NumberOfWorks] > L.Customer.Revenue - CostOfWorkers[L.NumberOfWorks]) break;
                        if (V.CurrentCapicity - S.Customer.Demand + L.Customer.Demand > Capicity) break;
                        int ServiceTime = L.Customer.CalcServiceTime(l);
                        if (L.Customer.Revenue - CostOfWorkers[l] < S.Customer.Revenue - CostOfWorkers[S.NumberOfWorks]) break;
                        L.EF = max(V.path.at(j-1).EF + TimeCost[V.path.at(j-1).Customer.id][L.Customer.id] + ServiceTime,
                                    L.Customer.TimeWindowL + ServiceTime);
                        if (L.EF > L.Customer.TimeWindowR) continue;
                        int TD = V.path[j+1].LS - V.path[j-1].EF;
                        int TC = TimeCost[V.path[j-1].Customer.id][L.Customer.id] + ServiceTime + TimeCost[L.Customer.id][V.path[j + 1].Customer.id];
                        int TW = max(0, L.Customer.TimeWindowL - TimeCost[V.path[j-1].Customer.id][L.Customer.id] - V.path[j-1].EF);
                        if (TD >= TC + TW){
                            L.NumberOfWorks = l;
                            V.CurrentCapicity += - S.Customer.Demand + L.Customer.Demand;
                            swap(V.path[j], L);
                            SuccessSwap = 1;
                            break;
                        }
                    }
                }
            }

            NotVisited.push_back(L);

        }
        CalcTotalPorfit();
    }

    void OneExchange(bool SA = 0, bool Tabu = 0, int iterT = 0){
        for (int v = 0; v < Vehicles.size(); v ++){
            Vehicle &V = Vehicles[v];
            int MinDelta = INF;
            int LocationI, LocationJ;
            for (int i = 1; i < V.path.size() - 2; i ++){
                for (int j = i + 1; j < V.path.size()-1; j ++){
                    ServeNode &L1 = V.path[i];
                    ServeNode &L2 = V.path[j];

                    int delta = - Distances[V.path.at(i-1).Customer.id][L1.Customer.id]
                                - Distances[V.path.at(i+1).Customer.id][L1.Customer.id]
                                + Distances[V.path.at(i-1).Customer.id][L2.Customer.id]
                                + Distances[V.path.at(i+1).Customer.id][L2.Customer.id]
                                - Distances[V.path.at(j-1).Customer.id][L2.Customer.id]
                                - Distances[V.path.at(j+1).Customer.id][L2.Customer.id]
                                + Distances[V.path.at(j-1).Customer.id][L1.Customer.id]
                                + Distances[V.path.at(j+1).Customer.id][L1.Customer.id];
                    if (delta >= 0 && !SA && !Tabu) continue;
                    if (delta < MinDelta){
                        MinDelta = delta;
                        LocationI = i;
                        LocationJ = j;
                    }
                }
            }
            if (MinDelta != INF && Tabu && (TabuTable[V.path.at(LocationI).Customer.id] > iterT || TabuTable[V.path.at(LocationJ).Customer.id] > iterT) && MinDelta >= 0) continue;

            if (Tabu && MinDelta != INF){
                TabuTable[V.path.at(LocationI).Customer.id] = iterT + TabuLength;
                TabuTable[V.path.at(LocationJ).Customer.id] = iterT + TabuLength;
            }
            if (MinDelta != INF){
                swap(V.path[LocationI], V.path[LocationJ]);
                if (!V.FeasibleCheck()) swap(V.path[LocationI], V.path[LocationJ]);
            }
        }
        CalcTotalPorfit();
    }


    void TwoExchange(bool SA = 0, bool Tabu = 0, int iterT = 0){
        for (int i = 0; i < Vehicles.size() - 1; i ++){
            for (int j = i + 1; j < Vehicles.size(); j ++){
                Vehicle &V1 = Vehicles[i];
                Vehicle &V2 = Vehicles[j];
                int MinDelta = INF;
                int LocationI, LocationJ;
                for (int u = 1; u < V1.path.size() - 1; u ++){
                    ServeNode &x = V1.path[u];
                    ServeNode &lx = V1.path[u-1], &rx = V1.path[u+1];
                    for (int v = 1; v < V2.path.size() - 1; v ++){
                        ServeNode &y = V2.path[v];
                        ServeNode &ly = V2.path[v-1], &ry = V2.path[v+1];
                        if (V1.CurrentCapicity + y.Customer.Demand - x.Customer.Demand > Capicity
                            || V2.CurrentCapicity - y.Customer.Demand + x.Customer.Demand > Capicity) continue;
                        int delta = -Distances[lx.Customer.id][x.Customer.id]
                                    -Distances[rx.Customer.id][x.Customer.id]
                                    -Distances[ly.Customer.id][y.Customer.id]
                                    -Distances[ry.Customer.id][y.Customer.id]
                                    +Distances[lx.Customer.id][y.Customer.id]
                                    +Distances[rx.Customer.id][y.Customer.id]
                                    +Distances[ly.Customer.id][x.Customer.id]
                                    +Distances[ry.Customer.id][x.Customer.id];
                        if (delta >= 0 && !SA && !Tabu) continue;
                        if (Tabu && TabuTable[x.Customer.id] > iterT || TabuTable[y.Customer.id] > iterT  && MinDelta >= 0) continue;
                        if (MinDelta > delta){
                            MinDelta = delta;
                            LocationI = u;
                            LocationJ = v;
                        }

                    }
                }

                if (MinDelta != INF){
                    V1.CalcEFandLS();
                    V2.CalcEFandLS();
                    ServeNode &x = V1.path[LocationI];
                    ServeNode &lx = V1.path[LocationI-1], &rx = V1.path[LocationI+1];
                    ServeNode &y = V2.path[LocationJ];
                    ServeNode &ly = V2.path[LocationJ-1], &ry = V2.path[LocationJ+1];

                    x.EF = min(ly.EF + TimeCost[x.Customer.id][ly.Customer.id] + x.Customer.serviceTime,
                                x.Customer.TimeWindowL + x.Customer.serviceTime);
                    y.EF = min(lx.EF + TimeCost[y.Customer.id][lx.Customer.id] + y.Customer.serviceTime,
                                y.Customer.TimeWindowL + y.Customer.serviceTime);
                    if (x.EF > x.Customer.TimeWindowR || y.EF > y.Customer.TimeWindowR) continue;
                    int TD1 = rx.LS - lx.EF;
                    int TC1 = TimeCost[lx.Customer.id][y.Customer.id] + y.Customer.serviceTime + TimeCost[rx.Customer.id][y.Customer.id];
                    int TW1 = max(0, y.Customer.TimeWindowL - TimeCost[lx.Customer.id][y.Customer.id] - lx.EF);
                    if (TD1 < TC1 + TW1) continue;
                    int TD2 = ry.LS - ly.EF;
                    int TC2 = TimeCost[ly.Customer.id][x.Customer.id] + x.Customer.serviceTime + TimeCost[ry.Customer.id][x.Customer.id];
                    int TW2 = max(0, x.Customer.TimeWindowL - TimeCost[ly.Customer.id][x.Customer.id] - ly.EF);
                    if (TD2 < TC2 + TW2) continue;
                    if (Tabu){
                        TabuTable[x.Customer.id] = iterT + TabuLength;
                        TabuTable[y.Customer.id] = iterT + TabuLength;
                    }
                    V1.CurrentCapicity -= x.Customer.Demand;
                    V1.CurrentCapicity += y.Customer.Demand;
                    V2.CurrentCapicity -= y.Customer.Demand;
                    V2.CurrentCapicity += x.Customer.Demand;
                    swap(x, y);

                }
            }
        }
        CalcTotalPorfit();
    }

    void TwoOpt(bool SA = 0){
        for (int i = 0; i < Vehicles.size() - 1; i ++){
            Vehicle &V1 = Vehicles[i];
            int Location1I, Location1J;
            int Location2I, Location2J;
            int VehicleID;
            int MinDelta = INF;
            for (int u1 = 1; u1 < V1.path.size() - 2; u1 ++){
                for (int v1 = u1 + 1; v1 < V1.path.size() - 1; v1 ++){
                    ServeNode &x1 = V1.path[u1], &y1 = V1.path[v1];
                    ServeNode &lx1 = V1.path[u1-1], &ry1 = V1.path[v1 + 1];
                    for (int j = i + 1; j < Vehicles.size(); j ++){
                        Vehicle &V2 = Vehicles[j];

                        for (int u2 = 1; u2 < V2.path.size() - 2; u2 ++){
                            for (int v2 = u2 + 1; v2 < V2.path.size() - 1; v2 ++){
                                ServeNode &x2 = V2.path[u2], &y2 = V2.path[v2];
                                ServeNode &lx2 = V2.path[u2-1], &ry2 = V2.path[v2 + 1];

                                int delta = - Distances[x1.Customer.id][lx1.Customer.id]
                                            - Distances[y1.Customer.id][ry1.Customer.id]
                                            + Distances[x2.Customer.id][lx1.Customer.id]
                                            + Distances[y2.Customer.id][ry1.Customer.id]
                                            - Distances[x2.Customer.id][lx2.Customer.id]
                                            - Distances[y2.Customer.id][ry2.Customer.id]
                                            + Distances[x1.Customer.id][lx2.Customer.id]
                                            + Distances[y1.Customer.id][ry2.Customer.id];
                                if (delta >= 0 && !SA) continue;
                                if (MinDelta > delta){
                                    MinDelta = delta;
                                    Location1I = u1; Location1J = v1;
                                    Location2I = u2; Location2J = v2;
                                    VehicleID = j;
                                }

                            }
                        }
                    }
                }
            }
            if (MinDelta != INF){
                vector<ServeNode> helps;
                helps.clear();
                helps.shrink_to_fit();
                Vehicle &V2 = Vehicles[VehicleID];
                Vehicle tmpV1 = V1, tmpV2 = V2;
                for (int l1 = Location1I; l1 <= Location1J; l1 ++){
                    helps.push_back(V1.path[l1]);
                }
                V1.path.erase(V1.path.begin() + Location1I, V1.path.begin() + Location1J + 1);
                V1.path.insert(V1.path.begin() + Location1I, V2.path.begin() + Location2I, V2.path.begin() + Location2J + 1);
                V2.path.erase(V2.path.begin() + Location2I, V2.path.begin() + Location2J + 1);
                V2.path.insert(V2.path.begin() + Location2I, helps.begin(), helps.end());
                if (!V1.FeasibleCheck() || !V2.FeasibleCheck()){
                    V1 = tmpV1;
                    V2 = tmpV2;
                }
            }
        }
        CalcTotalPorfit();
    }

    void AddWorkers(){
        if (NotVisited.empty()) return;
        double p_H = 0.7, p_L = 0.3;
        double AvgPro = TotalPorfit / (N - NotVisited.size());
        for (int i = 0; i < Vehicles.size(); i ++){
            Vehicle &L = Vehicles[i];
            for (int j = 0; j < L.path.size(); j ++){
                double rnd = (rand() % 10000) / 10000.;
                ServeNode &S = L.path.at(j);
                if (S.NumberOfWorks > M - 1) continue;
                if (S.Customer.Revenue - CostOfWorkers[S.NumberOfWorks] > AvgPro){
                    if (rnd > p_L){
                        S.NumberOfWorks ++;
                        S.Customer.CalcServiceTime(S.NumberOfWorks);
                    }
                }else{
                    if (rnd > p_H){
                        S.NumberOfWorks ++;
                        S.Customer.CalcServiceTime(S.NumberOfWorks);
                    }
                }
            }
        }
        CalcTotalPorfit();
    }

	void print(FILE *fp){
        fprintf(fp, "TotalPorfit:%lf\nTotalRevenue:%d\nTotalDistance:%d\nTotalChargeCost:%lf\nTotalEmpolyCost:%d\n",
                    TotalPorfit, TotalRevenue, TotalDistance, TotalChargeCost, TotalEmpolyCost);
        for (int i = 0; i < Vehicles.size(); i ++){
            Vehicles[i].print(fp);
        }
        return;
	}

	bool FeasibleCheck(bool debug = 0, bool ele = 0){
        vector <bool> vis;
        vis.clear();
        vis.shrink_to_fit();
        for (int i = 0; i <= N; i ++) vis.push_back(0);
        for (int i = 0; i < Vehicles.size(); i ++){
            for (int j = 1; j < Vehicles[i].path.size()-1; j ++){
                if (vis[Vehicles[i].path.at(j).Customer.id]) {printf("%d\n", Vehicles[i].path.at(j).Customer.id); return 0;}
                vis[Vehicles[i].path.at(j).Customer.id] = 1;

            }
            bool flag = Vehicles[i].FeasibleCheck(debug, ele);
            if (!flag) return 0;
        }
        return 1;
	}

};

Solution LSsolution, SAsolution, BestSolution, TabuSolution, ALNSsolution;
