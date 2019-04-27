#include<ilcplex/ilocplex.h>
#include<bits/stdc++.h>
#define sqr(x) ((x)*(x))
using namespace std;
typedef IloArray <IloNumVarArray> IloNumVarArray2;
typedef IloArray <IloNumVarArray2> IloNumVarArray3;
const int INF = 0x3f3f3f3f;
int **Distances;// Distance_from_to
int **TimeCost; // time
int *CostOfWorkers;
int N;// Number of customer
int M; // Maximum NumberOfWorks
int V;
double UnitChargeCost = 0.1; // 单位时间充电成本
int Capicity = 200;
int K = 1;
int speed = 250;
int H = 99999999;
double ChargeSpeed = 1000;
class Vertex{
public:
	int id;
	int Coord_x, Coord_y;
	int TimeWindowL, TimeWindowR;
	int Demand, Revenue;
	int serviceTime; //这个与需求和雇佣人数相关


	Vertex(int i, int x, int y, int TL, int TR, int D, int R){
		id = i;
		Coord_x = x; Coord_y = y;
		TimeWindowL = TL; TimeWindowR = TR;
		Demand = D; Revenue = R;
	}

	void ReadInfo(FILE *fp){
		fscanf(fp, "%d %d %d %d %d %d %d",&id, &Coord_x, &Coord_y, &TimeWindowL, &TimeWindowR, &Demand, &Revenue);
	}

	int CalcServiceTime(int NumberOfWorks){
		double decent = pow((double)NumberOfWorks, 0.9);
		return serviceTime = (int)((double)Demand / decent + 0.5);
	}// serviceTime = Demand / (N ^ 0.9)

	friend int CalcDistance(const Vertex &a, const Vertex &b);

	Vertex(){}

	~Vertex();
};
inline int CalcDistance(const Vertex &a, const Vertex &b){
	return (int)(sqrt(sqr(a.Coord_x - b.Coord_x) + sqr(a.Coord_y - b.Coord_y)) + 0.5);
}
Vertex *Customers;
