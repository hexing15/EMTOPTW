#include"base.h"
#include"io.h"
#include"Parameters.h"
#include"init_solution.h"
#include"SA.h"
#include"LS.h"
#include"ALNS.h"
#include"Tabu.h"
bool compare(const Vehicle &a, const Vehicle &b){
    return a.Porfit * b.path.size() > b.Porfit * a.path.size();
}

void OPTAlogrithm(){
    bool improve = 1;
    bad_times = 0;
    FILE *fp = fopen("best.txt", "w");
    while (improve){
      //  printf("11\n");
        SA();
      //  printf("22\n");
        LS();
   //    printf("33\n");
        ALNS();
//printf("44\n");
        TabuSearch();
        A = new int*[SolutionsPool.size()];
        Count = new int[N+1];
        for (int i = 1; i <= N; i ++) Count[i] = 0;
        for (int i = 0; i < SolutionsPool.size(); i ++){
            A[i] = new int[N+1];
            for (int j = 1; j <= N; j ++){
                A[i][j] = 0;
            }
            Vehicle L = SolutionsPool.at(i);
           // L.print(stdout);
            for (int j = 1; j < L.path.size() - 1; j ++){

                A[i][L.path.at(j).Customer.id] = 1;
                Count[L.path.at(j).Customer.id] = 1;

            }
        }

        int SizeOfPool = SolutionsPool.size();
        IloEnv env;
        try{
            IloModel model(env);
            //decision variables 决策变量
            IloNumVarArray X(env, SizeOfPool);

            for (int i = 0; i < SizeOfPool; i ++){
                X[i] = IloNumVar(env, 0, 1, ILOINT);
            }

            IloExpr expr_1(env);
            for (int i = 0; i < SizeOfPool; i ++){
                SolutionsPool[i].CalcPorfit();
                expr_1 += SolutionsPool[i].Porfit * X[i];
               // expr_12 += SolutionsPool[i].Distance * X[i];
            }
            IloObjective obj(env, expr_1, IloObjective::Maximize);
           // IloObjective obj1(env, expr_12, IloObjective::Minimize);
            model.add(obj);
        //    model.add(obj1);

            for (int i = 1; i <= N; i ++){
                if (Count[i] == 0) continue;
                IloExpr expr_2(env);

                for (int j = 0; j < SizeOfPool; j ++){
                    expr_2 += A[j][i] * X[j];
                }
                model.add(expr_2 <= 1);
            }

            IloExpr expr_3(env);
            for (int i = 0; i < SizeOfPool; i ++){
                expr_3 += X[i];
            }
            model.add(expr_3 <= K);
            IloCplex cplex(model);
            cplex.setOut(env.getNullStream());
            cplex.setWarning(env.getNullStream());
            cplex.solve();
            if (cplex.getStatus() == IloAlgorithm::Infeasible)
                env.out() << "No Solution" << endl;
            double Value = cplex.getObjValue();
            env.out() << "Solution status: " << cplex.getStatus() << endl;
            env.out() << "Solution Value " << Value << endl;
            cplex.exportModel("model.lp");
            bad_times ++;
            if (bad_times >= MAXBAD) improve = 0;

            if (Value - BestSolution.TotalPorfit > eps){
                bad_times = 0;
                improve = 1;
                int cnt = 0;
                Vehicle emptyV;
                emptyV.GenerateNewVeh();
                //printf("111\n");
                for (int i = 0; i < K; i ++) BestSolution.Vehicles[i] = emptyV;
                //printf("222\n");
                for (int i = 0; i < SizeOfPool; i ++){
                    if (cplex.getValue(X[i]) > 0.5){
                        BestSolution.Vehicles[cnt] = SolutionsPool[i];
                        cnt ++;
                    }
                }
                BestSolution.CalcTotalPorfit();
                BestSolution.buildUnvisit();
                printf("Best:\n");
                if (!BestSolution.FeasibleCheck(1)) printf("wrong\n");
                fprintf(fp, "%lf\n", BestSolution.TotalPorfit);
                BestSolution.print(stdout);
                if (LSsolution.TotalPorfit < BestSolution.TotalPorfit || LSsolution.TotalPorfit == BestSolution.TotalPorfit && LSsolution.TotalDistance < BestSolution.TotalDistance);
                    LSsolution = BestSolution;
                if (SAsolution.TotalPorfit < BestSolution.TotalPorfit|| SAsolution.TotalPorfit == BestSolution.TotalPorfit && SAsolution.TotalDistance < BestSolution.TotalDistance)
                    SAsolution = BestSolution;
                if (ALNSsolution.TotalPorfit < BestSolution.TotalPorfit|| ALNSsolution.TotalPorfit == BestSolution.TotalPorfit && ALNSsolution.TotalDistance < BestSolution.TotalDistance)
                    ALNSsolution = BestSolution;
                if (TabuSolution.TotalPorfit < BestSolution.TotalPorfit|| TabuSolution.TotalPorfit == BestSolution.TotalPorfit && TabuSolution.TotalDistance < BestSolution.TotalDistance)
                    TabuSolution = BestSolution;
            }

        }

        catch (IloException & ex){
            cerr << "Error: " << ex << endl;
        }
        catch (...){
            cerr << "Error " << endl;
        }

        env.end();


        for (int i = 0; i < SolutionsPool.size(); i ++){
            delete []A[i];
        }
        delete []A;
        if (SolutionsPool.size() > 200){
            sort(SolutionsPool.begin(), SolutionsPool.end(), compare);
            SolutionsPool.erase(SolutionsPool.begin() + 2 * K + 1, SolutionsPool.end());
        }
    }

}

bool operator <(Solution a, const Solution b){
    if (a.TotalPorfit < b.TotalPorfit || a.TotalPorfit == b.TotalPorfit && a.TotalDistance > b.TotalDistance) return 1;
    else return 0;
}

int main(int argc, char** argv){
	// ReadInstanceInfo
    int seed = unsigned (time(0));
	srand(seed);
	readCustInfo(argv[1]);
	readCostInfo(argv[2]);
	//CalcuateParameters
	SetAllParameters(argv);
	CalcALLDistance();

    //setAllParameters

    GenerateAllCustomers();
	// GenerateInitSolution
	ConstructInitSolution();
    double toc = clock();
    OPTAlogrithm();
    toc = clock() - toc;
    double TimeT = toc / CLOCKS_PER_SEC;
    if (!BestSolution.FeasibleCheck(1)) printf("wrong\n");

    printf("/***************************Final*****************/\n");
    BestSolution = max(BestSolution, max(LSsolution, max(TabuSolution, max(SAsolution, ALNSsolution))));
    BestSolution.print(stdout);
    char final_name[50];
    sprintf(final_name, "%s-%d:", argv[1], K);
    FILE *FinalResult = fopen("result.txt", "a");
    fprintf(FinalResult, "%s,%.2lf, %.3lf, %d\n", final_name, BestSolution.TotalPorfit, TimeT, seed);
    fclose(FinalResult);
}
