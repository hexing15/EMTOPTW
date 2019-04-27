#include"base.h"

ILOSTLBEGIN

void readCustInfo(char *FileName){
	FILE *fp = fopen(FileName, "r");
	fscanf(fp, "%d %d", &N, &M);
	Customers = new Vertex[N + 2];
	for (int i = 0; i <= N; i ++){
		Customers[i].ReadInfo(fp);
	}
	Customers[N + 1] = Customers[0];
	fclose(fp);
	return;
}

void readCostInfo(char *FileName){
	FILE *fp = fopen(FileName, "r");
	CostOfWorkers = new int[M + 1];
	CostOfWorkers[0] = 0;
	for (int i = 1; i <= M; i ++){
		fscanf(fp, "%d", &CostOfWorkers[i]);
	}
	fclose(fp);
	return;
}


void CalcALLDistance(){
	Distances = new int*[N + 1];
	TimeCost = new int*[N + 1];
	for (int i = 0; i <= N; i ++){
		Distances[i] = new int[N + 1];
		TimeCost[i] = new int[N + 1];
		for (int j = 0; j <= N; j ++){
			Distances[i][j] = CalcDistance(Customers[i], Customers[j]);
			TimeCost[i][j] = (int)((double)Distances[i][j] / (speed * 1.0) + 0.5);
		}
	}
	return;
}


void slove(){
	IloEnv env;
	try{
        char str[15];
        IloModel model(env);
        //set decision variable
        IloNumVarArray2 y(env, V), t(env, V), f(env, V), q(env, V), h(env, V);
        for (int i = 0; i < V; i ++){
            y[i] = IloNumVarArray(env, K);
            t[i] = IloNumVarArray(env, K);
            f[i] = IloNumVarArray(env, K);
            q[i] = IloNumVarArray(env, M + 1);
            h[i] = IloNumVarArray(env, K);
            for (int j = 0; j < K; j ++){
                string tmp, y_str, t_str, f_str, h_str;
                sprintf(str, "_{%d,%d}", i, j);
                tmp = str;
                y_str = "y" + tmp;
                t_str = "t" + tmp;
                f_str = "f" + tmp;
                h_str = "h" + tmp;
                y[i][j] = IloNumVar(env, 0, 1, ILOINT, y_str.c_str());
                t[i][j] = IloNumVar(env, Customers[0].TimeWindowL, Customers[0].TimeWindowR,  t_str.c_str());
                f[i][j] = IloNumVar(env, 0, INF,  f_str.c_str());
                h[i][j] = IloNumVar(env, 0, INF,  h_str.c_str());
            }
            for (int j = 1; j <= M; j ++){
                string tmp, q_str;
                sprintf(str, "_{%d,%d}", i, j);
                tmp = str;
                q_str = "q" + tmp;
                q[i][j] = IloNumVar(env, 0, 1, ILOINT, q_str.c_str());
            }
        }

        IloNumVarArray3 x(env, V);
        for (int i = 0; i < V; i ++){
            x[i] = IloNumVarArray2(env, V);
            for (int j = 0; j < V; j ++){
                x[i][j] = IloNumVarArray(env, K);
                for (int k = 0; k < K; k ++){
                    sprintf(str, "x_{%d,%d,%d}", i, j, k);
                    x[i][j][k] = IloNumVar(env, 0, 1, ILOINT, str);
                }
            }
        }
        // obj
        IloExpr expr_1(env);
        for (int i = 1; i <= N; i ++){
            for (int k = 0; k < K; k ++){
                expr_1 += Customers[i].Revenue * y[i][k];
                expr_1 -= UnitChargeCost * h[i][k];
            }
            for (int l = 1; l <= M; l ++){
                expr_1 -= CostOfWorkers[l] * q[i][l];
            }
        }
        IloObjective obj(env, expr_1, IloObjective::Maximize);
        model.add(obj);

        //eq(3.2)
        for (int j = 1; j < V; j ++){
            for (int k = 0; k < K; k ++){
                IloExpr expr_2(env);
                for (int i = 0; i <= N; i ++){
                    if (i == j) continue;
                    expr_2 += x[i][j][k];
                }
                model.add(expr_2 == y[j][k]);
            }
        }

        //eq(3.3)
        for (int i = 0; i <= N; i ++){
            for (int k = 0; k < K; k ++){
                IloExpr expr_3(env);
                for (int j = 1; j < V; j ++){
                    if (i == j) continue;
                    expr_3 += x[i][j][k];
                }
                model.add(expr_3 == y[i][k]);
            }
        }
        //eq(3.4)
        for (int k = 0; k < K; k ++){
            IloExpr expr_4(env);
            for (int j = 1; j < V; j ++){
                expr_4 += x[0][j][k];
            }
            model.add(expr_4 == 1);
        }

        //eq(3.5)
        for (int k = 0; k < K; k ++){
            IloExpr expr_5(env);
            for (int i = 0; i <= N; i ++){
                expr_5 += x[i][N + 1][k];
            }
            model.add(expr_5 == 1);
        }

		//eq(3.6)
		for (int k = 0; k < K; k ++){
			IloExpr expr_6(env);
			for (int i = 1; i <= N; i ++){
				expr_6 += Customers[i].Demand * y[i][k];
			}
			model.add(expr_6 <= Capicity);
		}

		//eq(3.7)
		for (int i = 1; i <= N; i ++){
			IloExpr expr_7(env);
			for (int l = 1; l <= M; l ++){
				expr_7 += q[i][l];
			}
			for (int k = 0; k < K; k ++){
                expr_7 -= y[i][k];
 			}
			model.add(expr_7 == 0);
		}

		//eq(3.8)
		for (int i = 1; i <= N; i ++){
			IloExpr expr_8(env);
			for (int k = 0; k < K; k ++){
				expr_8 += y[i][k];
			}
			model.add(expr_8 <= 1);
		}

		//eq(3.9)
		for (int k = 0; k < K; k ++){
			for (int i = 1; i <= N; i ++){
				IloExpr expr_9(env);
				expr_9 = t[i][k] + INF * (1 - y[i][k]);
				model.add(expr_9 >= Customers[i].TimeWindowL);
			}
		}


		//eq(3.10)
		for (int k = 0; k < K; k ++){
			for (int i = 1; i <= N; i ++){
				IloExpr expr_10(env);
				for (int l = 1; l <= M; l ++){
					expr_10 += q[i][l] * Customers[i].CalcServiceTime(l);
				}
				expr_10 += t[i][k] - INF * (1 - y[i][k]);
				model.add(expr_10 <= Customers[i].TimeWindowR);
			}
		}
		//eq(3.11)
		IloNumVarArray2 g(env, V);
		for (int i = 0; i < V; i ++){
			g[i] = IloNumVarArray(env, K);
			for (int k = 0; k < K; k ++){
				sprintf(str, "g_{%d,%d}", i, k);
				g[i][k] = IloNumVar(env, 0, INF, ILOINT, str);
			}
		}
		for (int i = 0; i <= N; i ++){
			for (int j = 1; j < V; j ++){
                if (i == j) continue;
				for (int k = 0; k < K; k ++){
					int j_1 = j;
					if (j_1 == N + 1) j_1 = 0;
					IloExpr expr_11(env);
					expr_11 += t[i][k];
					expr_11 += TimeCost[i][j_1];
					expr_11 -= INF * (1 - x[i][j][k]);
					expr_11 += g[i][k];
					model.add(t[j][k] >= expr_11);
				}
			}
		}
		for (int i = 0; i <= N; i ++){
			for (int k = 0; k < K; k ++){
				IloExpr expr_11(env);
				for (int l = 1; l <= M; l ++){
					expr_11 += q[i][l] * Customers[i].CalcServiceTime(l);
				}
				model.add(expr_11 <= g[i][k]);
				model.add(h[i][k] <= g[i][k]);
			}
		}
		//eq (3.12)
		for (int k = 0; k < K; k ++){
			IloExpr expr_12(env);
			expr_12 = f[0][k] - H;
			model.add(expr_12 == 0);
		}
		//eq(3.13)
		//IloNumVarArray2 G(env, V);
		/*for (int i = 0; i < V; i ++){
			G[i] = IloNumVarArray(env, K);
			for (int k = 0; k < K; k ++){
				sprintf(str, "G_{%d,%d}", i, k);
				G[i][k] = IloNumVar(env, 0, INF, ILOINT, str);
			}
		}*/
		for (int i = 0; i <= N; i ++){
			for (int j = 1; j < V; j ++){
				for (int k = 0; k < K; k ++){
                    int j_1 = j;
                    if (j_1 == N + 1) j_1 = 0;
                    if (i == j) continue;
                    IloExpr expr_13(env);
                    expr_13 += ChargeSpeed * h[i][k] + f[i][k];
                    expr_13 -= Distances[i][j_1];
                    expr_13 += INF * (1 - x[i][j][k]);
					model.add(f[j][k] <= expr_13);
					model.add(f[j][k] >= f[i][k] - Distances[i][j_1] - INF * (1 - x[i][j][k]));
				}
			}
		}
		for (int i = 0; i < V; i ++){
			for (int k = 0; k < K; k ++){
			//	model.add(G[i][k] <= ChargeSpeed * h[i][k] + f[i][k]);
				model.add(ChargeSpeed * h[i][k] + f[i][k] <= H);
				//model.add(f[i][k] <= H);
			}
		}

		IloCplex cplex(model);
		cplex.exportModel("model.lp");
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.solve();
		if (cplex.getStatus() == IloAlgorithm::Infeasible)
			env.out() << "No Solution" << endl;

		env.out() << "Solution status: " << cplex.getStatus() << endl;
		env.out() << "Solution Value " << cplex.getObjValue() << endl;
        printf("q:\n");
		for (int i = 1; i <= N; i ++){
            for (int l = 1; l <= M; l ++){
                if (cplex.getValue(q[i][l]) > 0.5)
                    cout << 1 << " ";
                else
                    cout << 0 << " ";
            }
            printf("\n");
		}

		for (int k = 0; k < K; k ++){
            printf("\ny\n");
            for (int i = 0; i < V; i ++){
                if (cplex.getValue(y[i][k]) > 0.5)
                    cout << 1 << " ";
                else
                    cout << 0 << " ";
            }

            printf("\nh\n");
            for (int i = 0; i < V; i ++){
                if (cplex.getValue(y[i][k]) > 0.5)
                 cout << cplex.getValue(h[i][k]) << " ";
            }

            printf("\nf\n");
            for (int i = 0; i < V; i ++){
                if (cplex.getValue(y[i][k]) > 0.5)
                 cout << cplex.getValue(f[i][k]) << " ";
            }


            printf("\npath\n");
            int X = 0;
            cout << X << " ";
            while (X != N + 1){
                for (int i = 1; i < V; i ++){
                    if (i == X) continue;
                    if (cplex.getValue(x[X][i][k])>0.5) {
                        cout << i << " ";
                        X = i;
                        break;
                    }
                }
            }
            printf("\n");
		}

	}

	catch (IloException & ex){
		cerr << "Error: " << ex << endl;
	}
	catch (...){
		cerr << "Error " << endl;
	}

	env.end();
}

int main(int argc, char *argv[]){
    //Read Data
    readCustInfo(argv[1]);
    readCostInfo(argv[2]);
    CalcALLDistance();
    V = N + 2;

	//Creating the environment
	slove();
}
