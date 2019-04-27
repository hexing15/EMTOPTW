
Solution CurrentSolution, NeighbourSolution;

double *OperatorWeight;
int TotalOperator = 9;
double lambda = 0.95;
const double UpdateWeight[] = {1., 0.75, 0.5, 0.25};
typedef pair<int, double> pair_int_double;
pair_int_double OperatorChoose(int OperatorID){
    pair_int_double ans = make_pair(0, 0.);
    double totalProb = 0;
    double rnd = (rand() % 10000) / 10000.;
    for (int i = 0; i < TotalOperator; i ++){
        totalProb += OperatorWeight[i];
    }
    totalProb *= rnd;
    double curProb = 0;
    for (int i = 0; i < TotalOperator; i ++){
        curProb += OperatorWeight[i];
        if (curProb >= totalProb) {
            OperatorID = i;
            break;
        }
    }
    switch (OperatorID){
        case 0:
            NeighbourSolution.Eliminator();
            break;
        case 1:
            NeighbourSolution.ZeroRelocate();
            break;
        case 2:
            NeighbourSolution.AddWorkers();
            NeighbourSolution.ZeroExchange();
            NeighbourSolution.ZeroRelocate();
            NeighbourSolution.OneRelocate();
            NeighbourSolution.TwoRelocate();
            break;
        case 3:
            NeighbourSolution.ZeroExchange();
            break;
        case 4:
            NeighbourSolution.OneExchange();
            NeighbourSolution.ZeroExchange();
            NeighbourSolution.ZeroRelocate();
            break;
        case 5:
            NeighbourSolution.OneRelocate();
            NeighbourSolution.ZeroExchange();
            NeighbourSolution.ZeroRelocate();
            break;
        case 6:
            NeighbourSolution.TwoRelocate();
            NeighbourSolution.ZeroExchange();
            NeighbourSolution.ZeroRelocate();
            break;
        case 7:
            NeighbourSolution.TwoExchange();
            NeighbourSolution.ZeroExchange();
            NeighbourSolution.ZeroRelocate();
            break;
        case 8:
            NeighbourSolution.TwoOpt();
            NeighbourSolution.ZeroRelocate();
            NeighbourSolution.ZeroExchange();
            break;
        default:break;
    }
    ans.first = NeighbourSolution.TotalDistance;
    ans.second = NeighbourSolution.TotalPorfit;
}

void ALNS(){
    OperatorWeight = new double[TotalOperator];
    for (int i = 0; i < TotalOperator; i ++){
        OperatorWeight[i] = 0.5;
    }

    CurrentSolution = ALNSsolution;
    NeighbourSolution = ALNSsolution;

    for (int iter = 0; iter < MAXIteration * 3; iter ++){
        int OperatorID = 0;
        double operdelta = 0;
        pair_int_double NewSolutionValue = OperatorChoose(OperatorID);
        if (!NeighbourSolution.FeasibleCheck()) continue;
        if (NewSolutionValue.second > ALNSsolution.TotalPorfit
            || (NewSolutionValue.second == ALNSsolution.TotalPorfit && NewSolutionValue.first < ALNSsolution.TotalDistance)){
            ALNSsolution = NeighbourSolution;
            CurrentSolution = NeighbourSolution;
            operdelta = UpdateWeight[0];
        }else{
            if (NewSolutionValue.second > CurrentSolution.TotalPorfit){
                CurrentSolution = NeighbourSolution;
                operdelta = UpdateWeight[1];
            }else{
                if (NewSolutionValue.first < CurrentSolution.TotalDistance){
                    CurrentSolution = NeighbourSolution;
                    operdelta = UpdateWeight[2];
                }else{
                    NeighbourSolution = CurrentSolution;
                    operdelta = UpdateWeight[3];
                }
            }
        }
        OperatorWeight[OperatorID] = lambda * OperatorWeight[OperatorID] + (1 - lambda) * operdelta;
    }
    PutInPool(CurrentSolution);
    if (ALNSsolution.TotalPorfit < CurrentSolution.TotalPorfit || (ALNSsolution.TotalPorfit == CurrentSolution.TotalPorfit && ALNSsolution.TotalDistance > CurrentSolution.TotalDistance))
        ALNSsolution = CurrentSolution;
    printf("ALNS:\n");
    if (!ALNSsolution.FeasibleCheck(1)) printf("wrong\n");
    ALNSsolution.print(stdout);
    delete []OperatorWeight;
}
