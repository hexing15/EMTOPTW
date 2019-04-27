
double Temp;
double alpha = 0.8;
double prob;

void SA(){
//FILE *SAfp = fopen("sa.txt", "a");
    Solution CurrentSolution, NeighborSolution;
    CurrentSolution = SAsolution;
    NeighborSolution = CurrentSolution;
    Temp = CurrentSolution.TotalDistance * 0.1;

    NeighborSolution.Eliminator();
    if (NeighborSolution.TotalPorfit > CurrentSolution.TotalPorfit || (NeighborSolution.TotalPorfit == CurrentSolution.TotalPorfit && NeighborSolution.TotalDistance < CurrentSolution.TotalDistance)
            && NeighborSolution.FeasibleCheck()) CurrentSolution = NeighborSolution;
    if (CurrentSolution.TotalPorfit > SAsolution.TotalPorfit
            || (CurrentSolution.TotalDistance < SAsolution.TotalDistance && CurrentSolution.TotalPorfit == SAsolution.TotalPorfit))SAsolution = CurrentSolution;
    //    PutInPool(CurrentSolution);

    if (bad_times > MAXBAD / 2) {
        NeighborSolution.AddWorkers();
        NeighborSolution.ZeroExchange();
        NeighborSolution.ZeroRelocate();
        NeighborSolution.OneRelocate();
        NeighborSolution.TwoRelocate();
        if (NeighborSolution.TotalPorfit > CurrentSolution.TotalPorfit ||  (NeighborSolution.TotalPorfit == CurrentSolution.TotalPorfit && NeighborSolution.TotalDistance < CurrentSolution.TotalDistance) && NeighborSolution.FeasibleCheck()) CurrentSolution = NeighborSolution;
        if (CurrentSolution.TotalPorfit > SAsolution.TotalPorfit) SAsolution = CurrentSolution;
     //   PutInPool(CurrentSolution);
    }
    for (int iter = 0; iter < MAXIteration; iter ++){
        NeighborSolution.TwoExchange(1,0,0);
        NeighborSolution.TwoOpt(1);
        NeighborSolution.TwoRelocate(1,0,0);

        prob = - NeighborSolution.TotalDistance + CurrentSolution.TotalDistance;
        prob = exp(prob / (1.0 * Temp));
        double rnd = (rand() % 10000) / 10000.;

        if (prob < rnd  || !NeighborSolution.FeasibleCheck()){
            NeighborSolution = CurrentSolution;
        }else CurrentSolution = NeighborSolution;

        NeighborSolution.OneExchange(1,0,0);

        NeighborSolution.OneRelocate(1,0,0);

        prob = - NeighborSolution.TotalDistance + CurrentSolution.TotalDistance;
        prob = exp(prob / (1.0 * Temp));
        rnd = (rand() % 10000) / 10000.;
        if (prob < rnd || NeighborSolution.TotalPorfit < CurrentSolution.TotalPorfit
                || (NeighborSolution.TotalPorfit == CurrentSolution.TotalPorfit && NeighborSolution.TotalDistance >  CurrentSolution.TotalDistance) || !NeighborSolution.FeasibleCheck()){
            NeighborSolution = CurrentSolution;
        }else CurrentSolution = NeighborSolution;

        NeighborSolution.ZeroExchange();
        NeighborSolution.ZeroRelocate();
        if (NeighborSolution.TotalPorfit > CurrentSolution.TotalPorfit
        || (NeighborSolution.TotalPorfit == CurrentSolution.TotalPorfit && NeighborSolution.TotalDistance < CurrentSolution.TotalDistance)
        && NeighborSolution.FeasibleCheck()) CurrentSolution = NeighborSolution;

        if (CurrentSolution.TotalPorfit > SAsolution.TotalPorfit || (CurrentSolution.TotalDistance == SAsolution.TotalDistance && CurrentSolution.TotalPorfit < SAsolution.TotalPorfit)) SAsolution = CurrentSolution;
        assert(CurrentSolution.FeasibleCheck(1) == 1);
        Temp *= alpha;
    }
    PutInPool(CurrentSolution);
    printf("SA:\n");
    if (!SAsolution.FeasibleCheck(1)) printf("wrong\n");
    SAsolution.print(stdout);
   // fprintf(SAfp, "%.2lf\n", SAsolution.TotalPorfit);
    //fclose(SAfp);
}
