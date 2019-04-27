

int totalTabu = 0;

void TabuSearch(){
    Solution CurrentSolution = TabuSolution;

    CurrentSolution.Eliminator();
    if (CurrentSolution.TotalPorfit >  TabuSolution.TotalPorfit || TabuSolution.TotalPorfit == CurrentSolution.TotalPorfit && TabuSolution.TotalDistance < CurrentSolution.TotalDistance  && CurrentSolution.FeasibleCheck()) TabuSolution = CurrentSolution;
        //PutInPool(CurrentSolution);
 //   printf("22\n");
    if (bad_times > MAXBAD / 2) {
        CurrentSolution.AddWorkers();

        CurrentSolution.ZeroExchange();

        CurrentSolution.ZeroRelocate();

        CurrentSolution.OneRelocate();
        CurrentSolution.TwoRelocate();

    }

    //for (int i = 1; i <= N; i ++) TabuTable[i] = 0;
    for (int iter = 0; iter < MAXIteration; iter ++){
        totalTabu ++;
        CurrentSolution.TwoExchange(0,1,totalTabu);
        CurrentSolution.TwoOpt();

        CurrentSolution.TwoRelocate(0,1,totalTabu);

        CurrentSolution.OneExchange(0,1,totalTabu);
        CurrentSolution.OneRelocate(0,1,totalTabu);
        CurrentSolution.ZeroExchange();
        CurrentSolution.ZeroRelocate();

        bool flag = CurrentSolution.FeasibleCheck();
        if (flag){
           if (CurrentSolution.TotalPorfit > TabuSolution.TotalPorfit || TabuSolution.TotalPorfit == CurrentSolution.TotalPorfit && TabuSolution.TotalDistance < CurrentSolution.TotalDistance){
                TabuSolution = CurrentSolution;
            }
        }else CurrentSolution = TabuSolution;

        assert(CurrentSolution.FeasibleCheck(1) == 1);
    }
    PutInPool(CurrentSolution);
   // printf("44\n");
    printf("Tabu:\n");
    if (!TabuSolution.FeasibleCheck(1)) printf("wrong\n");
    TabuSolution.print(stdout);
}
