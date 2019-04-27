


void LS(){
  //  printf("11\n");
    Solution CurrentSolution = LSsolution;

        CurrentSolution.Eliminator();
        if (CurrentSolution.TotalPorfit > LSsolution.TotalPorfit || LSsolution.TotalPorfit == CurrentSolution.TotalPorfit && LSsolution.TotalDistance < CurrentSolution.TotalDistance&& CurrentSolution.FeasibleCheck()) LSsolution = CurrentSolution;
        //PutInPool(CurrentSolution);
 //   printf("22\n");
    if (bad_times > MAXBAD / 2) {
        //printf("11\n");
        CurrentSolution.AddWorkers();
        //printf("22\n");
        CurrentSolution.ZeroExchange();
        //printf("33\n");
        CurrentSolution.ZeroRelocate();

        CurrentSolution.OneRelocate();
        CurrentSolution.TwoRelocate();
        //printf("44\n");
      //  PutInPool(CurrentSolution);
        //printf("55\n");
    }
   // printf("33\n");
    for (int iter = 0; iter < MAXIteration; iter ++){
        CurrentSolution.TwoExchange();
        CurrentSolution.TwoOpt();
        CurrentSolution.TwoRelocate();
        CurrentSolution.OneExchange();
        CurrentSolution.OneRelocate();
        CurrentSolution.ZeroExchange();
        CurrentSolution.ZeroRelocate();
        bool flag = CurrentSolution.FeasibleCheck();
        if (flag){
            if (CurrentSolution.TotalPorfit > LSsolution.TotalPorfit || LSsolution.TotalPorfit == CurrentSolution.TotalPorfit && LSsolution.TotalDistance < CurrentSolution.TotalDistance){
                LSsolution = CurrentSolution;
            }
        }else CurrentSolution = LSsolution;

        assert(CurrentSolution.FeasibleCheck(1) == 1);

    }
    PutInPool(CurrentSolution);
   // printf("44\n");
    printf("LS:\n");
    if (!LSsolution.FeasibleCheck(1)) printf("wrong\n");
    LSsolution.print(stdout);
   // FILE *LSfp = fopen("LSfp.txt", "a");
   // fprintf(LSfp, "%.2lf\n", LSsolution.TotalPorfit);
   // fclose(LSfp);
}
