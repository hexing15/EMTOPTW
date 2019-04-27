void PutInPool(const Vehicle &L){
    if (L.path.size() != 2)
        SolutionsPool.push_back(L);
}

void PutInPool(const Solution &L){
    for (int i = 0; i < L.Vehicles.size(); i ++){
        if (L.Vehicles[i].path.size() != 2)
            SolutionsPool.push_back(L.Vehicles.at(i));
    }
}


bool comp(const ServeNode &a, const ServeNode &b){
	return (a.Customer.Revenue * (1.0 * Distances[0][b.Customer.id])) >
	 	   (b.Customer.Revenue * (Distances[0][a.Customer.id] * 1.0));
}

bool comp1(const ServeNode &a, const ServeNode &b){
	return a.Customer.id < b.Customer.id;
}

void ConstructInitSolution(){
	sort(AllCustomers.begin() + 1, AllCustomers.end() - 1, comp);
	SolutionsPool.clear();
	SolutionsPool.shrink_to_fit();

	printf("LS Init:\n");
	LSsolution.GenerateInitSolution();
	LSsolution.CalcTotalPorfit();
	LSsolution.print(stdout);
	assert(LSsolution.FeasibleCheck(1) == 1);
    PutInPool(LSsolution);

	printf("SA Init:\n");
	SAsolution.GenerateInitSolution();
	SAsolution.CalcTotalPorfit();
	SAsolution.print(stdout);
	assert(SAsolution.FeasibleCheck(1) == 1);
    PutInPool(SAsolution);

    printf("ALNS Init:\n");
    ALNSsolution.GenerateInitSolution();
    ALNSsolution.CalcTotalPorfit();
    ALNSsolution.print(stdout);
    assert(ALNSsolution.FeasibleCheck(1) == 1);
    PutInPool(ALNSsolution);

    printf("Tabu Init:\n");
    TabuSolution.GenerateInitSolution();
    TabuSolution.CalcTotalPorfit();
    TabuSolution.print(stdout);
    assert(TabuSolution.FeasibleCheck(1) == 1);
    PutInPool(TabuSolution);

    if (LSsolution.TotalPorfit > SAsolution.TotalPorfit) BestSolution = LSsolution;
    else BestSolution = SAsolution;
    if (BestSolution.TotalPorfit < ALNSsolution.TotalPorfit) BestSolution = ALNSsolution;
    if (BestSolution.TotalPorfit < TabuSolution.TotalPorfit) BestSolution = TabuSolution;
    sort(AllCustomers.begin()+1, AllCustomers.end()-1, comp1);
    return;
}
