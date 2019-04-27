
void CalcALLDistance(){
	Distances = new int*[N + 1];
	TimeCost = new int*[N + 1];
	for (int i = 0; i <= N; i ++){
		Distances[i] = new int[N + 1];
		TimeCost[i] = new int[N + 1];
		for (int j = 0; j <= N; j ++){
			Distances[i][j] = CalcDistance(Customers[i], Customers[j]);
			TimeCost[i][j] = (int)((Distances[i][j] / (speed * 1.0)));
		}
	}
	return;
}

void SetAllParameters(char** str){
	UnitChargeCost = atof(str[3]);
    Capicity = atoi(str[4]);
    K = atoi(str[5]);
    speed = atoi(str[6]);
    H = atoi(str[7]);
    ChargeSpeed = atoi(str[8]);
    MAXIteration = 40;
    MAXBAD = 20;
    TabuTable = new int[N + 1];
    for (int i = 0; i <= N; i ++){
        TabuTable[i] = 0;
    }
    return;
}

void GenerateAllCustomers(){
    for (int i = 0; i <= N + 1; i ++){
        ServeNode CurrentCustomer;
        CurrentCustomer.Customer = Customers[i];
        CurrentCustomer.EF = CurrentCustomer.LS = 0;
        CurrentCustomer.NumberOfWorks = CurrentCustomer.TimeOfCharge = 0;
        AllCustomers.push_back(CurrentCustomer);
    }
}
