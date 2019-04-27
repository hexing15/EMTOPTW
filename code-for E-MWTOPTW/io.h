
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
