#include "Trivium.h"


void FilterTermBDPT(MapPoly& RetSuperpoly, Poly& ExpandingANF, int TargetRound, Ivterm &CubeFlag, int OutPutRound, string files, string appfiles)
{
	const int dim = CubeFlag.count();
	int State[288] = { 0 };
	int SState[288] = { 0 };

	int TempRound = OutPutRound - TargetRound;
	cout << "Round = " << TempRound << endl;
	Numeric_Mapping(State, SState, CubeFlag, TempRound);


	Poly TempANF(ExpandingANF.Size);
	for (uint32_t pt = 0; pt < ExpandingANF.Size; pt++)
	{
		int TestDeg = 0;
		for (int pd = 0; pd < 288; pd++)
			if ((ExpandingANF.poly[pt].pterm[(pd >> 5)] >> (31 - (pd & 0x1f))) & 1)
			{
				if ((pd != 287) && ((ExpandingANF.poly[pt].pterm[((pd + 1) >> 5)] >> (31 - ((pd + 1) & 0x1f))) & 1))
					TestDeg += SState[pd++];
				else
					TestDeg += State[pd];
			}
		if (TestDeg >= dim)
			TempANF.poly[TempANF.Size++] = ExpandingANF.poly[pt];
	}
	cout << "Save Terms(Pass Numeric Mapping): " << TempANF.Size << endl;
	ExpandingANF.SetPolyLen(TempANF.Size);


	int AllThread = std::thread::hardware_concurrency();
	cout << "The total number of threads available to this computer: " << AllThread << endl;
	omp_set_num_threads(AllThread);
	omp_set_nested(1); 
	int ThreadNumForLoop = AllThread >> 3;
	int count = 0;


#pragma omp parallel for schedule(dynamic,10) num_threads(ThreadNumForLoop)   // num_threads(5)
	for (int64_t pd = 0; pd < TempANF.Size; pd++)
	{
		clock_t ModelStartTime = clock();

		MapPoly TempPoly;
		bool TimeLimteFalg = TriviumEvalBDPT(TempPoly, TempANF.poly[pd], CubeFlag, TempRound, true, 360);
#pragma omp critical
		{
			if (TimeLimteFalg)
				ExpandingANF.poly[ExpandingANF.Size++] = TempANF.poly[pd];
			else
				RetSuperpoly.Merge(TempPoly);
		}
		clock_t ModelEndTime = clock();

#pragma omp atomic
		count += 1;

		cout << "Round: " << TempRound << ", ThreadID: " << omp_get_thread_num() << ", Have Tested Terms: (" << pd + 1 << ", "
			<< count << " / " << TempANF.Size << "), Superpoly Size: " << RetSuperpoly.size() << ", Accepted Terms: "
			<< ExpandingANF.Size << ", RunningTime: " << double(ModelEndTime - ModelStartTime) / CLOCKS_PER_SEC << 's' << endl;
	}


	cout << "Save Terms(Pass CBDP_MILP): " << ExpandingANF.Size << endl;

	char OutFile[100], OutFile2[100], AppOutFile[100];

	files.copy(OutFile, files.length(), 0);
	*(OutFile + files.length()) = '\0';

	files = files.substr(0, files.length() - 4) + '(' + to_string(TempRound) + ").txt";
	files.copy(OutFile2, files.length(), 0);
	*(OutFile2 + files.length()) = '\0';

	appfiles.copy(AppOutFile, appfiles.length(), 0);
	*(AppOutFile + appfiles.length()) = '\0';


	ExpandingANF.WriteValueToFile(TempRound, OutFile, "Terms: " + to_string(ExpandingANF.Size), true);
	RetSuperpoly.WritePolyValueToFile(OutFile, "Recovered Superpoly", true);
	RetSuperpoly.WritePolyToFile(OutFile, "Recovered Superpoly");

	ExpandingANF.WriteValueToFile(TempRound, OutFile2, "Terms: " + to_string(ExpandingANF.Size), true);
	RetSuperpoly.WritePolyValueToFile(OutFile2, "Recovered Superpoly", true);
	RetSuperpoly.WritePolyToFile(OutFile2, "Recovered Superpoly");


	ExpandingANF.WriteOutputToFile(TempRound, AppOutFile, "Terms: " + to_string(ExpandingANF.Size));
	ExpandingANF.WriteValueToFile(TempRound, AppOutFile, "Terms: " + to_string(ExpandingANF.Size));
	RetSuperpoly.WritePolyValueToFile(AppOutFile, "Recovered Superpoly", true);
	RetSuperpoly.WritePolyToFile(AppOutFile, "Recovered Superpoly");
}



bool FilterPartTermBDPT(MapPoly& RetSuperpoly, int& InPutRound, Poly& ANF, Ivterm &CubeFlag, string files, string appfiles)
{
	const int dim = CubeFlag.count();
	Poly InPutANF, OutPutANF;
	Poly UpdateFunc[8], TempPoly1(SIZEBw), TempPoly2(SIZEBw), TempPoly3(SIZEBw);
	Poly SNull(0), S1(4), S94(4), S178(4), S1S94(16), S1S178(16), S94S178(16), SAll(64);
	InvExpressInit(TempPoly1, S1, S94, S178, S1S94, S1S178, S94S178, SAll);
	UpdateFunc[0] = SNull;
	UpdateFunc[1] = S178;
	UpdateFunc[2] = S94;
	UpdateFunc[3] = S94S178;
	UpdateFunc[4] = S1;
	UpdateFunc[5] = S1S178;
	UpdateFunc[6] = S1S94;
	UpdateFunc[7] = SAll;

	InPutANF.PolyCopy(ANF);
	uint32_t SaveTermNum = 1;

	if (InPutANF.Size == 0)
		SaveTermNum = 0;

	int AllThread = std::thread::hardware_concurrency();
	cout << "The total number of threads available to this computer: " << AllThread << endl;
	omp_set_num_threads(AllThread);
	omp_set_nested(1);
	system("mkdir RoundLog");

	int TimeLimit = 60, TimeLimitQ = 60, MaxDeg = 80;
	int ThreadNumForLoop = AllThread >> 2;
	Poly *ThreadTempPoly = new Poly[AllThread];
	MapPoly *ThreadTempSP = new MapPoly[AllThread];
	const int LEN = 9;
	const int TRound[LEN] = { 20, 50, 100, 150, 200, 250, 300, 350, 400 };
	const int Ttime[LEN] = { 360, 300, 270, 240, 210, 180, 150, 120, 90 };
	const int TimeLimitRatio = 1;
	const int ParallelNum = 8;
	const double ThresholdCoeff = 0.50;
	const int BackExpandNumThr = 10000;
	const int RecoveryNumThr = 5000;


	while ((SaveTermNum <= 100000) && (SaveTermNum))
	{
		int TESTFLAG = true;
		clock_t TimeOutStart = clock();

		while (((OutPutANF.Size <= RecoveryNumThr)||(TESTFLAG)) && (InPutRound) )
		{
			TempPoly1.Size = 0;
			TempPoly2.Size = 0;
			TempPoly3.Size = 0;


			for (uint32_t pt = 0; pt < InPutANF.Size; pt++)
			{
				uint8_t flag = (((InPutANF.poly[pt].pterm[0] >> 31) & 1) << 2) | (((InPutANF.poly[pt].pterm[2] >> 2) & 1) << 1) | ((InPutANF.poly[pt].pterm[5] >> 14) & 1);
				uint8_t flagweight = ((InPutANF.poly[pt].pterm[0] >> 31) & 1) + ((InPutANF.poly[pt].pterm[2] >> 2) & 1) + ((InPutANF.poly[pt].pterm[5] >> 14) & 1);
				for (int i = 0; i < 8; i++)
					InPutANF.poly[pt].pterm[i] = (InPutANF.poly[pt].pterm[i] << 1) | ((InPutANF.poly[pt].pterm[i + 1] >> 31) & 1);
				InPutANF.poly[pt].pterm[8] <<= 1;
				InPutANF.poly[pt].deg -= flagweight;
				if (flag)
				{
					InPutANF.poly[pt].pterm[2] &= 0xfffffff7;
					InPutANF.poly[pt].pterm[5] &= 0xffff7fff;
					PolyMul(TempPoly2, UpdateFunc[flag], InPutANF.poly[pt]);
				}
				else
					TempPoly1.poly[TempPoly1.Size++] = InPutANF.poly[pt];
			}

			for (uint32_t pt = 0; pt < OutPutANF.Size; pt++)
			{
				uint8_t flag = (((OutPutANF.poly[pt].pterm[0] >> 31) & 1) << 2) | (((OutPutANF.poly[pt].pterm[2] >> 2) & 1) << 1) | ((OutPutANF.poly[pt].pterm[5] >> 14) & 1);
				uint8_t flagweight = ((OutPutANF.poly[pt].pterm[0] >> 31) & 1) + ((OutPutANF.poly[pt].pterm[2] >> 2) & 1) + ((OutPutANF.poly[pt].pterm[5] >> 14) & 1);
				for (int i = 0; i < 8; i++)
					OutPutANF.poly[pt].pterm[i] = (OutPutANF.poly[pt].pterm[i] << 1) | ((OutPutANF.poly[pt].pterm[i + 1] >> 31) & 1);
				OutPutANF.poly[pt].pterm[8] <<= 1;
				OutPutANF.poly[pt].deg -= flagweight;
				if (flag)
				{
					OutPutANF.poly[pt].pterm[2] &= 0xfffffff7;
					OutPutANF.poly[pt].pterm[5] &= 0xffff7fff;
					PolyMul(TempPoly2, UpdateFunc[flag], OutPutANF.poly[pt]);
				}
				else
					TempPoly2.poly[TempPoly2.Size++] = OutPutANF.poly[pt];
			}
			InPutRound -= 1;
			OutPutANF.Size = 0;
			InPutANF.PolyCopy(TempPoly1);
			cout << endl << "Before Degest: OutPutANF Terms: " << TempPoly2.Size << ", SaveTerms: " << InPutANF.Size << ", InPutRound = " << InPutRound << endl;
			// ---------------------------------------------------------------
			TempPoly2.RemoveDup();
			cout << "DegEst1(Pass RemoveDup): OutPutANF Terms: " << TempPoly2.Size << ", SaveTerms: " << InPutANF.Size << ", InPutRound = " << InPutRound << endl;
			// ---------------------------------------------------------------


			int State[288] = { 0 };
			int SState[288] = { 0 };
			Numeric_Mapping(State, SState, CubeFlag, InPutRound);
			uint32_t pf = 0;
			for (uint32_t pt = 0; pt < TempPoly2.Size; pt++)
			{
				int TestDeg = 0;
				for (int pd = 0; pd < 288; pd++)
					if ((TempPoly2.poly[pt].pterm[(pd >> 5)] >> (31 - (pd & 0x1f))) & 1)
					{
						if ((pd != 287) && ((TempPoly2.poly[pt].pterm[((pd + 1) >> 5)] >> (31 - ((pd + 1) & 0x1f))) & 1))
							TestDeg += SState[pd++];
						else
							TestDeg += State[pd];
					}
				if (TestDeg >= dim)
					TempPoly2.poly[pf++] = TempPoly2.poly[pt];
			}
			TempPoly2.Size = pf;
			cout << "DegEst1(Pass Numeric Mapping): OutPutANF Terms: " << TempPoly2.Size << ", SaveTerms: " << InPutANF.Size << ", InPutRound = " << InPutRound << endl;
			// ---------------------------------------------------------------


			int countflag = 0;
			for (int i = 0; i < TempPoly2.Size; i++)
				if (TempPoly2.poly[i].getdegk() < 0)
					countflag += 1;
			cout << "The true number of monomials whose algebraic degrees need to be evaluted: " << countflag << endl;
			if ( ((countflag < BackExpandNumThr) && (InPutRound > 200)) || (countflag < 500) && (InPutRound <= 200))
			{
				OutPutANF.PolyCopy(TempPoly2);
				TESTFLAG = true;
				continue;
			}
			TESTFLAG = false;


			string files1 = "RoundLog//TempFile(" + to_string(InPutRound) + ").txt";
			char filec[100];
			files1.copy(filec, files1.length(), 0);
			*(filec + files1.length()) = '\0';

			TempPoly2.WriteValueToFile(InPutRound, filec, "TermNum: " + to_string(TempPoly2.Size), 1);
			InPutANF.WriteValueToFile(InPutRound, filec, "SaveTerm: " + to_string(InPutANF.Size));
			RetSuperpoly.WritePolyValueToFile(filec, "Recovered Superpoly", 1);
			RetSuperpoly.WritePolyToFile(filec, "Recovered Superpoly");


			int CCount = 0;
#pragma omp parallel for schedule(dynamic,100) num_threads(ThreadNumForLoop)
			for (int64_t pd = 0; pd < TempPoly2.Size; pd++)
			{
				if (TempPoly2.poly[pd].getdegk() < 0)
				{
					MapPoly TempPoly;
					clock_t ModelStartTime = clock();
					int TermDeg = TriviumEvalBDPT(TempPoly, TempPoly2.poly[pd], CubeFlag, InPutRound, false, (InPutRound >= 400) * 30);
					TempPoly2.poly[pd].setdegk(TermDeg);
					clock_t ModelEndTime = clock();

#pragma omp atomic
					CCount += 1;

					if (TermDeg > 0)
					{
						string message = "Round: " + to_string(InPutRound) + ", ThreadId: " + to_string(omp_get_thread_num()) + ", Have Tested Terms(" + to_string(TermDeg) + "): " +
							to_string(CCount) + "/" + to_string(TempPoly2.Size) + ", RunningTime: " + to_string(double(ModelEndTime - ModelStartTime) / CLOCKS_PER_SEC) + "s \n";
						cout << message;
					}
				}
				else
					continue;
			}
			
			for (int64_t pd = 0; pd < TempPoly2.Size; pd++)
			{
				int deg = TempPoly2.poly[pd].getdegk();
				if (deg > 79)
					InPutANF.poly[InPutANF.Size++] = TempPoly2.poly[pd];
				else if (deg > 0)
					TempPoly3.poly[TempPoly3.Size++] = TempPoly2.poly[pd];
			}
			quickSortDeg(TempPoly3.poly, int64_t(0), int64_t(TempPoly3.Size) - 1);
			if (TempPoly3.Size > 0)
				MaxDeg = TempPoly3.poly[TempPoly3.Size - 1].getdegk();

			
			int *TermNum = new int[MaxDeg];
			for (int i = 0; i < MaxDeg; i++)
				TermNum[i] = 0;
			for (int i = 0; i < TempPoly3.Size; i++)
				TermNum[TempPoly3.poly[i].getdegk() - 1] += 1;
			cout << "(deg, MonomialNumber): ";
			for (int i = 0; i < MaxDeg; i++)
				cout << '(' << i << ", " << TermNum[i] << ") ";
			cout << endl;


			for (int i = 0; i < LEN; i++)
				if (InPutRound >= TRound[i])
					TimeLimit = Ttime[i] * TimeLimitRatio;
			int AllTermNum = TempPoly3.Size, pd = 0, dd = 0;

			for (; dd < MaxDeg; dd++)
				if (TermNum[dd] > 0)
					break;

			while (dd < MaxDeg)
			{
				int CCOUNT = TermNum[dd];
				int CCount = 0;

				for (int i = 0; i < ParallelNum; i++)
				{
					ThreadTempPoly[i].Size = 0;
					ThreadTempSP[i].clear();
				}

#pragma omp parallel for schedule(dynamic, 5) num_threads(ParallelNum)
				for (int ppd = pd; ppd < pd + CCOUNT; ppd++)
				{
					clock_t ModelStartTime = clock();
					MapPoly TempPoly; TempPoly.clear();
					int loc = omp_get_thread_num();
					bool TimeLimteFlag = TriviumEvalBDPT(TempPoly, TempPoly3.poly[ppd], CubeFlag, InPutRound, true, TimeLimit);

					if (TimeLimteFlag)
						ThreadTempPoly[loc].poly[ThreadTempPoly[loc].Size++] = TempPoly3.poly[ppd];
					else
						ThreadTempSP[loc].Merge(TempPoly);
					clock_t ModelEndTime = clock();

#pragma omp atomic
					CCount += 1;

					string message = "Round: " + to_string(InPutRound) + ", ThreadId: " + to_string(loc) + ", Have Tested Terms(" + to_string(TempPoly3.poly[ppd].getdegk())
						+ "): (" + to_string(ppd + 1) + ", " + to_string(CCount) + "/" + to_string(TermNum[dd]) + "/" + to_string(AllTermNum) + "), The size of the Superpoly: "
						+ to_string(ThreadTempSP[loc].size()) + ", Accepted Terms: " + to_string(ThreadTempPoly[loc].Size) + ", RunningTime: "
						+ to_string(double(ModelEndTime - ModelStartTime) / CLOCKS_PER_SEC) + "s \n";
					cout << message;
				}

				AllTermNum -= TermNum[dd];
				pd += CCOUNT; dd += 1;
				for (int i = 0; i < ParallelNum; i++)
					for (int j = 0; j < ThreadTempPoly[i].Size; j++)
						InPutANF.poly[InPutANF.Size++] = ThreadTempPoly[i].poly[j];
				for (int i = 0; i < ParallelNum; i++)
					RetSuperpoly.Merge(ThreadTempSP[i]);


				TempPoly3.WriteValueToFile(InPutRound, filec, to_string(pd) + "/" + to_string(TempPoly3.Size), 1);
				InPutANF.WriteValueToFile(InPutRound, filec, "SaveTerm: " + to_string(InPutANF.Size));
				RetSuperpoly.WritePolyValueToFile(filec, "Recovered Superpoly", 1);
				RetSuperpoly.WritePolyToFile(filec, "Recovered Superpoly");
				cout << "Round: " << InPutRound << ", The size of the Superpoly: " << RetSuperpoly.size() << ", Accepted Terms: " << InPutANF.Size << endl;

				
				int TimeLimitTermNum = 0;
				for (int i = 0; i < ParallelNum; i++)
					TimeLimitTermNum += ThreadTempPoly[i].Size;
				if (TimeLimitTermNum > ParallelNum*5)
					break;
			}
			for (; pd < TempPoly3.Size; pd++)
				OutPutANF.poly[OutPutANF.Size++] = TempPoly3.poly[pd];
			delete[] TermNum;

			if ( ((OutPutANF.Size < 3000) && (InPutRound >= 200)) || ((OutPutANF.Size < 500) && (InPutRound < 200)) )
			{
				InPutANF.Merge(OutPutANF);
				OutPutANF.Size = 0;
			}
			
			OutPutANF.WriteValueToFile(InPutRound, filec, "OutPutANF: " + to_string(OutPutANF.Size), 1);
			InPutANF.WriteValueToFile(InPutRound, filec, "InPutANF: " + to_string(InPutANF.Size));
			RetSuperpoly.WritePolyValueToFile(filec, "Recovered Superpoly", 1);
			RetSuperpoly.WritePolyToFile(filec, "Recovered Superpoly");
			cout << "DegEst2(Pass MILP Degest): OutPutANF Terms: " << OutPutANF.Size << ", SaveTerms: " << InPutANF.Size << ", InPutRound = " << InPutRound << endl;
		}

		if ((InPutRound == 0) || ((InPutANF.Size + OutPutANF.Size) == 0))
		{
			char OutFile1[100], AppOutFile1[100];
			string files2 = files.substr(0, files.length() - 4) + '(' + to_string(InPutRound) + ").txt";
			files2.copy(OutFile1, files2.length(), 0);
			*(OutFile1 + files2.length()) = '\0';


			appfiles.copy(AppOutFile1, appfiles.length(), 0);
			*(AppOutFile1 + appfiles.length()) = '\0';


			InPutANF.WriteValueToFile(InPutRound, OutFile1, "Terms: " + to_string(InPutANF.Size), true);
			RetSuperpoly.WritePolyValueToFile(OutFile1, "Recovered Superpoly", 1);
			RetSuperpoly.WritePolyToFile(OutFile1, "Recovered Superpoly");


			InPutANF.WriteValueToFile(InPutRound, AppOutFile1, "Terms: " + to_string(InPutANF.Size));
			RetSuperpoly.WritePolyValueToFile(AppOutFile1, "Recovered Superpoly", 1);
			RetSuperpoly.WritePolyToFile(AppOutFile1, "Recovered Superpoly");
			break;
		}

		quickSortDeg(OutPutANF.poly, int64_t(0), int64_t(OutPutANF.Size) - 1);
		for (int i = 0; i < OutPutANF.Size; i++)
			OutPutANF.poly[i].show("degk: " + to_string(OutPutANF.poly[i].getdegk()));
		MaxDeg = OutPutANF.poly[OutPutANF.Size - 1].getdegk();

		clock_t TotalStartTime = clock();

		int *TermNum = new int[MaxDeg];
		for (int i = 0; i < MaxDeg; i++)
			TermNum[i] = 0;
		for (int i = 0; i < OutPutANF.Size; i++)
			TermNum[OutPutANF.poly[i].getdegk() - 1] += 1;
		cout << "(deg, MonomialNumber): ";
		for (int i = 0; i < MaxDeg; i++)
			cout << '(' << i << ", " << TermNum[i] << ") ";
		cout << endl;

		string files1 = "RoundLog//TempFile(" + to_string(InPutRound) + ").txt";
		char filec[100];
		files1.copy(filec, files1.length(), 0);
		*(filec + files1.length()) = '\0';

		for (int i = 0; i < LEN; i++)
			if (InPutRound >= TRound[i])
				TimeLimit = Ttime[i] * TimeLimitRatio;
		TimeLimitQ *= TimeLimitRatio;
		int AllTermNum = OutPutANF.Size, pd = 0, dd = 0;

		for (; dd < MaxDeg; dd++)
			if (TermNum[dd] > 0)
				break;
		
		while (dd < MaxDeg)
		{
			int CCOUNT = TermNum[dd];
			int CCount = 0;

			for (int i = 0; i < ParallelNum; i++)
			{
				ThreadTempPoly[i].Size = 0;
				ThreadTempSP[i].clear();
			}

#pragma omp parallel for schedule(dynamic, 5) num_threads(ParallelNum)
			for (int ppd = pd; ppd < pd + CCOUNT; ppd++)
			{
				clock_t ModelStartTime = clock();
				MapPoly TempPoly; TempPoly.clear();
				bool TimeLimteFlag = 0;

				TimeLimteFlag = TriviumEvalBDPT(TempPoly, OutPutANF.poly[ppd], CubeFlag, InPutRound, true, TimeLimit);

				int loc = omp_get_thread_num();

				if (TimeLimteFlag)
					ThreadTempPoly[loc].poly[ThreadTempPoly[loc].Size++] = OutPutANF.poly[ppd];
				else
					ThreadTempSP[loc].Merge(TempPoly);
				clock_t ModelEndTime = clock();

#pragma omp atomic
				CCount += 1;

				string message = "Round: " + to_string(InPutRound) + ", ThreadId: " + to_string(loc) + ", Have Tested Terms(" + to_string(OutPutANF.poly[ppd].getdegk())
					+ "): (" + to_string(ppd + 1) + ", " + to_string(CCount) + "/" + to_string(TermNum[dd]) + "/" + to_string(AllTermNum) + "), The size of the Superpoly: "
					+ to_string(ThreadTempSP[loc].size()) + ", Accepted Terms: " + to_string(ThreadTempPoly[loc].Size) + ", RunningTime: "
					+ to_string(double(ModelEndTime - ModelStartTime) / CLOCKS_PER_SEC) + "s \n";
				cout << message;
			}

			AllTermNum -= TermNum[dd];
			for (int i = 0; i < ParallelNum; i++)
				for (int j = 0; j < ThreadTempPoly[i].Size; j++)
					InPutANF.poly[InPutANF.Size++] = ThreadTempPoly[i].poly[j];
			for (int i = 0; i < ParallelNum; i++)
				RetSuperpoly.Merge(ThreadTempSP[i]);

			pd += CCOUNT; dd += 1;
			OutPutANF.WriteValueToFile(InPutRound, filec, to_string(pd) + "/" + to_string(OutPutANF.Size), 1);
			InPutANF.WriteValueToFile(InPutRound, filec, "SaveTerm: " + to_string(InPutANF.Size));
			RetSuperpoly.WritePolyValueToFile(filec, "Recovered Superpoly", 1);
			RetSuperpoly.WritePolyToFile(filec, "Recovered Superpoly");
			cout << "Round: " << InPutRound << ", The size of the Superpoly: " << RetSuperpoly.size() << ", Accepted Terms: " << InPutANF.Size << endl;

			int TimeLimitTermNum = 0;
			for (int i = 0; i < ParallelNum; i++)
				TimeLimitTermNum += ThreadTempPoly[i].Size;
			if (TimeLimitTermNum > TermNum[dd - 1] * ThresholdCoeff)
				break;
		}

		while (dd < MaxDeg)
		{
			int CCOUNT = TermNum[dd];
			int CCount = 0;

			for (int i = 0; i < ParallelNum; i++)
			{
				ThreadTempPoly[i].Size = 0;
				ThreadTempSP[i].clear();
			}

#pragma omp parallel for schedule(dynamic, 5) num_threads(ParallelNum)
			for (int ppd = pd; ppd < pd + CCOUNT; ppd++)
			{
				clock_t ModelStartTime = clock();
				MapPoly TempPoly; TempPoly.clear();
				bool TimeLimteFlag = 0;

				TimeLimteFlag = TriviumEvalBDPT(TempPoly, OutPutANF.poly[ppd], CubeFlag, InPutRound, true, TimeLimitQ);

				int loc = omp_get_thread_num();

				if (TimeLimteFlag)
					ThreadTempPoly[loc].poly[ThreadTempPoly[loc].Size++] = OutPutANF.poly[ppd];
				else
					ThreadTempSP[loc].Merge(TempPoly);
				clock_t ModelEndTime = clock();

#pragma omp atomic
				CCount += 1;

				string message = "Round: " + to_string(InPutRound) + ", ThreadId: " + to_string(loc) + ", Have Tested Terms(" + to_string(OutPutANF.poly[ppd].getdegk())
					+ "): (" + to_string(ppd + 1) + ", " + to_string(CCount) + "/" + to_string(TermNum[dd]) + "/" + to_string(AllTermNum)
					+ "), The size of the Superpoly: " + to_string(ThreadTempSP[loc].size()) + ", Accepted Terms: " + to_string(ThreadTempPoly[loc].Size)
					+ ", RunningTime: " + to_string(double(ModelEndTime - ModelStartTime) / CLOCKS_PER_SEC) + "s \n";
				cout << message;
			}

			AllTermNum -= TermNum[dd];
			for (int i = 0; i < ParallelNum; i++)
				for (int j = 0; j < ThreadTempPoly[i].Size; j++)
					InPutANF.poly[InPutANF.Size++] = ThreadTempPoly[i].poly[j];
			for (int i = 0; i < ParallelNum; i++)
				RetSuperpoly.Merge(ThreadTempSP[i]);


			pd += CCOUNT; dd += 1;
			OutPutANF.WriteValueToFile(InPutRound, filec, to_string(pd) + "/" + to_string(OutPutANF.Size), 1);
			InPutANF.WriteValueToFile(InPutRound, filec, "SaveTerm: " + to_string(InPutANF.Size));
			RetSuperpoly.WritePolyValueToFile(filec, "Recovered Superpoly", 1);
			RetSuperpoly.WritePolyToFile(filec, "Recovered Superpoly");
			cout << "Round: " << InPutRound << ", The size of the Superpoly: " << RetSuperpoly.size() << ", Accepted Terms: " << InPutANF.Size << endl;

			
			int TimeLimitTermNum = 0;
			for (int i = 0; i < ParallelNum; i++)
				TimeLimitTermNum += ThreadTempPoly[i].Size;
			if (TimeLimitTermNum > TermNum[dd - 1] * 0.7)
				break;
		}
		for (; pd < OutPutANF.Size; pd++)
			InPutANF.poly[InPutANF.Size++] = OutPutANF.poly[pd];
		cout << "Round: " << InPutRound << ", The size of the Superpoly: " << RetSuperpoly.size() << ", Accepted Terms: " << InPutANF.Size << endl;

		clock_t TotalEndTime = clock();
		delete[] TermNum;
		// -----------------------------------------------------------------------------------------------------

		InPutANF.RemoveDup();
		cout << "Save Terms(Pass BDPT_MILP): " << InPutANF.Size << ", Rumtime: " << double(TotalEndTime - TotalStartTime) / CLOCKS_PER_SEC << endl;


		char OutFile[100], AppOutFile[100];
        string files2 = files.substr(0, files.length() - 4) + '(' + to_string(InPutRound) + ").txt";
		files2.copy(OutFile, files2.length(), 0);
		*(OutFile + files2.length()) = '\0';

		appfiles.copy(AppOutFile, appfiles.length(), 0);
		*(AppOutFile + appfiles.length()) = '\0';


		InPutANF.WriteValueToFile(InPutRound, OutFile, "Terms: " + to_string(InPutANF.Size), true);
		RetSuperpoly.WritePolyValueToFile(OutFile, "Recovered Superpoly", 1);
		RetSuperpoly.WritePolyToFile(OutFile, "Recovered Superpoly");

		InPutANF.WriteValueToFile(InPutRound, AppOutFile, "Terms: " + to_string(InPutANF.Size));
		RetSuperpoly.WritePolyValueToFile(AppOutFile, "Recovered Superpoly", 1);
		RetSuperpoly.WritePolyToFile(AppOutFile, "Recovered Superpoly");

		SaveTermNum = InPutANF.Size;
		OutPutANF.Size = 0;
		if (SaveTermNum == 0)
			return true;
	}

	delete[] ThreadTempPoly;
	delete[] ThreadTempSP;

	if (SaveTermNum != 0)
		return false;
	else
		return true;
}



bool FilterPartTermBDPT0(MapPoly& RetSuperpoly, int& InPutRound, Poly& ANF, Ivterm &CubeFlag, string files, string appfiles)
{
	const int dim = CubeFlag.count();

	cout << "ReadANF.Size = " << ANF.Size << ", ";
	ANF.RemoveDup();
	cout << "RemoveDup(ReadANF).Size = " << ANF.Size << endl;

	Poly InPutANF(ANF.Size);
	
	int AllThread = std::thread::hardware_concurrency();
	cout << "The total number of threads available to this computer: " << AllThread << endl;
	omp_set_num_threads(AllThread);
	omp_set_nested(1);
	system("mkdir RoundLog");

	int TimeLimit = 60, TimeLimitQ = 60;
	int ThreadNumForLoop = AllThread >> 2;
	int CCount = 0;

#pragma omp parallel for schedule(dynamic,100) num_threads(ThreadNumForLoop)
	for (int pd = 0; pd < ANF.Size; pd++)
	{
		if (ANF.poly[pd].getdegk() < 0)
		{
			MapPoly TempPoly;
			clock_t ModelStartTime = clock();
			int TermDeg = TriviumEvalBDPT(TempPoly, ANF.poly[pd], CubeFlag, InPutRound, false, (InPutRound>=450)*30);
			ANF.poly[pd].setdegk(TermDeg);
			clock_t ModelEndTime = clock();

#pragma omp atomic
			CCount += 1;

			if (TermDeg > 0)
			{
				string message = "Round: " + to_string(InPutRound) + ", ThreadId: " + to_string(omp_get_thread_num()) + ", Have Tested Terms(" + to_string(TermDeg) + "): " +
					to_string(CCount) + "/" + to_string(ANF.Size) + ", RunningTime: " + to_string(double(ModelEndTime - ModelStartTime) / CLOCKS_PER_SEC) + "s \n";
				cout << message;
			}
		}
		else
			continue;
	}

	for (int pd = 0; pd < ANF.Size; pd++)
		if (ANF.poly[pd].getdegk() > 0)
			InPutANF.poly[InPutANF.Size++] = ANF.poly[pd];
	ANF.Size = 0;

	cout << "Start Sort!" << endl;
	quickSortDeg(InPutANF.poly, int64_t(0), int64_t(InPutANF.Size) - 1);


	clock_t TotalStartTime = clock();
	int MaxDeg = InPutANF.poly[InPutANF.Size - 1].getdegk();
	int *TermNum = new int[MaxDeg];
	for (int i = 0; i < MaxDeg; i++)
		TermNum[i] = 0;
	for (int i = 0; i < InPutANF.Size; i++)
		TermNum[InPutANF.poly[i].getdegk() - 1] += 1;
	cout << "(deg, MonomialNumber): ";
	for (int i = 0; i < MaxDeg; i++)
		cout << '(' << i << ", " << TermNum[i] << ") ";
	cout << endl;


	const int LEN = 7;
	int TRound[LEN] = { 100, 150, 200, 250, 300, 350, 400 };
	int Ttime[LEN] = { 270, 240, 210, 180, 150, 120, 90 };
	int TimeLimitRatio = 1;
	int SkipDeg = min(0, MaxDeg);
	int ParallelNum = 8;
	double ThresholdCoeff1 = 0.50;

	int pd = 0, dd = 0;
	string files1 = "RoundLog//TempFile(" + to_string(InPutRound) + ").txt";
	char filec[100];
	files1.copy(filec, files1.length(), 0);
	*(filec + files1.length()) = '\0';
	Poly *ThreadTempPoly = new Poly[AllThread];
	MapPoly *ThreadTempSP = new MapPoly[AllThread];
	int AllTermNum = InPutANF.Size;
	for (int i = 0; i < LEN; i++)
		if (InPutRound >= TRound[i])
			TimeLimit = Ttime[i];
	TimeLimit *= TimeLimitRatio;
	TimeLimitQ *= TimeLimitRatio;


	while (dd < SkipDeg)
	{
		int CCOUNT = TermNum[dd];
		for (int ppd = pd; ppd < pd + CCOUNT; ppd++)
			ANF.poly[ANF.Size++] = InPutANF.poly[ppd];
		AllTermNum -= CCOUNT;
		pd += CCOUNT;
		dd += 1;
		cout << "Round: " << InPutRound << ", Termdeg: " << dd << ", The size of the Superpoly: " << RetSuperpoly.size() << ", Accepted Terms: " << ANF.Size << endl;
	}
	InPutANF.WriteValueToFile(InPutRound, filec, to_string(pd) + "/" + to_string(InPutANF.Size), 1);
	ANF.WriteValueToFile(InPutRound, filec, "SaveTerm: " + to_string(ANF.Size));
	RetSuperpoly.WritePolyValueToFile(filec, "Recovered Superpoly", 1);
	RetSuperpoly.WritePolyToFile(filec, "Recovered Superpoly");


	for (; dd < MaxDeg; dd++)
		if (TermNum[dd] > 0)
			break;
	while (dd < MaxDeg)
	{
		int CCOUNT = TermNum[dd];
		int CCount = 0;

		for (int i = 0; i < ParallelNum; i++)
		{
			ThreadTempPoly[i].Size = 0;
			ThreadTempSP[i].clear();
		}

#pragma omp parallel for schedule(dynamic, 20) num_threads(ParallelNum)
		for (int ppd = pd; ppd < pd + CCOUNT; ppd++)
		{
			MapPoly TempPoly; TempPoly.clear();
			int loc = omp_get_thread_num();

			clock_t ModelStartTime = clock();
			bool TimeLimteFlag = TriviumEvalBDPT(TempPoly, InPutANF.poly[ppd], CubeFlag, InPutRound, true, TimeLimit);  // 120  // 还得再改回来
			if (TimeLimteFlag)
				ThreadTempPoly[loc].poly[ThreadTempPoly[loc].Size++] = InPutANF.poly[ppd];			
			else
				ThreadTempSP[loc].Merge(TempPoly);
			clock_t ModelEndTime = clock();

#pragma omp atomic
			CCount += 1;

			string message = "Round: " + to_string(InPutRound) + ", ThreadId: " + to_string(loc) + ", Have Tested Terms(" + to_string(InPutANF.poly[ppd].getdegk()) + "): ("
				+ to_string(ppd + 1) + ", " + to_string(CCount) + "/" + to_string(TermNum[dd]) + "/" + to_string(AllTermNum) + "), The size of the Superpoly: " + to_string(ThreadTempSP[loc].size())
				+ ", Accepted Terms: " + to_string(ThreadTempPoly[loc].Size) + ", RunningTime: " + to_string(double(ModelEndTime - ModelStartTime) / CLOCKS_PER_SEC) + "s\n";
			cout << message;
		}

		AllTermNum -= TermNum[dd];
		for (int i = 0; i < ParallelNum; i++)
			for (int j = 0; j < ThreadTempPoly[i].Size; j++)
				ANF.poly[ANF.Size++] = ThreadTempPoly[i].poly[j];
		for (int i = 0; i < ParallelNum; i++)
			RetSuperpoly.Merge(ThreadTempSP[i]);

		cout << "Round: " << InPutRound << ", deg: " << dd+1 << ", The size of the Superpoly: " << RetSuperpoly.size() << ", Accepted Terms: " << ANF.Size << endl;

		pd += CCOUNT; 
		InPutANF.WriteValueToFile(InPutRound, filec, to_string(pd) + "/" + to_string(InPutANF.Size), 1);
		ANF.WriteValueToFile(InPutRound, filec, "SaveTerm: " + to_string(ANF.Size));
		RetSuperpoly.WritePolyValueToFile(filec, "Recovered Superpoly", 1);
		RetSuperpoly.WritePolyToFile(filec, "Recovered Superpoly");

		dd += 1;
		int TimeLimitTermNum = 0;
		for (int i = 0; i < ParallelNum; i++)
			TimeLimitTermNum += ThreadTempPoly[i].Size;
		if (TimeLimitTermNum > TermNum[dd - 1] * ThresholdCoeff1)
			break;
	}

	for (; dd < MaxDeg; dd++)
	{
		int CCOUNT = TermNum[dd];
		int CCount = 0;

		for (int i = 0; i < ParallelNum; i++)
		{
			ThreadTempPoly[i].Size = 0;
			ThreadTempSP[i].clear();
		}

#pragma omp parallel for schedule(dynamic, 5) num_threads(ParallelNum)
		for (int ppd = pd; ppd < pd + CCOUNT; ppd++)
		{
			MapPoly TempPoly; TempPoly.clear();
			int loc = omp_get_thread_num();

			clock_t ModelStartTime = clock();
			bool TimeLimteFlag = TriviumEvalBDPT(TempPoly, InPutANF.poly[ppd], CubeFlag, InPutRound, true, TimeLimitQ);
			if (TimeLimteFlag)
				ThreadTempPoly[loc].poly[ThreadTempPoly[loc].Size++] = InPutANF.poly[ppd];
			else
				ThreadTempSP[loc].Merge(TempPoly);
			clock_t ModelEndTime = clock();


#pragma omp atomic
			CCount += 1;

			string message = "Round: " + to_string(InPutRound) + ", ThreadId: " + to_string(loc) + ", Have Tested Terms(" + to_string(InPutANF.poly[ppd].getdegk()) + "): ("
				+ to_string(ppd + 1) + ", " + to_string(CCount) + "/" + to_string(TermNum[dd]) + "/" + to_string(AllTermNum) + "), The size of the Superpoly: " + to_string(ThreadTempSP[loc].size())
				+ ", Accepted Terms: " + to_string(ThreadTempPoly[loc].Size) + ", RunningTime: " + to_string(double(ModelEndTime - ModelStartTime) / CLOCKS_PER_SEC) + "s\n";
			
			cout << message;
		}

		for (int i = 0; i < ParallelNum; i++)
			for (int j = 0; j < ThreadTempPoly[i].Size; j++)
				ANF.poly[ANF.Size++] = ThreadTempPoly[i].poly[j];
		for (int i = 0; i < ParallelNum; i++)
			RetSuperpoly.Merge(ThreadTempSP[i]);
		AllTermNum -= TermNum[dd];

		cout << "Round: " << InPutRound << ", The size of the Superpoly: " << RetSuperpoly.size() << ", Accepted Terms: " << ANF.Size << endl;

		pd += CCOUNT;
		InPutANF.WriteValueToFile(InPutRound, filec, to_string(pd) + "/" + to_string(InPutANF.Size), 1);
		ANF.WriteValueToFile(InPutRound, filec, "SaveTerm: " + to_string(ANF.Size));
		RetSuperpoly.WritePolyValueToFile(filec, "Recovered Superpoly", 1);
		RetSuperpoly.WritePolyToFile(filec, "Recovered Superpoly");
	}
	InPutANF.Size = 0;

	delete[] ThreadTempPoly;
	delete[] ThreadTempSP;
	delete[] TermNum;

	cout << "Round: " << InPutRound << ", The size of the Superpoly: " << RetSuperpoly.size() << ", Accepted Terms: " << ANF.Size << endl;

	clock_t TotalEndTime = clock();
	// -----------------------------------------------------------------------------------------------------

	cout << "Save Terms(Pass BDPT_MILP): " << InPutANF.Size << ", Rumtime: " << double(TotalEndTime - TotalStartTime) / CLOCKS_PER_SEC << endl;


	char OutFile[100], AppOutFile[100];
	files = files.substr(0, files.length() - 4) + '(' + to_string(InPutRound) + ").txt";
	files.copy(OutFile, files.length(), 0);
	*(OutFile + files.length()) = '\0';

	appfiles.copy(AppOutFile, appfiles.length(), 0);
	*(AppOutFile + appfiles.length()) = '\0';

	ANF.WriteValueToFile(InPutRound, OutFile, "Terms: " + to_string(InPutANF.Size), true);
	RetSuperpoly.WritePolyValueToFile(OutFile, "Recovered Superpoly", 1);
	RetSuperpoly.WritePolyToFile(OutFile, "Recovered Superpoly");

	ANF.WriteValueToFile(InPutRound, AppOutFile, "Terms: " + to_string(InPutANF.Size));
	RetSuperpoly.WritePolyValueToFile(AppOutFile, "Recovered Superpoly", 1);
	RetSuperpoly.WritePolyToFile(AppOutFile, "Recovered Superpoly");
}



int TriviumEvalBDPT(MapPoly& ResPoly, const Term &Pterm, Ivterm &Cube, int Round, bool SolutFlag, int TimeLimit, int ThreadNum)
{
	Round -= 64;
	GRBEnv Env = GRBEnv();
	Env.set(GRB_IntParam_LogToConsole, 0);
	if (ThreadNum != 0)
		Env.set(GRB_IntParam_Threads, ThreadNum);
	GRBModel Model = GRBModel(Env);
	if (TimeLimit != 0)
		Model.set(GRB_DoubleParam_TimeLimit, TimeLimit);

	if (SolutFlag)
	{
		Model.set(GRB_IntParam_PoolSolutions, 2000000000);  // Limit how many solutions to collect
		Model.set(GRB_IntParam_PoolSearchMode, 2);  // do a systematic search for the k- best solutions
		Model.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);  // Limit the search space by setting a gap for the worst possible solution that will be acce
		Model.set(GRB_IntParam_MIPFocus, 3);
	}

	int Loc1[5] = { 65, 170, 90, 91, 92 };
	int Loc2[5] = { 161, 263, 174, 175, 176 };
	int Loc3[5] = { 242, 68, 285, 286, 287 };

	GRBVar* kvar = Model.addVars(MAPPOLYLEN, 'B');
	GRBVar* vvar = Model.addVars(80, 'B');
	GRBVar* svar = Model.addVars(288, 'B');
	for (int i = 69; i < 80; i++)
		Model.addConstr(kvar[i] == 0);
	GRBLinExpr key_sum;
	for (int i = 0; i < MAPPOLYLEN; i++)
		key_sum += kvar[i];
	Model.setObjective(key_sum, GRB_MAXIMIZE);


	for (int i = 0; i < 80; i++)
		if (Cube.test(i))
			Model.addConstr(vvar[i] == 1, "Iv_constr_" + to_string(i));
		else
			Model.addConstr(vvar[i] == 0, "Iv_constr_" + to_string(i));

	Initilization(Model, svar, kvar, vvar);


	for (int r = 0; r < Round; r++)
	{
		TriviumCoreBDPT(Model, svar, Loc1, r);
		TriviumCoreBDPT(Model, svar, Loc2, r);
		TriviumCoreBDPT(Model, svar, Loc3, r);
		GRBVar temp = svar[287];
		for (int i = 287; i > 0; i--)
			svar[i] = svar[i - 1];
		svar[0] = temp;
	}

	// GRBLinExpr TermConstr;
	for (int pd = 0; pd < 288; pd++)
	{
		if ((Pterm.pterm[pd >> 5] >> (31 - (pd & 0x1f))) & 1)
			Model.addConstr(svar[pd] == 1, "FinalRoundState" + to_string(pd) + " to 1");
		else
			Model.addConstr(svar[pd] == 0, "FinalRoundState" + to_string(pd) + " to 0");
	}

	Model.optimize();

	if (Model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
	{
		if (SolutFlag == false)
			return int(Model.get(GRB_DoubleAttr_ObjVal)) + 1;

		int nSolutions = Model.get(GRB_IntAttr_SolCount);
		cout << "- Number of solutions found : " << nSolutions << ", RunningTime: " << Model.get(GRB_DoubleAttr_Runtime) << 's' << endl;

		for (int res = 0; res < nSolutions; res++)
		{
			KeyTerm p_term;
			Model.set(GRB_IntParam_SolutionNumber, res);

			for (int loc = 0; loc < MAPPOLYLEN; loc++)
			{
				if (kvar[loc].get(GRB_DoubleAttr_Xn) >= 0.5)
					p_term.set(loc);
			}
			ResPoly.InSert(p_term);
		}
		return false;
	}
	else if (Model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
	{
		if (SolutFlag)
			return true;
		else
			return MAPPOLYLEN;
	}
	else
		return false;
}


void Initilization(GRBModel &Model, GRBVar* sVar, GRBVar *kVar, GRBVar *vVar)
{
	// anf key(0-131), v(132-212)
	vector<vector<int>> anfs[288];

	for (int i = 0; i < 64; i++)
		anfs[i] = { { 5 + i } };
	for (int i = 64; i < 93; i++)
		anfs[i] = { { i - 64 } };

	for (int i = 93; i < 145; i++)
		anfs[i] = { { i - 93 + MAPPOLYLEN + 14 },{ i - 93 + 80 } };
	for (int i = 145; i < 157; i++)
		anfs[i] = { { i - 145 + MAPPOLYLEN + 66 },{ i - 145 + 54 } };
	for (int i = 157; i < 177; i++)
		anfs[i] = { { i - 157 + MAPPOLYLEN } };

	for (int i = 177; i < 237; i++)
		anfs[i] = { { i - 177 + MAPPOLYLEN + 18, i - 177 + MAPPOLYLEN + 19 },{ i - 177 + MAPPOLYLEN + 5 },{ i - 177 + MAPPOLYLEN + 20 } };
	anfs[237] = { { 78 + MAPPOLYLEN, 79 + MAPPOLYLEN },{ 65 + MAPPOLYLEN } };
	anfs[238] = { { 66 + MAPPOLYLEN } };
	anfs[239] = { { 67 + MAPPOLYLEN } };
	anfs[240] = { { 68 + MAPPOLYLEN } };


	vector<int> varcount1, varcount2;
	for (int i = 0; i < MAPPOLYLEN + 80; i++)
	{
		varcount1.push_back(0);
		varcount2.push_back(0);
	}

	for (int i = 0; i < 288; i++)
		for (int j = 0; j < anfs[i].size(); j++)
			for (vector<int>::iterator pt = anfs[i][j].begin(); pt != anfs[i][j].end(); pt++)
				varcount1[*pt] += 1;

	GRBVar **newvar = new GRBVar*[MAPPOLYLEN + 80];
	for (int i = 0; i < MAPPOLYLEN; i++)
		if (varcount1[i] <= 1)
			newvar[i] = &kVar[i];
		else
		{
			newvar[i] = Model.addVars(varcount1[i], 'B');
			ModelCopy(Model, kVar[i], newvar[i], varcount1[i]);
		}
	for (int i = MAPPOLYLEN; i < MAPPOLYLEN + 80; i++)
		if (varcount1[i] <= 1)
			newvar[i] = &vVar[i - MAPPOLYLEN];
		else
		{
			newvar[i] = Model.addVars(varcount1[i], 'B');
			ModelCopy(Model, vVar[i - MAPPOLYLEN], newvar[i], varcount1[i]);
		}

	for (int i = 0; i < 288; i++)
	{
		if (anfs[i].size() == 0)
		{
			Model.addConstr(sVar[i] == 0);
			continue;
		}

		GRBLinExpr sstemp;
		for (int j = 0; j < anfs[i].size(); j++)
		{
			if (anfs[i][j].size() == 2)
				Model.addConstr(newvar[anfs[i][j][0]][varcount2[anfs[i][j][0]]] == newvar[anfs[i][j][1]][varcount2[anfs[i][j][1]]++]);
			sstemp += newvar[anfs[i][j][0]][varcount2[anfs[i][j][0]]++];
		}
		Model.addConstr(sVar[i] == sstemp);
	}


#if 0  // print state
	for (int i = 0; i < MAPPOLYLEN; i++)
		cout << 'k' << i << '(' << varcount1[i] << ") ";
	cout << endl;
	for (int j = MAPPOLYLEN; j < MAPPOLYLEN + 80; j++)
		cout << 'v' << j - MAPPOLYLEN << '(' << varcount1[j] << ") ";
	cout << endl;

	for (int i = 0; i < 288; i++)
	{
		cout << 's' << i << ": ";

		if (anfs[i].size() == 0)
		{
			cout << 0 << endl;
			continue;
		}

		for (int j = 0; j < anfs[i].size(); j++)
		{
			if (j > 0)
				cout << " + ";
			for (vector<int>::iterator pt = anfs[i][j].begin(); pt != anfs[i][j].end(); pt++)
				if (*pt < MAPPOLYLEN)
					cout << 'k' << *pt;
				else
					cout << 'v' << *pt - MAPPOLYLEN;
		}
		cout << endl;
	}
#endif
}


void ModelCopy(GRBModel &Model, GRBVar &var, GRBVar* newvar, int len)
{
	GRBLinExpr varsum;
	for (int i = 0; i < len; i++)
		varsum += newvar[i];
	Model.addConstr(varsum >= var);
	for (int i = 0; i < len; i++)
		Model.addConstr(newvar[i] <= var);
}


void TriviumCoreBDPT(GRBModel &Model, GRBVar* Var, int loc[], int round)
{
	GRBVar* var_new;
	string var_name[10];
	char var_type[10] = { 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B' };
	double UB[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	double LB[10] = { 0 };

	for (int i = 0; i < 5; i++)
		var_name[i] = "y_" + to_string(round) + '_' + to_string(loc[i]);
	for (int i = 5; i < 9; i++)
		var_name[i] = "z_" + to_string(round) + '_' + to_string(loc[i]) + '_' + to_string(i - 4);
	var_name[9] = "a_" + to_string(round) + '_' + to_string(loc[0]);
	var_new = Model.addVars(LB, UB, LB, var_type, var_name, 10);

	// Model for Copy
	for (int i = 0; i < 4; i++)
	{
		Model.addConstr(Var[loc[i]] <= var_new[i] + var_new[5 + i], "Copy1 Constr_" + to_string(round) + '_' + to_string(i));
		Model.addConstr(Var[loc[i]] >= var_new[5 + i], "Copy2 Constr_" + to_string(round) + '_' + to_string(i));
		Model.addConstr(Var[loc[i]] >= var_new[i], "Copy3 Constr_" + to_string(round) + '_' + to_string(i));
	}

	// Model for And
	Model.addConstr(var_new[9] == var_new[7], "And1 Constr_" + to_string(round) + "_4");
	Model.addConstr(var_new[9] == var_new[8], "And2 Constr_" + to_string(round) + "_4");

	// Model for Xor
	Model.addConstr(var_new[4] == Var[loc[4]] + var_new[9] + var_new[5] + var_new[6], "Xor Constr_" + to_string(round) + "_5");

	// replace variables 
	for (int i = 0; i < 5; i++)
		Var[loc[i]] = var_new[i];
}


void Numeric_Mapping(int State[], int SState[], Ivterm& Cube, int Round)
{
    int dA[1152+93] = {0};
    int dB[1152+84] = {0};
    int dC[1152+111] = {0};

    int dAA[111] = { 0 };
    int dBB[111] = { 0 };
    int dCC[111] = { 0 };

    for (int i = 80; i < 93; i++)
        dA[92 - i] = -9999;
    for (int i = 0; i < 84; i++)
        dB[83 - i] = -9999;
    for (int i = 0; i < 108; i++)
        dC[110 - i] = -9999;

    for (int i = 0; i < 80; i++)
        if (Cube.test(i))
            dB[83 - i] = 1;

    for (int t = 1; t <= Round; t++)
    {
        int lA[3] = {dA[(t-1)+27], dA[t-1], dB[(t-1)+6]};
        int lB[3] = {dB[(t-1)+15], dB[t-1], dC[(t-1)+24]};
        int lC[3] = {dC[(t-1)+45], dC[t-1], dA[(t-1)+24]};

        int dA_arr[2] = {Deg_Mul(dA, dB, dC, t, 0), maxint(lC, 3)};
        dA[93 + t - 1] = maxint(dA_arr, 2);
        int dB_arr[2] = {Deg_Mul(dA, dB, dC, t, 1), maxint(lA, 3)};
        dB[84 + t - 1] = maxint(dB_arr, 2);
        int dC_arr[2] = {Deg_Mul(dA, dB, dC, t, 2), maxint(lB, 3)};
        dC[111 + t - 1] = maxint(dC_arr, 2);

        if (t == Round)
        {
            dAA[0] = dB_arr[0];
            dBB[0] = dC_arr[0];
            dCC[0] = dA_arr[0];
        }
    }
    for (int i = 0; i < 93; i++)
        State[i] = dA[Round + 92 - i];
    for (int i = 0; i < 84; i++)
        State[93 + i] = dB[Round + 83 - i];
    for (int i = 0; i < 111; i++)
        State[177 + i] = dC[Round + 110 - i];

    for (int t = Round+1; t <= Round+110; t++)
    {
        int lA[3] = { dA[(t - 1) + 27], dA[t - 1], dB[(t - 1) + 6] };
        int lB[3] = { dB[(t - 1) + 15], dB[t - 1], dC[(t - 1) + 24] };
        int lC[3] = { dC[(t - 1) + 45], dC[t - 1], dA[(t - 1) + 24] };

        int dA_arr[2] = { Deg_Mul(dA, dB, dC, t, 0), maxint(lC, 3) };
        dA[93 + t - 1] = maxint(dA_arr, 2);
        int dB_arr[2] = { Deg_Mul(dA, dB, dC, t, 1), maxint(lA, 3) };
        dB[84 + t - 1] = maxint(dB_arr, 2);
        int dC_arr[2] = { Deg_Mul(dA, dB, dC, t, 2), maxint(lB, 3) };
        dC[111 + t - 1] = maxint(dC_arr, 2);

        dAA[t - Round] = dB_arr[0];
        dBB[t - Round] = dC_arr[0];
        dCC[t - Round] = dA_arr[0];
    }

    for (int i = 0; i < 92; i++)
        SState[91 - i] = dAA[i];
    for (int i = 0; i < 83; i++)
        SState[93 + 82 - i] = dBB[i];
    for (int i = 0; i < 110; i++)
        SState[177 + 109 - i] = dCC[i];
    SState[92] = State[92] + State[93];
    SState[176] = State[176] + State[177];
}


int Deg_Mul(int dA[], int dB[], int dC[], int t, int flag)
{
    int *dX1, *dX2;
    if (flag == 0)
    {
        dX1 = dC;
        dX2 = dB;
    }
    else if (flag == 1)
    {
        dX1 = dA;
        dX2 = dC;
    }
    else if (flag == 2)
    {
        dX1 = dB;
        dX2 = dA;
    }
    else 
        return 0;

    int nn[3] = {93, 84, 111};
    int pho[3] = {2, 0, 1};

    int pho1 = pho[flag];
    int pho2 = pho[pho1];
    int t1 = t - nn[pho1] + 1;
    
    if (t1 <= 0)
        return dX1[nn[pho1] + t1 - 1] + dX1[nn[pho1] + t1];
    
    int t2 = t1 - nn[pho2] + 2;
    int temp_arr[3] = {dX2[nn[pho2] + t2 - 1] + dX1[nn[pho1] + t1],
                       dX2[nn[pho2] + t2 + 1] + dX1[nn[pho1] + t1 - 1],
                       dX2[nn[pho2] + t2 - 1] + dX2[nn[pho2] + t2] + dX2[nn[pho2] + t2 + 1]};

    int linear_f0[3][3] = {{dC[t1+45], dC[t1], dA[t1+24]}, 
                           {dA[t1+27], dA[t1], dB[t1+6]},
                           {dB[t1+15], dB[t1], dC[t1+24]}};
    int linear_f1[3][3] = {{dC[(t1-1)+45], dC[t1-1], dA[(t1-1)+24]},
                           {dA[(t1-1)+27], dA[t1-1], dB[(t1-1)+6]},
                           {dB[(t1-1)+15], dB[t1-1], dC[(t1-1)+24]}};

    int max_arr[3] = {minint(temp_arr, 3), 
                      maxint(linear_f0[pho1], 3) + dX1[nn[pho1] + t1 - 1], 
                      maxint(linear_f1[pho1], 3) + dX1[nn[pho1] + t1]};
                      
    int d = maxint(max_arr, 3);

    return d;
}

int minint(int arr[], int len)
{
    if (len == 0)
        return 0;
    int temp = arr[0];
    for (int i = 1; i < len; i++)
        if (temp > arr[i])
            temp = arr[i];
    return temp;
}

int maxint(int arr[], int len)
{
    if (len == 0)
        return 0;
    int temp = arr[0];
    for (int i = 1; i < len; i++)
        if (temp < arr[i])
            temp = arr[i];
    return temp;
}


void ExpressRecursivelyForTrivium(const int OutPutRound, const int TargetRound, char file[])
{
    Poly UpdateFunction[8], ZR(6);
    Poly SNull(0), S1(4), S94(4), S178(4), S1S94(16), S1S178(16), S94S178(16), SAll(64);
    InvExpressInit(ZR, S1, S94, S178, S1S94, S1S178, S94S178, SAll);
    UpdateFunction[0] = SNull;
    UpdateFunction[1] = S178;
    UpdateFunction[2] = S94;
    UpdateFunction[3] = S94S178;
    UpdateFunction[4] = S1;
    UpdateFunction[5] = S1S178; 
    UpdateFunction[6] = S1S94;
    UpdateFunction[7] = SAll;

    Poly* OutPutANF = new Poly;
    Poly* InPutANF  = new Poly;
    Poly* TempPolyPt;
    Poly TempPoly(100*SIZEBw);
    (*OutPutANF).SetPolyLen(100* SIZEBw);
    (*InPutANF).SetPolyLen(100* SIZEBw);
    (*InPutANF).PolyCopy(ZR);

    int CountRound = OutPutRound;
    int FinalRound = OutPutRound - TargetRound + 1;

    while((CountRound--) > FinalRound)
    {
        clock_t start = clock();
        ExpressOneRound(*OutPutANF, *InPutANF, UpdateFunction, TempPoly);
        clock_t end = clock();
        TempPolyPt = InPutANF;
        InPutANF = OutPutANF;
        OutPutANF = TempPolyPt;
        cout << "CountRound = " << CountRound << ", time = " << double(end - start) / CLOCKS_PER_SEC << endl;
    }
    ExpressOneRound(*OutPutANF, *InPutANF, UpdateFunction, TempPoly);
    OutPutANF->WriteValueToFile(TargetRound, file, "OutPutANF");
    // OutPutANF->WriteOutputToFile(TargetRound, file, "OutPutANF");
    delete OutPutANF;
    delete InPutANF;
}


void ExpressOneRound(Poly& OutPutANF, Poly& InPutANF, Poly UpdateFunc[], Poly& TempPoly)
{
    OutPutANF.Size = 0;
    for (uint32_t pt = 0; pt < InPutANF.Size; pt++)
    {
        uint8_t flag = ( ((InPutANF.poly[pt].pterm[0]>>31)&1) << 2) | ( ((InPutANF.poly[pt].pterm[2]>>2)&1) << 1 ) | ( (InPutANF.poly[pt].pterm[5]>>14)&1 );
        uint8_t flagweight = ((InPutANF.poly[pt].pterm[0]>>31)&1) + ((InPutANF.poly[pt].pterm[2]>>2)&1) + ( (InPutANF.poly[pt].pterm[5]>>14)&1 );
        for (int i = 0; i < 8; i++)
            InPutANF.poly[pt].pterm[i] = (InPutANF.poly[pt].pterm[i] << 1) | ( (InPutANF.poly[pt].pterm[i+1] >> 31) & 1);
        InPutANF.poly[pt].pterm[8] <<= 1;
        InPutANF.poly[pt].deg -= flagweight;

        if (flag)
        {
            InPutANF.poly[pt].pterm[2] &= 0xfffffff7;
            InPutANF.poly[pt].pterm[5] &= 0xffff7fff;
            PolyMul(OutPutANF, UpdateFunc[flag], InPutANF.poly[pt]);
        }
        else
            OutPutANF.poly[OutPutANF.Size++] = InPutANF.poly[pt];
    }

    OutPutANF.RemoveDup();
}


void InvExpressInit(Poly& ZR, Poly& S1, Poly& S94, Poly& S178, Poly& S1S94, Poly& S1S178, Poly& S94S178, Poly& SAll)
{
    S1.poly[0].deg = S94.poly[0].deg = S178.poly[0].deg = 2;
    S1.poly[1].deg = S94.poly[1].deg = S178.poly[1].deg = 1;
    S1.poly[2].deg = S94.poly[2].deg = S178.poly[2].deg = 1;
    S1.poly[3].deg = S94.poly[3].deg = S178.poly[3].deg = 1;
    S1.Size = S94.Size = S178.Size = 4;

    S1.poly[0].pterm[8] = (1 << (31 ^ 285&0x1f)) ^ (1 << (31 ^ 286&0x1f));    // 286 287 //0x6;
    S1.poly[1].pterm[8] = (1 << (31 ^ 287&0x1f));                             // 288
    S1.poly[2].pterm[7] = (1 << (31 ^ 242&0x1f));                             // 243
    S1.poly[3].pterm[2] = (1 << (31 ^ 68 &0x1f));                              // 69

    S94.poly[0].pterm[2] = (1 << (31 ^ 90 &0x1f)) ^ (1 << (31 ^ 91&0x1f));     // 91 92
    S94.poly[1].pterm[2] = (1 << (31 ^ 92 &0x1f));                             // 93
    S94.poly[2].pterm[2] = (1 << (31 ^ 65 &0x1f));                             // 66
    S94.poly[3].pterm[5] = (1 << (31 ^ 170&0x1f));                            // 171

    S178.poly[0].pterm[5] = (1 << (31 ^ 174&0x1f)) ^ (1 << (31 ^ 175&0x1f));  // 175 176
    S178.poly[1].pterm[5] = (1 << (31 ^ 176&0x1f));                           // 177
    S178.poly[2].pterm[5] = (1 << (31 ^ 161&0x1f));                           // 162
    S178.poly[3].pterm[8] = (1 << (31 ^ 263&0x1f));                           // 264


    PolyMul(S1S94, S1, S94);
    PolyMul(S1S178, S1, S178);
    PolyMul(S94S178, S94, S178);
    PolyMul(SAll, S1, S94, S178);

    ZR.poly[0].pterm[2] = (1 << (31 ^ 65 & 0x1f));   // 66  (2 1 1073741824)
    ZR.poly[1].pterm[2] = (1 << (31 ^ 92 & 0x1f));   // 93  (2 28 8)
    ZR.poly[2].pterm[5] = (1 << (31 ^ 161 & 0x1f));  // 162 (5 1 1073741824)
    ZR.poly[3].pterm[5] = (1 << (31 ^ 176 & 0x1f));  // 171 (5 16 32768)
    ZR.poly[4].pterm[7] = (1 << (31 ^ 242 & 0x1f));  // 243 (7 18 8192)
    ZR.poly[5].pterm[8] = (1 << (31 ^ 287 & 0x1f));  // 288 (8 31 1)
    ZR.poly[0].deg = ZR.poly[1].deg = ZR.poly[2].deg = ZR.poly[3].deg = ZR.poly[4].deg = ZR.poly[5].deg = 1;
    ZR.Size = 6;

#if 0
    char file[] = "Init.txt";
    ZR.write_output(0, file, "ZR");
    S1.write_output(0, file, "S1");
    S94.write_output(0, file, "S94");
    S178.write_output(0, file, "S178");
    S1S94.write_output(0, file, "S1S94");
    S1S178.write_output(0, file, "S1S178");
    S94S178.write_output(0, file, "S94S178");
    SAll.write_output(0, file, "SAll");
#endif
}


void ReadANF(int &TargetRound, Poly& InPutANF, char file[], bool degkflag)
{
    fstream fin;
    string line1, c1;
    uint32_t temp = 0;
	uint32_t tempdegk = 0;

    fin.open(file, ios_base::in);
    getline(fin, line1);
    fin >> c1 >> TargetRound >> temp;
    cout << line1 << endl << c1 << ' ' << TargetRound << ' ' << temp << endl;
    InPutANF.SetPolyLen(temp);
    InPutANF.Size = temp;

    clock_t starttime = clock();
	if (degkflag)
		for (uint32_t pt = 0; pt < InPutANF.Size; pt++)
		{
			fin >> hex >> tempdegk;
			InPutANF.poly[pt].setdegk(tempdegk);
			for (uint8_t pd = 0; pd < 9; pd++)
				fin >> hex >> InPutANF.poly[pt].pterm[pd];
			fin >> hex >> temp;
			InPutANF.poly[pt].deg = temp;
		}
	else
		for (uint32_t pt = 0; pt < InPutANF.Size; pt++)
		{
			for (uint8_t pd = 0; pd < 9; pd++)
				fin >> hex >> InPutANF.poly[pt].pterm[pd];
			fin >> hex >> temp;
			InPutANF.poly[pt].deg = temp;
			InPutANF.poly[pt].setdegk(-1);
		}
        
    clock_t endtime = clock();
    cout << "Load ANF: " << double(endtime - starttime) / CLOCKS_PER_SEC << " s" << endl;

    fin.close();
    return ;
}


void ReadANF(int &TargetRound, Poly& InPutANF, MapPoly& SP, char file[], bool degkflag)
{
	fstream fin;
	string line1, c1, c2, c3;
	uint32_t temp = 0;
	uint32_t tempdegk = 0;

	fin.open(file, ios_base::in);
	getline(fin, line1);
	fin >> c1 >> TargetRound >> temp;
	cout << line1 << endl << c1 << ' ' << TargetRound << ' ' << temp << endl;
	InPutANF.SetPolyLen(temp);
	InPutANF.Size = temp;

	clock_t starttime = clock();
	if (degkflag)
		for (uint32_t pt = 0; pt < InPutANF.Size; pt++)
		{
			fin >> hex >> tempdegk;
			InPutANF.poly[pt].setdegk(tempdegk);
			for (uint8_t pd = 0; pd < 9; pd++)
				fin >> hex >> InPutANF.poly[pt].pterm[pd];
			fin >> hex >> temp;
			InPutANF.poly[pt].deg = temp;
		}
	else
		for (uint32_t pt = 0; pt < InPutANF.Size; pt++)
		{
			InPutANF.poly[pt].setdegk(-1);
			for (uint8_t pd = 0; pd < 9; pd++)
				fin >> hex >> InPutANF.poly[pt].pterm[pd];
			fin >> hex >> temp;
			InPutANF.poly[pt].deg = temp;
		}
	
	fin >> line1 >> c1 >> c2;
	int ssize = 0;
	ssize = stoi(line1);
	getline(fin, line1);
	cout << ssize << ' ' << c1 << ' ' << c2 << ' ' << line1 << endl;


	while (ssize--)
	{
		getline(fin, line1);
		KeyTerm termtemp; termtemp.reset();

		if (line1.length() > MAPPOLYLEN)
		{
			cout << "False!" << endl;
			return;
		}

		int len = MAPPOLYLEN;
		if (line1.length() < MAPPOLYLEN)
			len = line1.length();

		for (int i = 0; i < len; i++)
			termtemp.set(len - 1 - i, line1[i] - '0');
		SP.InSert(termtemp);
	}

	clock_t endtime = clock();
	cout << "Load ANF: " << double(endtime - starttime) / CLOCKS_PER_SEC << " s" << endl;

	fin.close();
	return;
}


void Convert132varsTo80vars(MapPoly &Ret, MapPoly &InSP)
{
	int Loc[52][4] = { {2, 27, 28, 29}, { 3, 28, 29, 30 }, { 4, 29, 30, 31 }, { 5, 30, 31, 32 }, { 6, 31, 32, 33 }, { 7, 32, 33, 34 }, { 8, 33, 34, 35 }, { 9, 34, 35, 36 }, { 10, 35, 36, 37 },
	{ 11, 36, 37, 38 }, { 12, 37, 38, 39 }, { 13, 38, 39, 40 }, { 14, 39, 40, 41 }, { 15, 40, 41, 42 }, { 16, 41, 42, 43 }, { 17, 42, 43, 44 }, { 18, 43, 44, 45 }, { 19, 44, 45, 46 },
	{ 20, 45, 46, 47 }, { 21, 46, 47, 48 }, { 22, 47, 48, 49 }, { 23, 48, 49, 50 }, { 24, 49, 50, 51 }, { 25, 50, 51, 52 }, { 26, 51, 52, 53 }, { 27, 52, 53, 54 }, { 28, 53, 54, 55 },
	{ 29, 54, 55, 56 }, { 30, 55, 56, 57 }, { 31, 56, 57, 58 }, { 32, 57, 58, 59 }, { 33, 58, 59, 60 }, { 34, 59, 60, 61 }, { 35, 60, 61, 62 }, { 36, 61, 62, 63 }, { 37, 62, 63, 64 },
	{ 38, 63, 64, 65 }, { 39, 64, 65, 66 }, { 40, 65, 66, 67 }, { 41, 66, 67, 68 }, { 42, 67, 68, 69 }, { 43, 68, 69, 70 }, { 44, 69, 70, 71 }, { 45, 70, 71, 72 }, { 46, 71, 72, 73 },
	{ 47, 72, 73, 74 }, { 48, 73, 74, 75 }, { 49, 74, 75, 76 }, { 50, 75, 76, 77 }, { 51, 76, 77, 78 }, { 52, 77, 78, 79 }, {53, 78, 79} };

	vector<KeyTerm>* sk = new vector<KeyTerm>[55];

	KeyTerm temp1(0), temp2(0), temp3(0);
	for (int i = 0; i < 51; i++)
	{
		temp1.set(Loc[i][0]); temp2.set(Loc[i][1]); temp2.set(Loc[i][2]); temp3.set(Loc[i][3]);
		sk[i].push_back(temp1); sk[i].push_back(temp2); sk[i].push_back(temp3);
		temp1.reset(); temp2.reset(); temp3.reset();
	}
	temp1.set(Loc[51][0]); temp2.set(Loc[51][1]); temp2.set(Loc[51][2]);
	sk[51].push_back(temp1); sk[51].push_back(temp2);

	temp1.reset(); temp2.reset();
	temp1.set(66); temp2.set(67);
	sk[52].push_back(temp1); sk[52].push_back(temp3);
	sk[53].push_back(temp2); sk[53].push_back(temp3);

	sk[54].push_back(temp1); sk[54].push_back(temp2); sk[54].push_back(temp3);
	temp3.set(66); temp3.set(67);
	sk[54].push_back(temp3);


	int AllThread = std::thread::hardware_concurrency();
	cout << "The total number of threads available to this computer: " << AllThread << endl;
	omp_set_num_threads(AllThread);
	omp_set_nested(1);
	int ThreadNumForLoop = AllThread;

	MapPoly *TempPoly = new MapPoly[ThreadNumForLoop];

	vector<KeyTerm> TempKeyTerms;
	for (MapTerm::iterator pt = InSP.Mapterm.begin(); pt != InSP.Mapterm.end(); pt++)
		if (pt->second & 1)
			TempKeyTerms.push_back(pt->first);

	int SpLen = TempKeyTerms.size();
	vector<KeyTerm> *Arr[56];

	int Count = 0;

#pragma omp parallel for schedule(dynamic,10) num_threads(ThreadNumForLoop)
	for (int pf = 0; pf < SpLen; pf++)
	{
		KeyTerm Temp = TempKeyTerms[pf];

		int count = 0;
		for (int i = 80; i < 132; i++)
			if (Temp.test(i))
			{
				Arr[count++] = sk + (i - 80);
				Temp.reset(i);
			}
		if (Temp.test(66) && Temp.test(67))
		{
			Arr[count++] = sk + 54;
			Temp.reset(66); Temp.reset(67);
		}
		else if (Temp.test(66) && !Temp.test(67))
		{
			Arr[count++] = sk + 52;
			Temp.reset(66);
		}
		else if (!Temp.test(66) && Temp.test(67))
		{
			Arr[count++] = sk + 53;
			Temp.reset(67);
		}

		vector<KeyTerm> tempkeyterm = { Temp };
		Arr[count++] = &tempkeyterm;
		KeyTermsMulti(TempPoly[omp_get_thread_num()], Arr, count);

#pragma omp atomic
		Count += 1;

		if (Count % 1000 == 1)
			cout << "(Count / All): " << Count << " / " << SpLen << ") " << endl;
	}

	for (int pt = 0; pt < ThreadNumForLoop; pt++)
		Ret.Merge(TempPoly[pt]);
}


void KeyTermsMulti(MapPoly& Ret, vector<KeyTerm> *kpolys[], int len)
{
	vector<KeyTerm> InANF(*(kpolys[0])), OutANF;
	vector<KeyTerm> *InPt = &InANF;
	vector<KeyTerm> *OutPt = &OutANF;
	vector<KeyTerm> *TempPt = &InANF;

	for (int i = 1; i < len; i++)
	{
		OutPt->clear();
		for (vector<KeyTerm>::iterator pt = InPt->begin(); pt != InPt->end(); pt++)
			for (vector<KeyTerm>::iterator pd = (*(kpolys[i])).begin(); pd != (*(kpolys[i])).end(); pd++)
				OutPt->push_back(KeyTermsMulti(*pt, *pd));
		TempPt = InPt;
		InPt = OutPt;
		OutPt = TempPt;
	}

	for (vector<KeyTerm>::iterator pt = InPt->begin(); pt != InPt->end(); pt++)
		Ret.InSert(*pt);
}


KeyTerm KeyTermsMulti(const KeyTerm &p1, const KeyTerm &p2)
{
	KeyTerm Ret;
	for (int i = 0; i < 80; i++)
		Ret.set(i, p1.test(i) | p2.test(i));
	return Ret;
}


void quickSortDeg(Term s[], int64_t l, int64_t r)
{
	if (l < r)
	{
		int64_t i = l, j = r;
		Term x = s[l];

		while (i < j)
		{
			while ((i < j) && (Term_Greater(s[j], x)))
				j--;
			if (i < j)
			{
				s[i] = s[j];
				i += 1;
			}
			while ((i < j) && (!Term_Greater(s[i], x)))
				i++;
			if (i < j)
			{
				s[j] = s[i];
				j -= 1;
			}
		}
		s[i] = x;

		quickSortDeg(s, l, i - 1);
		quickSortDeg(s, i + 1, r);
	}
}