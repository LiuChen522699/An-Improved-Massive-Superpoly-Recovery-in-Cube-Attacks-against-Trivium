#include "Poly.h"
#include "Trivium.h"

bool RecoverySuperpoly(int OutPutRound, int Cube[], int dim, int count = 0);


int main()
{
	int Cube[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
		27, 28, 29, 30, 32, 34, 36, 38, 40, 42, 45, 47, 49, 51, 53, 55, 57, 60, 62, 64, 66, 68, 70, 72, 77, 75, 79 };


	bool Flag = RecoverySuperpoly(846, Cube, 54, 54);
	// bool Flag = RecoverySuperpoly(848, Cube, 54, 54);

    cout << endl << "OK!" << endl;
    system("pause");
    return 0;
}


bool RecoverySuperpoly(int OutPutRound, int Cube[], int dim, int count)
{
    const int TargetRoundMM = 300;
    string file_str = "ExpandingANF" + to_string(TargetRoundMM) + ".txt";
    char ANFFile[100];
    file_str.copy(ANFFile, file_str.length(), 0);
    *(ANFFile + file_str.length() ) = '\0';

    file_str = "mkdir FilteredANF" + to_string(OutPutRound) + '_' + to_string(count);
    system(file_str.c_str());

	string Outfile = "FilteredANF" + to_string(OutPutRound) + '_' + to_string(count) + "//FilteredANF" + to_string(OutPutRound) + '_' + to_string(count) + ".txt";
	string AppOutfile = "FilteredANF" + to_string(OutPutRound) + '_' + to_string(count) + "//FilteredANF" + to_string(OutPutRound) + '_' + to_string(count) + "(All results).txt";
	
	char ANFFile2[100];
	Outfile.copy(ANFFile2, Outfile.length(), 0);
	*(ANFFile2 + Outfile.length()) = '\0';

	file_str = "FilteredANF" + to_string(OutPutRound) + '_' + to_string(count) + "//TestResult.txt";
	char FileName[100];
	file_str.copy(FileName, file_str.length(), 0);
	*(FileName + file_str.length()) = '\0';

	file_str = "FilteredANF" + to_string(OutPutRound) + '_' + to_string(count) + "//SuperpolyValue.txt";
	char FileName2[100];
	file_str.copy(FileName2, file_str.length(), 0);
	*(FileName2 + file_str.length()) = '\0';


#if 1  
    ExpressRecursivelyForTrivium(OutPutRound, TargetRoundMM, ANFFile);
#endif


	MapPoly Superpoly;
	Poly ExpandingANF(10);
	int TargetRound = 0;
    Ivterm CubeFlag = 0;
    for (int i = 0; i < dim; i++)
        CubeFlag.set(Cube[i]);


#if 1
	ReadANF(TargetRound, ExpandingANF, ANFFile);
	FilterTermBDPT(Superpoly, ExpandingANF, TargetRound, CubeFlag, OutPutRound, Outfile, AppOutfile);
	string Cubestr1 = "- No." + to_string(count) + " Round(" + to_string(TargetRound) + ") Part" + to_string(0) + " Cube(" + to_string(CubeFlag.count()) + "): ";
	for (int i = 0; i < 80; i++)
		if (CubeFlag.test(i))
			Cubestr1 += 'v' + to_string(i);
	Superpoly.WritePolyToFile(FileName, Cubestr1);
	Superpoly.WritePolyValueToFile(FileName2, Cubestr1, true);
	cout << "Number of terms = " << Superpoly.size() << endl;
	Superpoly.clear();
#endif


#if 1
	ReadANF(TargetRound, ExpandingANF, Superpoly, ANFFile2);
	int TempTargetRound = TargetRound;

	if (ExpandingANF.Size > 10000)
		FilterPartTermBDPT0(Superpoly, TempTargetRound, ExpandingANF, CubeFlag, Outfile, AppOutfile);

	bool TimeFlag = FilterPartTermBDPT(Superpoly, TempTargetRound, ExpandingANF, CubeFlag, Outfile, AppOutfile);
	string Cubestr2 = "- No." + to_string(count) + " Round(" + to_string(TempTargetRound) + ") Part" + to_string(0) + " Cube(" + to_string(CubeFlag.count()) + "): ";
	for (int i = 0; i < 80; i++)
		if (CubeFlag.test(i))
			Cubestr2 += 'v' + to_string(i);
	Superpoly.WritePolyToFile(FileName, Cubestr2);
	Superpoly.WritePolyValueToFile(FileName2, Cubestr2);
	Superpoly.clear();
#endif


	cout << "OK!" << endl;
	return true;
}
