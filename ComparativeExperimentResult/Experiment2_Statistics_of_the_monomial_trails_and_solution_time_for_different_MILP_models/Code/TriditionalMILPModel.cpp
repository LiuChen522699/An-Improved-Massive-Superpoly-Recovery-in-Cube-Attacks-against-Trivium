#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <gurobi_c++.h>
#include <bitset>
#include <iomanip>
#include <ctime>
#include <map>

using namespace std;

typedef bitset<80> IvTerm;
const int KTERMLEN = 80;

class kTerm
{
public: 
	bitset<KTERMLEN> kterm;

	kTerm::kTerm()
	{
		kterm = 0;
	}
	kTerm::~kTerm()
	{
		kterm = 0;
	}

	bool kTerm::operator<(const kTerm& t2) const
	{
		if (kterm.count() < t2.kterm.count())
			return true;

		if (kterm.count() > t2.kterm.count())
			return false;

		for (int i = 0; i < KTERMLEN; i++)
		{
			bool vart1 = kterm.test(i);
			bool vart2 = t2.kterm.test(i);
			if (vart1 < vart2)
				return true;
			else if (vart1 > vart2)
				return false;
		}
		return false;
	}
};
typedef pair<const kTerm, unsigned int> MapPair;
typedef map<kTerm, unsigned int> MapTerm;

class SPoly
{
public:
	MapTerm Term;

	SPoly()
	{
		Term.clear();
	}
	~SPoly()
	{
		Term.clear();
	}

	void InSert(const kTerm &p)
	{
		MapTerm::iterator it;
		it = Term.find(p);
		if (it != Term.end())
			it->second += 1;
		else
			Term.insert(MapPair(p, 1));
	}

	void InSert(const MapPair &p)
	{
		MapTerm::iterator it;
		it = Term.find(p.first);
		if (it != Term.end())
			it->second += p.second;
		else
			Term.insert(p);
	}

	void show(string message)
	{
		int solnum = 0;
		for (MapTerm::iterator pt = Term.begin(); pt != Term.end(); pt++)
			solnum += pt->second;

		cout << "SPoly(" << Term.size() << ", " << solnum << "): " << message << endl;
		for (MapTerm::iterator pt = Term.begin(); pt != Term.end(); pt++)
		{
			if (!(pt->first).kterm.count())
			{
				cout << setw(12) << pt->second << " | " << 1 << endl;
				continue;
			}

			cout << setw(12) << pt->second << " | ";
			for (int i = 0; i < KTERMLEN; i++)
				if ((pt->first).kterm.test(i))
					cout << 'k' << i;
			cout << endl;
		}
	}

	void writetofile(char file[], string message = "")
	{
		bool OutPutTermNumFlag = true;
		fstream fout;
		fout.open(file, ios_base::app);


		int solnum = 0;
		for (MapTerm::iterator pt = Term.begin(); pt != Term.end(); pt++)
			solnum += pt->second;

		fout << "SPoly(" << Term.size() << ", " << solnum << "): " << message << endl << "f(k) = ";

		bool outputplusflag = false;

		for (MapTerm::iterator pt = Term.begin(); pt != Term.end(); pt++)
		{
			if ((pt->second & 1) == 0)
				continue;

			if (outputplusflag)
				fout << " + ";
			outputplusflag = true;

			if (!(pt->first).kterm.count())
			{
				fout << 1;
				continue;
			}
			
			for (int i = 0; i < KTERMLEN; i++)
				if ((pt->first).kterm.test(i))
					fout << 'k' << i;
		}

		fout << endl << " ------------------------------------- " << endl;

		if (OutPutTermNumFlag)
			for (MapTerm::iterator pt = Term.begin(); pt != Term.end(); pt++)
			{
				if (!(pt->first).kterm.count())
				{
					fout << setw(12) << pt->second << "  |  " << 1 << endl;
					continue;
				}

				fout << setw(12) << pt->second << "  |  ";
				for (int i = 0; i < 80; i++)
					if ((pt->first).kterm.test(i))
						fout << 'k' << i;
				fout << endl;
			}
		fout << " --------------------------------------------------------------- " << endl << endl << endl;
		fout.close();
	}

};


void ModelCopy(GRBModel &Model, GRBVar &var, GRBVar* newvar, int len);
void TriviumCore(GRBModel &Model, GRBVar* Var, int loc[], int round);
void TriviumEval(int Cube[], const int dim, int TargetRound, int iv_con[]);

int main()
{

	int TargetRound = 784;
	int iv_con[10] = { 0 };

	vector<vector<int>> vCube = {
		{ 2, 5, 6, 9, 10, 13, 21, 23, 25, 27, 29, 32, 34, 36, 38, 40, 42, 44, 45, 48, 51, 53, 55, 57, 59, 63, 65, 68, 73, 78}, // x2
		{ 2,5,6,9,13,16,19,21,23,25,27,29,30,32,34,36,38,40,42,44,45,48,53,57,59,61,65,68,73,75,78 }, // x4
		{ 2,3,4,6,8,9,10,14,16,19,21,22,23,25,28,30,32,36,37,39,41,49,51,56,59,64,68,71,74,76,79 }, // x7 + 1
		{ 2,4,6,8,10,13,15,19,24,28,29,31,32,34,37,38,40,41,44,47,49,51,53,55,57,59,62,68,70,73,76,78 }, // x9
		{ 1,3,6,8,11,14,15,18,22,25,27,29,34,37,40,42,46,48,50,52,55,57,59,61,66,68,69,71,74,79 }, // x11
		{ 1,6,8,10,15,19,20,22,24,26,29,31,34,37,38,40,42,44,47,49,51,53,55,57,59,62,68,70,76,78 }, // x19
		{ 2,4,6,8,10,13,15,19,22,24,28,29,32,34,37,38,40,41,44,47,49,51,53,55,57,59,62,70,73,78 }, // x20
		{ 2,3,4,6,8,9,10,14,16,19,21,22,23,25,28,29,30,34,36,37,39,41,46,47,48,51,56,59,64,68,71,74,79 }, // x21
		{ 2,5,6,9,11,13,16,19,21,23,24,25,27,29,30,32,34,36,38,40,42,44,45,48,51,53,55,57,59,68,73 }, // x22
		{ 2,3,6,8,11,13,14,16,17,20,22,24,27,30,32,35,37,39,42,44,46,47,48,49,50,53,55,59,64,68,70,72,78 }, // x23 + x68
		{ 2,4,6,8,10,12,14,16,19,21,22,23,25,29,30,32,34,39,41,46,47,48,49,51,56,59,64,67,68,71,79 }, // x24
		{ 1,3,6,11,14,16,20,22,24,27,30,32,35,37,39,42,44,46,47,48,49,50,53,55,59,64,68,70,72,78 }, // x25
		{ 2,5,6,9,11,13,16,19,21,23,25,27,29,32,34,36,38,40,42,44,45,48,51,55,57,65,68,70,75,78 }, // x26
		{ 1,3,6,8,11,13,14,17,20,22,24,27,30,32,34,37,39,42,44,46,47,49,50,53,55,57,59,62,68,70,72,78 }, // x35
		{ 1,2,3,6,8,11,14,17,20,22,24,27,32,34,37,39,42,44,46,47,48,50,53,55,57,59,62,64,68,70,72 }, // x37 + 1
		{ 2,3,6,8,11,13,14,16,20,22,24,27,30,32,34,35,37,39,42,44,46,48,50,53,55,57,59,62,68,72,78 }, // x38
		{ 2,4,6,8,10,12,14,16,21,23,25,30,32,34,36,37,39,41,46,47,48,49,51,56,59,64,67,68,71,79 }, // x39 + 1
		{ 2,5,6,10,13,16,21,23,25,27,29,34,36,38,40,42,44,45,48,51,53,55,59,61,65,68,70,73,75,78 }, // x41 + 1
		{ 2,4,6,8,10,12,13,15,19,20,22,24,26,31,34,37,38,42,44,47,49,53,55,57,59,68,70,73,76,78 }, // x43
		{ 4,6,8,10,13,15,19,20,24,26,28,31,34,37,38,40,41,42,44,47,49,51,53,55,57,59,68,70,73,76,78 }, // x44
		{ 1,4,6,8,10,13,15,19,20,22,24,26,28,29,32,34,38,40,41,42,49,51,53,55,57,59,62,68,70,76,78 }, // x47
		{ 2,4,6,8,10,13,15,19,20,22,24,26,28,29,31,34,37,38,40,42,44,53,55,57,59,62,68,70,76,78 }, // x49
		{ 2,4,6,8,13,15,19,20,22,24,26,28,32,34,37,38,40,42,47,49,51,53,55,57,59,62,70,73,76,78 }, // x50
		{ 2,5,9,10,13,16,21,23,25,27,29,30,32,34,36,38,40,42,44,45,48,53,55,57,59,63,65,68,75,78 }, // x51
		{ 2,4,6,8,10,13,15,19,20,24,26,31,34,37,38 ,40,42,44,47,49,51,53,55,57,59,68,70,73,76,78 }, // x52
		{ 2,5,6,9,13,16,23,25,27,29,30,34,36,38,40,42,44,45,48,51,53,55,57,59,61,63,65,68,70,78 }, // x53
		{ 1,4,6,8,10,12,13,15,19,20,22,24,26,28,31,34,38,40,41,42,44,47,49,55,57,59,68,70,73,76,78 }, // x54
		{ 2,4,6,8,10,13,15,19,22,24,26,31,32,37,38,40,42,44,47,49,51,53,55,57,59,68,70,73,76,78 }, // x55
		{ 2,4,6,8,10,13,15,20,22,24,26,28,31,34,37,38,40,42,44,47,49,51,53,55,57,59,62,68,70,73,78 }, // x56
		{ 2,4,6,8,10,13,15,19,20,22,24,26,28,29,31,32,34,37,38,40,42,44,47,51,53,57,59,62,70,73,76,78 }, // x57
		{ 1,4,6,8,10,13,15,20,22,24,26,28,31,32,34,37,38,41,42,47,49,53,55,57,59,68,70,73,76,78 }, // x58
		{ 2,4,6,8,10,12,14,16,19,21,22,23,25,28,29,30,34,36,37,39,41,46,48,51,56,59,64,71,76,79 }, // x59
		{ 1,4,6,8,10,13,15,19,20,22,24,26,28,29,31,34,38,40,42,44,49,51,53,55,57,59,68,70,76,78 }, // x60
		{ 1,4,6,8,10,12,13,15,19,22,24,26,29,31,34,37,38,42,44,47,49,51,53,55,57,59,68,70,73,78 }, // x61
		{ 2,4,6,8,15,19,20,22,24,26,29,31,32,34,37,38,40,42,44,47,49,51,53,55,57,59,68,70,76,78 }, // x62
		{ 2,4,6,8,13,15,19,20,22,24,26,28,31,32,34,40,41,42,47,49,51,53,55,57,59,68,70,73,76,78 }, // x63
		{ 2,4,6,8,10,13,15,19,20,22,24,26,28,29,31,32,34,40,42,44,47,49,53,55,57,59,68,70,76,78 }, // x64
		{ 2,4,6,8,10,12,15,19,20,22,24,26,28,29,32,34,37,40,42,44,47,51,53,55,57,59,62,68,70,78 }, // x65
		{ 2,4,6,8,10,13,15,19,20,22,24,26,28,29,32,34,37,42,44,47,49,51,53,55,57,59,68,70,76,78 }, // x66
		{ 2,3,4,6,8,10,14,16,19,21,22,23,25,28,30,32,34,36,37,39,41,46,51,56,59,64,68,71,74,76,79 }, // x67
		{ 2,6,8,10,13,15,19,20,22,24,26,28,29,31,32,34,37,38,40,41,42,44,51,53,55,57,59,68,70,76,78 }, // x68 + 1
		{ 1,2,6,8,11,13,14,17,20,22,24,27,32,34,37,39,42,44,46,47,48,50,53,55,57,59,62,64,68,70,78 }, // x69
		{ 2,4,6,8,10,12,13,15,19,22,24,28,29,32,34,37,38,40,41,44,47,49,51,53,55,57,65,68,70,73,74,76,78},  // k63k64 + k65 + k38
		{ 2,4,6,8,10,12,13,15,19,24,28,29,32,34,37,40,41,44,47,49,51,53,55,57,59,62,65,70,72,73,74,76,78 },  // k71k72 + k73 + k46
		{ 2,4,6,8,10,12,13,15,19,24,28,29,32,34,37,38,40,41,44,47,49,51,53,55,57,59,68,70,72,73,74,76,78 },  // k73k74 + k75 + k48
	};

	vector<vector<int>> vCube2 = {
		{ 0,2,4,5,6,7,9,11,13,14,15,18,20,22,24,26,32,35,37,39,42,44,46,48,53,55,57,61,68,70,72,79 },  //25
		{ 0,2,4,5,6,7,9,11,13,14,15,18,20,22,24,26,32,35,39,42,44,46,48,52,55,57,62,68,70,74,76,79 },  //25 + x40
		{ 0,2,4,5,6,7,9,11,13,14,15,18,20,22,24,26,32,35,37,39,42,44,46,48,52,53,55,57,61,62,68,70,79 },  //36
		{ 0,2,4,5,7,9,11,13,14,15,18,20,24,26,30,32,35,37,39,40,42,44,46,48,52,53,55,62,68,70,74,79 },  //38
		{ 0,4,5,6,7,9,11,13,14,15,18,20,22,24,26,30,32,35,37,39,40,44,46,48,52,55,57,62,68,70,74,79 },  //42
		{ 0,4,5,6,7,9,11,13,14,18,20,22,24,26,30,35,37,39,40,44,46,48,52,55,57,62,68,70,72,74,76,79 },  //53
		{ 0,2,4,5,6,9,11,13,14,15,17,18,19,20,22,24,26,30,32,35,39,40,44,48,53,55,57,61,62,70,74,76,79 },  //58
		{ 0,5,6,7,9,11,13,17,19,22,24,26,30,32,35,37,39,42,44,46,48,52,53,55,57,61,62,68,72,74,76,79 },  //60
		{ 0,2,4,5,6,7,9,11,13,15,17,18,19,22,24,26,30,32,37,39,42,44,46,52,53,57,61,62,68,74,76,79 },  //62
		{ 0,4,5,7,9,11,13,14,15,17,18,19,20,22,24,26,30,32,35,37,39,40,44,48,53,55,61,68,72,74,76,79 },  //64
		{ 0,4,5,6,7,9,11,13,15,18,20,22,24,26,30,32,35,37,39,40,42,44,46,48,55,57,62,68,70,72,76,79 },  //66
		{ 0,4,5,6,7,9,11,13,14,15,17,19,22,24,26,32,35,37,39,40,42,46,48,52,55,57,62,68,70,74,76,79 },  //67
		{ 0,2,4,6,8,11,13,16,19,21,23,26,28,30,32,34,36,38,40,42,44,46,49,50,53,56,62,64,69,72,74,75,77,79 },  //9 + x34x35 + x36
		{ 0,2,4,6,8,11,13,16,19,21,23,26,28,30,32,34,36,38,40,42,44,46,50,53,56,58,62,64,69,72,74,75,77,79 },  //22 + x47x48 + x49
		{ 0,2,4,6,8,11,13,16,19,21,23,26,28,30,32,34,36,38,40,42,44,46,49,50,53,56,58,59,62,69,71,74,75,79 },  //24 + x49x50 + x51
		{ 0,2,4,6,8,11,13,16,19,21,23,28,30,32,34,36,38,40,42,44,46,49,50,53,56,59,62,64,66,69,72,74,75,77,79 },  //11 + x36x37 + x38
		{ 0,2,4,6,8,11,13,16,19,21,23,26,28,32,34,36,38,40,42,44,46,50,53,56,58,59,62,64,66,69,71,72,74,75,77,79 },  //52 + x77x78 + x79
		{ 0,2,4,6,8,11,13,16,19,21,23,26,28,32,34,36,38,40,42,44,46,48,50,52,53,56,58,59,62,66,69,71,72,74,75,77,79 },  //9 + x34x35 + x36 + x61 + x17x18 + x19
		{ 0,2,4,5,6,7,9,11,13,14,15,18,20,22,24,26,32,35,37,39,42,44,48,52,53,55,57,61,62,68,70,74,79},  // k71k72 + k73 + k46
		{ 0,2,4,5,6,7,9,11,13,15,18,20,22,24,26,30,32,35,37,39,42,44,46,52,53,57,62,68,70,72,74,79 }  // k27x28 + k29 + k2
	};


	for (int i = 0; i < vCube.size(); i++)
	{
		cout << "- No." << i << ", Round: 784, Cube(" << vCube[i].size() << "): ";
		for (vector<int>::iterator pt = vCube[i].begin(); pt != vCube[i].end(); pt++)
			cout << 'v' << *pt << ' ';
		cout << endl;

		int *Cube = new int[vCube[i].size()];
		int pf = 0;
		for (vector<int>::iterator pt = vCube[i].begin(); pt != vCube[i].end(); pt++)
			Cube[pf++] = *pt;
		TriviumEval(Cube, int(vCube[i].size()), 784, iv_con);
		delete[] Cube;
	}

	for (int i = 0; i < vCube2.size(); i++)
	{
		cout << "- No." << i << ", Round: 799, Cube(" << vCube[i].size() << "): ";
		for (vector<int>::iterator pt = vCube[i].begin(); pt != vCube[i].end(); pt++)
			cout << 'v' << *pt << ' ';
		cout << endl;

		int *Cube = new int[vCube2[i].size()];
		int pf = 0;
		for (vector<int>::iterator pt = vCube2[i].begin(); pt != vCube2[i].end(); pt++)
			Cube[pf++] = *pt;
		TriviumEval(Cube, int(vCube2[i].size()), 799, iv_con);
		delete[] Cube;
	}


    system("pause");
    return 0;
}



void TriviumEval(int Cube[], const int dim, int TargetRound, int iv_con[])
{

	char cfile[] = "SuperpolyRecovery(TraditionModel).txt";
	int Round = TargetRound;
	GRBEnv env = GRBEnv();
	GRBModel Model = GRBModel(env);
	Model.set(GRB_IntParam_LogToConsole, 0);
	Model.set(GRB_IntParam_PoolSearchMode, 2);
	Model.set(GRB_IntParam_PoolSolutions, 2000000000);
	Model.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);
	// Model.set(GRB_IntParam_Threads, 11);
	// Model.set(GRB_DoubleParam_TimeLimit, 24 * 3600);

	IvTerm CubeF(0), Ivterm(0);
	for (int i = 0; i < 80; i++)
		if ((iv_con[i >> 3] >> (7 ^ 7 & i)) & 1)
			Ivterm.set(i);
	for (int i = 0; i < dim; i++)
		CubeF.set(Cube[i]);
	Ivterm |= CubeF;
	Ivterm.flip();

	string sCube = "Round: " + to_string(TargetRound) + ", cube(" + to_string(dim) + "): ";
	for (int i = 0; i < dim; i++)
		sCube += "v" + to_string(Cube[i]);

	// Trivium ????????????bit??????????
	int Loc1[5] = { 65, 170, 90, 91, 92 };
	int Loc2[5] = { 161, 263, 174, 175, 176 };
	int Loc3[5] = { 242, 68, 285, 286, 287 };

	// ????????, ??????????????????????????, ????????????
	GRBVar* kvar = Model.addVars(80, 'B');
	GRBVar* vvar = Model.addVars(80, 'B');
	GRBVar* svar = Model.addVars(288, 'B');

	for (int i = 0; i < 80; i++)
		if (CubeF.test(i))
			Model.addConstr(vvar[i] == 1, "Iv_constr_" + to_string(i));
		else if (Ivterm.test(i))
			Model.addConstr(vvar[i] == 0, "Iv_constr_" + to_string(i));

	GRBLinExpr key_sum;
	for (int i = 0; i < KTERMLEN; i++)
		key_sum += kvar[i];
	Model.setObjective(key_sum, GRB_MAXIMIZE);

	for (int i = 0; i < 80; i++)
		Model.addConstr(svar[i] == kvar[i]);
	for (int i = 80; i < 93; i++)
		Model.addConstr(svar[i] == 0);
	for (int i = 0; i < 80; i++)
		Model.addConstr(svar[93 + i] == vvar[i]);
	for (int i = 80; i < 84; i++)
		Model.addConstr(svar[93 + i] == 0);
	for (int i = 0; i < 108; i++)
		Model.addConstr(svar[177 + i] == 0);

	// ????????
	for (int r = 0; r < Round; r++)
	{
		TriviumCore(Model, svar, Loc1, r);
		TriviumCore(Model, svar, Loc2, r);
		TriviumCore(Model, svar, Loc3, r);
		GRBVar temp = svar[287];
		for (int i = 287; i > 0; i--)
			svar[i] = svar[i - 1];
		svar[0] = temp;
	}

	for (int i = 0; i < 288; i++)
	{
		if ((i == 65) | (i == 92) | (i == 161) | (i == 176) | (i == 242) | (i == 287))
			continue;
		else
			Model.addConstr(svar[i] == 0, "fin_round_con_" + to_string(i));
	}
	Model.addConstr(svar[65] + svar[92] + svar[161] + svar[176] + svar[242] + svar[287] == 1, "final_constr");


	Model.update();
	Model.optimize();
	SPoly SP;
	sCube += ", RunningTime: " + to_string(Model.get(GRB_DoubleAttr_Runtime)) + " s";

	if (Model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
	{
		int nSolutions = Model.get(GRB_IntAttr_SolCount);
		cout << "- Number of solutions found : " << nSolutions << ", RunTime: " << Model.get(GRB_DoubleAttr_Runtime) << 's' << endl;

		for (int res = 0; res < nSolutions; res++)
		{
			kTerm p_term;
			Model.set(GRB_IntParam_SolutionNumber, res);

			for (int loc = 0; loc < KTERMLEN; loc++)
			{
				if (kvar[loc].get(GRB_DoubleAttr_Xn) >= 0.5)
					p_term.kterm.set(loc);
			}
			SP.InSert(p_term);
		}
		SP.show(sCube);
	}
	else
	{
		cout << "Infeasible System! RunTime: " << Model.get(GRB_DoubleAttr_Runtime) << 's' << endl;
	}

	SP.writetofile(cfile, sCube);
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


void TriviumCore(GRBModel &Model, GRBVar* Var, int loc[], int round)
{
    GRBVar* var_new;
    string var_name[10];
    char var_type[10] = {'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B'};
    double UB[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
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
        Model.addConstr(Var[loc[i]] <= var_new[i] + var_new[5 + i] , "Copy1 Constr_" + to_string(round) + '_' + to_string(i));
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


