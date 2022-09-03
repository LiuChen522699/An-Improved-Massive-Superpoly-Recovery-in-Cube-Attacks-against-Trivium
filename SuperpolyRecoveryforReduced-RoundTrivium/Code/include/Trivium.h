#pragma once
#include <gurobi_c++.h>
#include <bitset>
#include <map>
#include <omp.h>
#include <iomanip>
#include <thread>
#include <ctime>
#include "Poly.h"


#ifndef TRIVIUM_H_
#define TRIVIUM_H_

const int MAPPOLYLEN = 132;
typedef bitset<80> Ivterm;
typedef bitset<MAPPOLYLEN> BitKeyterm;

class KeyTerm
{
private:
    BitKeyterm Keyterm;

public:
    bool operator<(const KeyTerm& t2) const
    {
        if (Keyterm.count() < t2.Keyterm.count())
            return true;
        
        if (Keyterm.count() > t2.Keyterm.count())
            return false;
        
        for (int i = MAPPOLYLEN-1; i >= 0; i--)
        {
            bool vart1 = Keyterm.test(i);
            bool vart2 = t2.Keyterm.test(i);
            if (vart1 < vart2)
                return true;
            else if (vart1 > vart2)
                return false;
        }
        return false;
    }

	KeyTerm operator&(const KeyTerm& t2) const
	{
		KeyTerm ret(0);
		for (int i = 0; i < MAPPOLYLEN; i++)
			ret.set(i, Keyterm[i] & t2.Keyterm[i]);
		return ret;
	}

	KeyTerm And(const BitKeyterm& t2) const
	{
		KeyTerm ret(0);
		for (int i = 0; i < MAPPOLYLEN; i++)
			ret.set(i, Keyterm[i] & t2[i]);
		return ret;
	}

    int count() const
    {
        return Keyterm.count();
    }

    void show() const
    {
		if (!Keyterm.count())
		{
			cout << 1 << endl;
			return;
		}
        for (int i = 0; i < MAPPOLYLEN; i++)
            if (Keyterm.test(i))
                cout << 'k' << i;
		cout << endl;
    }

    bool test(int loc) const
    {
        return Keyterm.test(loc);
    }

    void set(int loc)
    {
        Keyterm.set(loc);
    }
    void set(int loc, bool val)
    {
        Keyterm.set(loc, val);
	}

	void reset()
	{
		Keyterm.reset();
	}

	void reset(int loc)
	{
		Keyterm.reset(loc);
	}

    KeyTerm()
    {
        Keyterm = 0;
    }

    KeyTerm(int loc)
    {
        if (loc < MAPPOLYLEN)
            Keyterm.set(loc);
    }

    KeyTerm(int Locs[], int len)
    {
        for (int i = 0; i < len; i++)
            if (Locs[i] < MAPPOLYLEN)
                Keyterm.set(Locs[i]);
    }

    KeyTerm(const KeyTerm& t1)
    {
        Keyterm = t1.Keyterm;
    }

    ~KeyTerm()
    {
        Keyterm = 0;
    }

    BitKeyterm get() const
    {
        return Keyterm;
    }
};


typedef pair<const KeyTerm, unsigned int> MapPair;
typedef map<KeyTerm, unsigned int> MapTerm;


class MapPoly
{
public:
	MapTerm Mapterm;

	MapPoly()
	{
		Mapterm.clear();
	}

	MapPoly(MapPoly& p1)
	{
		Mapterm.insert(p1.Mapterm.begin(), p1.Mapterm.end());
	}

	~MapPoly()
	{
		Mapterm.clear();
	}

	bool InSert(const MapPair& t1)
	{
		MapTerm::iterator it;
		it = Mapterm.find(t1.first);
		if (it != Mapterm.end())
		{
			it->second += t1.second;
			return false;
		}
		else
		{
			Mapterm.insert(t1);
			return true;
		}
	}

	bool InSert(const KeyTerm& t1)
	{
		MapTerm::iterator it;
		it = Mapterm.find(t1);
		if (it != Mapterm.end())
		{
			it->second += 1;
			return false;
		}
		else
		{
			Mapterm.insert(MapPair(t1, 1));
			return true;
		}
	}

	bool Merge(MapPoly& p1)
	{
		MapTerm::iterator pt;
		for (pt = p1.Mapterm.begin(); pt != p1.Mapterm.end(); pt++)
			this->InSert(*pt);
		return true;
	}

	bool show()
	{
		cout << "The Constent of this dictionary: " << endl;
		if (!Mapterm.size())
		{
			cout << "NULL!";
			return false;
		}

		MapTerm::iterator pt;
		for (pt = Mapterm.begin(); pt != Mapterm.end(); pt++)
		{
			(pt->first).show();
			cout << " => " << pt->second << endl;
		}
		return true;
	}

	bool showANF(string message = "")
	{
		cout << message << endl << "# p(key) = ";

		if (Mapterm.size() == 0)
		{
			cout << "0 " << endl << "- The number of terms which are contained in superpoly: 0"
				<< endl << "- The number of feasible trails is 0" << endl;
			return false;
		}

		int CoutCount = 0;
		uint64_t FeasibleSolutions = 0;
		uint64_t TrueTermsNum = 0;

		for (MapTerm::reverse_iterator pr = Mapterm.rbegin(); pr != Mapterm.rend(); pr++)
		{
			FeasibleSolutions += pr->second;
			TrueTermsNum += ((pr->second) & 1);

			if (((pr->second) & 1) == 0)
				continue;

			if (CoutCount)
				cout << " + ";
			if ((pr->first.count()) == 0)
				cout << '1';
			else
				for (int i = 0; i < MAPPOLYLEN; i++)
					if (pr->first.test(i))
						cout << 'k' << i + 1;

			CoutCount += 1;
		}
		cout << endl << "- The number of terms which are contained in superpoly: " << TrueTermsNum
			<< endl << "- The number of feasible trails if " << FeasibleSolutions << endl;
		return true;
	}

	bool empty()
	{
		return Mapterm.empty();
	}

	bool clear()
	{
		Mapterm.clear();
		return true;
	}

	int size()
	{
		return Mapterm.size();
	}


	void WritePolyToFile(char file[], string message)
	{
		fstream fout;
		fout.open(file, ios_base::app);

		fout << "- The superpoly of the cube ( With " << MAPPOLYLEN << " vars; " << message << ") ";

		if (Mapterm.size() == 0)
		{
			fout << "with degree = -1, number of terms = 0; "
				<< endl << "# p(key) = 0 " << endl << endl;
			fout.close();
			return;
		}

		fout << "with degree = " << int(Mapterm.rbegin()->first.count())
			<< ", number of terms = " << Mapterm.size() << endl << "# p(key) = ";


		int CoutCount = 0;
		uint64_t FeasibleSolutions = 0;
		uint64_t TrueTermsNum = 0;

		for (MapTerm::reverse_iterator pr = Mapterm.rbegin(); pr != Mapterm.rend(); pr++)
		{
			FeasibleSolutions += pr->second;
			TrueTermsNum += ((pr->second) & 1);

			if (((pr->second) & 1) == 0)
				continue;

			if (CoutCount)
				fout << " + ";
			if ((pr->first.count()) == 0)
				fout << '1';
			else
				for (int i = 0; i < MAPPOLYLEN; i++)
					if (pr->first.test(i))
						fout << 'k' << i;

			CoutCount += 1;
		}
		fout << endl << "- The number of terms which are contained in superpoly: " << TrueTermsNum
			<< endl << "- The number of feasible trails if " << FeasibleSolutions << endl;

#if 0 
		fout << "------------------------------------------------\n";
		for (MapTerm::reverse_iterator pr = Mapterm.rbegin(); pr != Mapterm.rend(); pr++)
		{
			fout << "- " << setw(8) << pr->second;
			fout << "  |    " << (pr->second & 0x1) << "    |  ";
			if (pr->first.get() == 0)
				fout << 1 << endl;
			else
			{
				for (int i = 0; i < MAPPOLYLEN; i++)
					if (pr->first.test(i))
						fout << 'k' << i;
				fout << endl;
			}
		}
		fout << "------------------------------------------------\n";
#endif

		fout << endl;
		fout.close();
	}


	void WritePolyValueToFile(char file[], string message, int flag = 0)
	{
		fstream fout;
		fout.open(file, ios_base::app);
		if (Mapterm.size() == 0)
		{
			fout << 0 << ' ' << message << endl;
			fout.close();
			return;
		}

		int TermNum = 0, AllMonomials = 0;
		for (MapTerm::reverse_iterator pr = Mapterm.rbegin(); pr != Mapterm.rend(); pr++)
		{
			AllMonomials += pr->second;
			if (pr->second & 0x1)
				TermNum += 1;
		}

		if (flag)
			fout << TermNum << ' ' << message << " All Trails: " << AllMonomials << endl;
		else
			fout << TermNum << ' ' << message << endl;

		for (MapTerm::reverse_iterator pr = Mapterm.rbegin(); pr != Mapterm.rend(); pr++)
			if (pr->second & 0x1)
				fout << (pr->first.get()) << endl;

		fout << endl;
		fout.close();
	}


	void ReadPolyValueFromFile(char file[])
	{
		Mapterm.clear();

		fstream fin;
		string line;
		int ssize = 0;

		fin.open(file, ios_base::in);
		fin >> ssize;
		getline(fin, line);
		cout << ssize << ' ' << line << endl;

		while (ssize--)
		{
			getline(fin, line);
			KeyTerm termtemp; termtemp.reset();

			if (line.length() > MAPPOLYLEN)
			{
				cout << "False!" << endl;
				return;
			}

			int len = MAPPOLYLEN;
			if (line.length() < MAPPOLYLEN)
				len = line.length();

			for (int i = 0; i < len; i++)
				termtemp.set(len - 1 - i, line[i] - '0');
			this->InSert(termtemp);
		}
		this->showANF();
		cout << "Read All Date!" << endl;
	}


	void ReadAndWrite80VarsPolyValueFromFile(char file[], char ofile[])
	{
		MapPoly Temp; Temp.clear();
		Temp.ReadPolyValueFromFile(file);

		for (MapTerm::reverse_iterator pr = Temp.Mapterm.rbegin(); pr != Temp.Mapterm.rend(); pr++)
			if (pr->second & 0x1)
				this->InSert(pr->first);

		bool OutPutTermNumFlag = true;
		fstream fout;
		fout.open(ofile, ios_base::app);

		fout << "- The superpoly of the cube ( 80 vars ) ";

		if (Mapterm.size() == 0)
		{
			fout << "with degree = -1, number of terms = 0; "
				<< endl << "# p(key) = 0 " << endl << endl;
			fout.close();
			return;
		}

		fout << "with degree = " << int(Mapterm.rbegin()->first.count())
			<< ", number of terms = " << Mapterm.size() << endl << "# p(key) = ";


		int CoutCount = 0;
		uint64_t TrueTermsNum = 0;

		for (MapTerm::reverse_iterator pr = Mapterm.rbegin(); pr != Mapterm.rend(); pr++)
		{
			TrueTermsNum += ((pr->second) & 1);

			if (((pr->second) & 1) == 0)
				continue;

			if (CoutCount)
				fout << " + ";
			if ((pr->first.count()) == 0)
				fout << '1';
			else
				for (int i = 0; i < 80; i++)
					if (pr->first.test(i))
						fout << 'k' << i;

			CoutCount += 1;
		}
		fout << endl << "- The number of terms which are contained in superpoly: " << TrueTermsNum << endl << endl;

		if (OutPutTermNumFlag)
		{
			fout << "------------------------------------------------\n";
			for (MapTerm::reverse_iterator pr = Mapterm.rbegin(); pr != Mapterm.rend(); pr++)
			{
				fout << "- " << setw(8) << pr->second << "  |    "; //  << (pr->second & 0x1) << "    |  ";
				if (pr->first.get() == 0)
					fout << 1 << endl;
				else
				{
					for (int i = 0; i < 80; i++)
						if (pr->first.test(i))
							fout << 'k' << i;
					fout << endl;
				}
			}
			fout << "------------------------------------------------" << endl;
		}

		fout << endl;
		fout.close();
	}


	void ReadAndWrite132VarsPolyValueFromFile(char file[], char ofile[])
	{
		MapPoly Temp; Temp.clear();
		Temp.ReadPolyValueFromFile(file);

		for (MapTerm::reverse_iterator pr = Temp.Mapterm.rbegin(); pr != Temp.Mapterm.rend(); pr++)
			if (pr->second & 0x1)
				this->InSert(pr->first);

		bool OutPutTermNumFlag = true;

		fstream fout;
		fout.open(ofile, ios_base::app);
		fout << "- The superpoly of the cube ( 132 vars ) ";

		if (Mapterm.size() == 0)
		{
			fout << "with degree = -1, number of terms = 0; "
				<< endl << "# p(key) = 0 " << endl << endl;
			fout.close();
			return;
		}

		fout << "with degree = " << int(Mapterm.rbegin()->first.count())
			<< ", number of terms = " << Mapterm.size() << endl << "# p(key) = ";

		string sk[52];
		sk[0] = "(k2 + k27*k28 + k29)";
		sk[1] = "(k3 + k28*k29 + k30)";
		sk[2] = "(k4 + k29*k30 + k31)";
		sk[3] = "(k5 + k30*k31 + k32)";
		sk[4] = "(k6 + k31*k32 + k33)";
		sk[5] = "(k7 + k32*k33 + k34)";
		sk[6] = "(k8 + k33*k34 + k35)";
		sk[7] = "(k9 + k34*k35 + k36)";
		sk[8] = "(k10 + k35*k36 + k37)";
		sk[9] = "(k11 + k36*k37 + k38)";
		sk[10] = "(k12 + k37*k38 + k39)";
		sk[11] = "(k13 + k38*k39 + k40)";
		sk[12] = "(k14 + k39*k40 + k41)";
		sk[13] = "(k15 + k40*k41 + k42)";
		sk[14] = "(k16 + k41*k42 + k43)";
		sk[15] = "(k17 + k42*k43 + k44)";
		sk[16] = "(k18 + k43*k44 + k45)";
		sk[17] = "(k19 + k44*k45 + k46)";
		sk[18] = "(k20 + k45*k46 + k47)";
		sk[19] = "(k21 + k46*k47 + k48)";
		sk[20] = "(k22 + k47*k48 + k49)";
		sk[21] = "(k23 + k48*k49 + k50)";
		sk[22] = "(k24 + k49*k50 + k51)";
		sk[23] = "(k25 + k50*k51 + k52)";
		sk[24] = "(k26 + k51*k52 + k53)";
		sk[25] = "(k27 + k52*k53 + k54)";
		sk[26] = "(k28 + k53*k54 + k55)";
		sk[27] = "(k29 + k54*k55 + k56)";
		sk[28] = "(k30 + k55*k56 + k57)";
		sk[29] = "(k31 + k56*k57 + k58)";
		sk[30] = "(k32 + k57*k58 + k59)";
		sk[31] = "(k33 + k58*k59 + k60)";
		sk[32] = "(k34 + k59*k60 + k61)";
		sk[33] = "(k35 + k60*k61 + k62)";
		sk[34] = "(k36 + k61*k62 + k63)";
		sk[35] = "(k37 + k62*k63 + k64)";
		sk[36] = "(k38 + k63*k64 + k65)";
		sk[37] = "(k39 + k64*k65 + k66)";
		sk[38] = "(k40 + k65*k66 + k67)";
		sk[39] = "(k41 + k66*k67 + k68)";
		sk[40] = "(k42 + k67*k68 + k69)";
		sk[41] = "(k43 + k68*k69 + k70)";
		sk[42] = "(k44 + k69*k70 + k71)";
		sk[43] = "(k45 + k70*k71 + k72)";
		sk[44] = "(k46 + k71*k72 + k73)";
		sk[45] = "(k47 + k72*k73 + k74)";
		sk[46] = "(k48 + k73*k74 + k75)";
		sk[47] = "(k49 + k74*k75 + k76)";
		sk[48] = "(k50 + k75*k76 + k77)";
		sk[49] = "(k51 + k76*k77 + k78)";
		sk[50] = "(k52 + k77*k78 + k79)";
		sk[51] = "(k53 + k78*k79)";

		uint64_t TrueTermsNum = 0;
		bool outputplusflag = false;
		for (MapTerm::reverse_iterator pr = Mapterm.rbegin(); pr != Mapterm.rend(); pr++)
		{
			if ((pr->second & 1) == 0)
				continue;

			if (outputplusflag)
				fout << " + ";
			outputplusflag = true;

			if (!(pr->first).count())
			{
				fout << 1;
				TrueTermsNum += 1;
				continue;
			}

			int solnum = 1;

			for (int i = 0; i < 66; i++)
				if ((pr->first).test(i))
					fout << 'k' << i;
			if ((pr->first).test(66))
			{
				fout << "(k66+1)";
				solnum *= 2;
			}
			if ((pr->first).test(67))
			{
				fout << "(k67+1)";
				solnum *= 2;
			}
			for (int i = 68; i < 80; i++)
				if ((pr->first).test(i))
					fout << 'k' << i;
			for (int i = 80; i < 130; i++)
				if ((pr->first).test(i))
				{
					fout << sk[i - 80];
					solnum *= 3;
				}
			if ((pr->first).test(131))
			{
				fout << sk[51];
				solnum *= 2;
			}

			TrueTermsNum += solnum;
		}

		fout << endl << " TermNum: " << TrueTermsNum << endl << " ------------------------------------- " << endl;

		if (OutPutTermNumFlag)
			for (MapTerm::reverse_iterator pr = Mapterm.rbegin(); pr != Mapterm.rend(); pr++)
			{
				fout << setw(12) << pr->second << "  |    "; // << (pr->second & 0x1) << "    |  ";

				if (!(pr->first).count())
				{
					fout << 1 << endl;
					continue;
				}

				for (int i = 0; i < 66; i++)
					if ((pr->first).test(i))
						fout << 'k' << i;
				if ((pr->first).test(66))
					fout << "(k66+1)";
				if ((pr->first).test(67))
					fout << "(k67+1)";
				for (int i = 68; i < 80; i++)
					if ((pr->first).test(i))
						fout << 'k' << i;
				for (int i = 80; i < 132; i++)
					if ((pr->first).test(i))
						fout << sk[i - 80];
				fout << endl;
			}
		fout << " --------------------------------------------------------------- " << endl << endl << endl;
		fout.close();
	}

};


void TriviumCoreBDPT(GRBModel &Model, GRBVar* Var, int loc[], int round);
int TriviumEvalBDPT(MapPoly& ResPoly, const Term &Pterm, Ivterm &Cube, int Round, bool SolutFlag = true, int TimeLimit = 0, int ThreadNum = 0);
void FilterTermBDPT(MapPoly& RetSuperpoly, Poly& ExpandingANF, int TargetRound, Ivterm &CubeFlag, int OutPutRound, string files, string appfiles);
bool FilterPartTermBDPT(MapPoly& RetSuperpoly, int& InPutRound, Poly& ANF, Ivterm &CubeFlag, string files, string appfiles);
bool FilterPartTermBDPT0(MapPoly& RetSuperpoly, int& InPutRound, Poly& ANF, Ivterm &CubeFlag, string files, string appfiles);
void ModelCopy(GRBModel &Model, GRBVar &var, GRBVar* newvar, int len);
void Initilization(GRBModel &Model, GRBVar* sVar, GRBVar *kVar, GRBVar *vVar);



void Convert132varsTo80vars(MapPoly &Ret, MapPoly &InSP);
void KeyTermsMulti(MapPoly& Ret, vector<KeyTerm> *kpolys[], int len);
KeyTerm KeyTermsMulti(const KeyTerm &p1, const KeyTerm &p2);



void Numeric_Mapping(int State[], int SState[], Ivterm& Cube, int Round);
int Deg_Mul(int dA[], int dB[], int dC[], int t, int flag);
int minint(int arr[], int len);
int maxint(int arr[], int len);



void InvExpressInit(Poly& ZR, Poly& S1, Poly& S94, Poly& S178, Poly& S1S94, Poly& S1S178, Poly& S94S178, Poly& SAll);
void ExpressOneRound(Poly& OutPutANF, Poly& InPutANF, Poly UpdateFunc[], Poly& TempPoly);
void ExpressRecursivelyForTrivium(const int OutPutRound, const int TargetRound, char file[]);
void ReadANF(int &TargetRound, Poly& InPutANF, char file[], bool degkflag = true);
void ReadANF(int &TargetRound, Poly& InPutANF, MapPoly& SP, char file[], bool degkflag = true);



void quickSortDeg(Term s[], int64_t l, int64_t r);


#endif