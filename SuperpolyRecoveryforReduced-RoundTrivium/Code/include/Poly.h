#pragma once
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>

using namespace std;

const int SIZEBw = 1000000;

#ifndef POLY_H_
#define POLY_H_


class Term
{
private:
	int degk = -1;

public:
	uint32_t pterm[9] = {0};
	uint8_t deg = 0;

    void write(char file[], string message)
    {
        fstream fout;
        fout.open(file, ios_base::app);
        fout << "- " << message << ": ";

        if (deg == 0)
            fout << '1' << endl;
        else 
            for (int i = 0; i < 9; i++)
                for (int j = 0; j < 32; j++)
                    if ((pterm[i] >> (31 ^ j)) & 1)
                        fout << 's' << ((i << 5) + j + 1);
        fout << endl;
        fout.close();
    }

	void show(string message)
	{
		cout << "- " << message << ": ";

		if (deg == 0)
			cout << '1';
		else
			for (int i = 0; i < 9; i++)
				for (int j = 0; j < 32; j++)
					if ((pterm[i] >> (31 ^ j)) & 1)
						cout << 's' << ((i << 5) + j + 1);
		cout << endl;
	}

    Term& operator=(const Term &p)
    {
        if (this == &p)
            return *this;
        
        for (int i = 0; i < 9; i++)
            this->pterm[i] = p.pterm[i];
        this->deg = p.deg;
		this->degk = p.degk;
        return *this;
    }

	bool operator!=(const Term &p)
	{
		for (int i = 0; i < 9; i++)
			if (pterm[i] != p.pterm[i])
				return true;
		return false;
	}

	void setdegk(int n)
	{
		if (n >= 80)
			degk = 80;
		else
			degk = n;
	}

	int getdegk() const
	{
		return degk;
	}

};

class Poly
{
public:
	Term* poly;
	uint32_t Size = 0;

	Poly();
	Poly(uint32_t size);
	Poly(const Poly &p);
	~Poly();

	Poly & operator=(const Poly &p);	
	Poly & PolyCopy(const Poly &p);
	Poly & SetPolyLen(const uint32_t len);
    Poly & RemoveDup();
	Poly & Merge(const Poly &p);

	void WriteValueToFile(int n, char file[], string message = " ", bool ios_base_flag = 0);
	void WriteOutputToFile(int n, char file[], string message = " ", bool ios_base_flag = 0);
    void show(int Round, string message = " ");
};

Term operator* (const Term &pt1, const Term &pt2);
bool operator< (const Term &p1, const Term &p2);
bool operator<= (const Term &p1, const Term &p2);
bool operator> (const Term &p1, const Term &p2);
bool operator== (const Term &p1, const Term &p2);
bool Term_Greater(const Term &p1, const Term &p2);
bool Divisibility(const Term &BigP, const Term &SmallP);

void PolyAdd(Poly &result, const Poly &p1, const Poly &p2);
void PolyMul(Poly &result, Poly &p1, Term &pt1);
void PolyMul(Poly &result, Poly &p1, Poly &p2);
void PolyMul(Poly &result, Poly& p1, Poly& p2, Poly& p3);



uint8_t Weight(uint32_t n);
uint8_t degree(uint32_t pt[], bool flag = 0);

template <class T>
void quickSort(T s[], int64_t l, int64_t r);

#endif