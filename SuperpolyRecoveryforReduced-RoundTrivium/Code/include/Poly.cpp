#include "Poly.h"


Poly::Poly()
{
	poly = new Term[SIZEBw];
	Size = 0;
}


Poly::Poly(uint32_t size)
{
	poly = new Term[size + 1];
	Size = 0;
}


Poly::Poly(const Poly &p)
{
	Size = p.Size;
	poly = new Term[Size + 1];
	for (uint32_t pd = 0; pd < Size; pd++)
		poly[pd] = p.poly[pd];
}


Poly & Poly::operator=(const Poly &p)
{
	if (this == &p)
		return *this;
    delete[] poly;
    Size = p.Size;
    poly = new Term [Size + 1];

    if (Size != 0)
	    for (uint32_t pd = 0; pd < Size; pd++)
		    poly[pd] = p.poly[pd];
	return *this;
}


Poly & Poly::SetPolyLen(const uint32_t len)
{
    delete [] poly;
    poly = new Term [len +1];
    Size = 0;
    return *this;
}


Poly::~Poly()
{
	delete [] poly;
}


Poly & Poly::PolyCopy(const Poly &p)
{
	if (this == &p)
		return *this;
	Size = p.Size;
	for (uint32_t pd = 0; pd < Size; pd++)
		poly[pd] = p.poly[pd];
	return *this;
}


Poly & Poly::RemoveDup()
{
    if (this->Size <= 1)
        return *this;
    
    quickSort(this->poly, int64_t(0), int64_t(this->Size)-1);

    uint32_t pd = 0, pf = 0, TempSize = this->Size;
    while (pd < TempSize-1)
    {
        if (this->poly[pd] == this->poly[pd+1])
            pd += 2;
        else
            this->poly[pf++] = this->poly[pd++];
    }
    if (pd == TempSize-1)
        this->poly[pf++] = this->poly[pd];
    this->Size = pf;

    return *this;
}


Poly & Poly::Merge(const Poly &p)
{
	if (p.Size == 0)
		return *this;

	for (int i = 0; i < p.Size; i++)
		this->poly[this->Size++] = p.poly[i];
	this->RemoveDup();

	return *this;
}


void Poly::WriteOutputToFile(int Round, char file[], string message, bool ios_base_flag)
{
    fstream fout;
	if (ios_base_flag == 0)
		fout.open(file, ios_base::app);
	else
		fout.open(file, ios_base::out);

    if (Size == 0)
    {
        fout << "Round_" << Round << " OutPutFunction: " << endl << "# 0" << message << endl;
        fout.close();
        return ;
    }

    uint32_t length = Size;
    fout << "Round_" << Round << " OutPutFunction with degree = "  
         << int(poly[0].deg) << ", number of terms = "
         << Size << "; " << message << endl << "- Poly = ";
        
    for (uint32_t pr = 0; pr < Size; pr++)
    {
        if (pr)
            fout << " + ";
        for (int i = 0; i < 9; i++)
            for (int j = 0; j < 32; j++)
                if ( (poly[pr].pterm[i] >> (31 ^ j)) & 1 )
                    fout << 's' << ( (i<<5) + j + 1);
    }
    fout << endl << endl;
    fout.close();
}


void Poly::WriteValueToFile(int Round, char file[], string message, bool ios_base_flag)
{
    fstream fout;
	if (ios_base_flag == 0)
		fout.open(file, ios_base::app);
	else
		fout.open(file, ios_base::out);

    if (Size == 0)
    {
		fout << "Round_of_Expanding_Recursively_" << Round << " OutPutFunction with degree = " << int(poly[0].deg)
			<< ", number of terms = " << Size << "; " << message << endl << "- " << Round << ' ' << Size << endl << endl;
        fout.close();
        return ;
    }

    uint32_t length = Size;
    fout << "Round_of_Expanding_Recursively_" << Round << " OutPutFunction with degree = "  << int(poly[0].deg)
         << ", number of terms = " << Size << "; " << message << endl;
    fout << "- " << Round << ' ' << Size << endl;
        
    for (uint32_t pr = 0; pr < Size; pr++)
	{
		fout << hex << "0x" << poly[pr].getdegk() << ' ';
        for (int i = 0; i < 9; i++)
            fout << hex << "0x" << poly[pr].pterm[i] << ' ';
        fout << hex << "0x" << int(poly[pr].deg) << endl;
    }
    fout << endl << endl;
    fout.close();
}



void Poly::show(int Round, string message)
{
    if (Size == 0)
    {
        cout << "Round_" << Round << " OutPutFunction: " << endl << "- 0" << message << endl;
        return ;
    }

    uint32_t length = Size;
    cout << "Round_" << Round << " OutPutFunction with degree = "  
         << int(poly[0].deg) << ", number of terms = "
         << Size << "; " << message << endl << "- Poly = ";
        
    for (uint32_t pr = 0; pr < Size; pr++)
    {
        if (pr)
            cout << " + ";
        for (int i = 0; i < 9; i++)
            for (int j = 0; j < 32; j++)
                if ( (poly[pr].pterm[i] >> (31 ^ j)) & 1 )
                    cout << 's' << ( (i<<5) + j + 1);
    }
    cout << endl;
}


void PolyAdd(Poly &result, const Poly &p1, const Poly &p2)
{
    uint32_t len1 = p1.Size, len2 = p2.Size;
	result.Size = 0;

    if (len1 == 0)
        for (uint32_t pd = 0; pd < len2; pd++)
        {
            result.poly[result.Size] = p2.poly[pd];
            result.Size += 1;
        }
    else if (len2 == 0)
        for (uint32_t pt = 0; pt < len1; pt++)
        {
            result.poly[result.Size] = p1.poly[pt];
            result.Size += 1;
        }
    else if ( (len1 != 0) && (len2 != 0) )
    {
        uint32_t pt = 0; 
        uint32_t pd = 0;
        while( (pt < len1) && (pd < len2) )
        {
            if (p1.poly[pt] > p2.poly[pd])
                result.poly[result.Size++] = p1.poly[pt++];
            else if (p1.poly[pt] < p2.poly[pd])
                result.poly[result.Size++] = p2.poly[pd++];
            else
            {
                pd += 1;
                pt += 1;
            }
        }
        for (; pt < len1; pt++)
            result.poly[result.Size++] = p1.poly[pt];
        for (; pd < len2; pd++)
            result.poly[result.Size++] = p2.poly[pd];
    }
}


void PolyMul(Poly &result, Poly &p1, Term &pt1)
{
    for (uint32_t i = 0; i < p1.Size; i++)
        result.poly[result.Size++] = p1.poly[i] * pt1;
}


void PolyMul(Poly &result, Poly &p1, Poly &p2)
{
	uint32_t len1 = p1.Size;
	uint32_t len2 = p2.Size;
	result.Size = 0;

	if ((len1 == 0) || (len2 == 0))
		return;

	if ((len1 != 0) && (len2 != 0))
		for (uint32_t pd = 0; pd < len2; pd++)
			PolyMul(result, p1, p2.poly[pd]);

	result.RemoveDup();
}


void PolyMul(Poly &result, Poly& p1, Poly& p2, Poly& p3)
{
    uint32_t TempSize = p1.Size * p2.Size;
    Poly TempPoly(TempSize);
    PolyMul(TempPoly, p1, p2);
    PolyMul(result, TempPoly, p3);
}


// Term Order and Operation
Term operator* (const Term &pt1, const Term &pt2)
{
	Term result;
    if (pt1.deg == 0)
        result = pt2;
    else if (pt2.deg == 0)
        result = pt1;
    else
    {
        for (int i = 0; i < 9; i++)
		    result.pterm[i] = (pt1.pterm[i]) | (pt2.pterm[i]);
	    result.deg = degree(result.pterm);
		result.setdegk(-1);
	}
	return result;
}


bool operator< (const Term &p1, const Term &p2)
{
	if (p1.deg < p2.deg)
		return true;

	else if (p1.deg > p2.deg)
		return false;

	else
	{
		for (int i = 0; i < 9; i++)
		{
			if (p1.pterm[i] < p2.pterm[i])
				return true;
			else if (p1.pterm[i] > p2.pterm[i])
				return false;
		}
	}
	return false;
}


bool operator<= (const Term &p1, const Term &p2)
{
	if (p1 > p2)
        return false;
    else
        return true;
}


bool operator> (const Term &p1, const Term &p2)
{
	if (p1.deg > p2.deg)
		return true;

	else if (p1.deg < p2.deg)
		return false;

	else
	{
		for (int i = 0; i < 9; i++)
		{
			if (p1.pterm[i] > p2.pterm[i])
				return true;
			else if (p1.pterm[i] < p2.pterm[i])
				return false;
		}
	}
	return false;
}


bool Term_Greater(const Term &p1, const Term &p2)
{
	int dk1 = p1.getdegk();
	int dk2 = p2.getdegk();

	if (dk1 != dk2)
		return (dk1 > dk2);

	int d1 = 0, d2 = 0, len = min(p1.deg, p2.deg);
	int *listd1 = new int[p1.deg];
	int *listd2 = new int[p2.deg];
	for (int i = 0; i < 93; i++)
	{
		if ((p1.pterm[i >> 5] >> (31 - (i & 31))) & 1)
			listd1[d1++] = i;
		if ((p2.pterm[i >> 5] >> (31 - (i & 31))) & 1)
			listd2[d2++] = i;
	}
	for (int i = 93; i < 177; i++)
	{
		if ((p1.pterm[i >> 5] >> (31 - (i & 31))) & 1)
			listd1[d1++] = i - 93;
		if ((p2.pterm[i >> 5] >> (31 - (i & 31))) & 1)
			listd2[d2++] = i - 93;
	}
	for (int i = 177; i < 288; i++)
	{
		if ((p1.pterm[i >> 5] >> (31 - (i & 31))) & 1)
			listd1[d1++] = i - 177;
		if ((p2.pterm[i >> 5] >> (31 - (i & 31))) & 1)
			listd2[d2++] = i - 177;
	}
	sort(listd1, listd1 + d1);
	sort(listd2, listd2 + d2);

	d1 = d2 = 0;
	for (int i = 0; i < len; i++)
	{
		d1 += listd1[i];
		d2 += listd2[i];
	}
	delete[] listd1;
	delete[] listd2;

	if (d1 < d2)
		return true;
	if (d1 > d2)
		return false;

	if (p1.deg >= p2.deg)
		return true;
	else
		return false;
}


bool operator== (const Term &p1, const Term &p2)
{

	if (p1.deg != p2.deg)
		return false;

	for (int i = 0; i < 9; i++)
		if (p1.pterm[i] != p2.pterm[i])
			return false;

	return true;
}


bool Divisibility(const Term &BigP, const Term &SmallP)
{
    if (BigP.deg < SmallP.deg)
        return false;

    Term temp;
    for (int i = 0; i < 9; i++)
        temp.pterm[i] = BigP.pterm[i] | SmallP.pterm[i];
    if (temp == BigP)
        return true;
    else
        return false;
}


template <class T>
void quickSort(T s[], int64_t l, int64_t r)
{
	if (l < r)
	{      
		int64_t i = l, j = r;
        T x = s[l];
		while (i < j)
		{
			while((i < j) && (s[j] <= x))
				j--; 
			if(i < j)
				s[i++] = s[j];
			while(i < j && s[i] > x)
				i++; 
			if(i < j)
				s[j--] = s[i];
		}
		s[i] = x;
		quickSort(s, l, i - 1);
		quickSort(s, i + 1, r);
	}
}


uint8_t Weight(uint32_t n)
{
	n = n - ((n >> 1) & 0x55555555);
	n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
	n = (n & 0x0f0f0f0f) + ((n >> 4) & 0x0f0f0f0f);
	n = n + (n >> 8);
	n = n + (n >> 16);
	return uint8_t(n & 0x0000003f);
}


uint8_t degree(uint32_t pt[], bool flag)
{
    uint8_t value = 0;
    if (!flag)
        for (int i = 0; i < 9; i++)
            value += Weight(pt[i]);
    else
        value = Weight(pt[0]) + Weight(pt[1]) + Weight(pt[2] & 0xffff0000);
    
    return value;
}