#include "physicalstuff.h"
#include <stdio.h>

PhysicalStuff *Take;

void newValues(int);
void save(std::vector <double>, int);

int main()
{
    Take = new PhysicalStuff();

    //newValues(0);
    //save(Take->IVCurve(), Take->mode);
    
    newValues(1);
    save(Take->IVCurve(), 7); //Take->mode);
    
    Take->~PhysicalStuff();
    
    return 0;
}

void newValues(int mode)
{
    Take->Te            = 10;           //[eV]
    Take->n             = 2e17;        //[m^-3]
    Take->I             = 2.0;           //[A]
    Take->Ee            = 1000;        //[V]
    Take->inf           = 3000;        //[V]
    Take->ExtraI        = 0;           //[A]
    Take->ExtraEe       = 500;         //[V]
    Take->R             = 1 - 0.01*20; //[%]
    Take->filmThickness = 35;          //[nm]
    Take->mode          = mode;        //'1' - with tunnel amplifier
    Take->SetFECurrentLimit();
    return;
}

void save(std::vector <double> XY, int mode)
{
    std::vector <double> X, Y;
    for (int i=0; i<XY.size(); i+=2)
    {
        X.push_back(XY[i]);
        Y.push_back(XY[i+1]);
    }
    
    FILE *fp;
    char name[] = "modeX.txt";
    name[4] = mode + '0';
    fp = fopen(name, "w");
    for (int i = 0; i < X.size(); i++)
        fprintf(fp, "%f\t%f\n", X[i], Y[i]);
    fclose(fp);
    return;
}
