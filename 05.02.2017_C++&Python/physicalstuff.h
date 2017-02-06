#ifndef PHYSICALSTUFF_H
#define PHYSICALSTUFF_H

#include <vector>

class PhysicalStuff
{
public:
    PhysicalStuff();
    ~PhysicalStuff() {}

public:
    double Te, I, n, Ee, ExtraI, ExtraEe, R, filmThickness, mode, inf;

private:
    double RightFELimit;
    static const double dE  = 1.0;
    static const double e   = 1.6e-19;
    static const double m_e = 9.11e-31;
    static const double m_i = 1.67e-27;
    static const double S   = 0.0009;
    static const double Vp  = 0;
    static const double WorkFunction_Al2O3 = 4.28;

private:
    double dU(double Current);
    double FilmsCurrent(double Phy);

    //----    Maxwellian electrons    ----
    double Electrons(double U);
    double eDist(double E);

    //----    Ions saturation current    ----
    double Ions(double U);

    //----    Electron gun current    ----
    double Gun(double U);
    double gunDist(double E);

    //----    Additional electron component    ----
    double ExtraElectrons(double U);
    double extraDist(double E);

    //----    SEE current    ----
    double SEE(double U); double actualU;
    double seeIntegrand(double U);
    double Sigma(double);

    //----    SIE current    ----
    double SIE(double);
    double Gamma(double);

    //----    FE current    ----
    double FECurrent(double U, double Phy);

private:
    //========    Auxiliary functions    ========
    double Integration(double (PhysicalStuff::*)(double), double LeftLimit, double RightLimit, double Step = dE);

public:
    //========    VAC    ========
    void SetFECurrentLimit();
    std::vector <double> IVCurve();
};

#endif // PHYSICALSTUFF_H
