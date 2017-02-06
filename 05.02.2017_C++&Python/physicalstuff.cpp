#include "physicalstuff.h"
#include <cmath>

PhysicalStuff::PhysicalStuff()
{
}

//===========================
//========    VAC    ========
//===========================
std::vector <double> PhysicalStuff::IVCurve()
{
    std::vector <double> result;

    if (!mode)
    {
        for (double U = Vp; U < inf; U++)
        {
            result.push_back(U);
            result.push_back(FilmsCurrent(U));
        }
    }
    else
    {
        for (double Phy = Vp; Phy < inf; Phy += 0.01)
        {
            double U = dU(FilmsCurrent(Phy));
            result.push_back(Phy + U);
            result.push_back(-FECurrent(Phy + U, Phy));
        }
    }

    return result;
}

double PhysicalStuff::dU(double Current)
{
    double U = 0, Step = RightFELimit;

    if (Current > 0)
        U += Step;
    else
        U -= Step;

    for (int i = 0; i < 20; i++)
    {
        Step /= 2.0;
        if (std::fabs(Current + FECurrent(U+Step, 0)*R) < std::fabs(Current + FECurrent(U, 0)*R))
            U += Step;
        if (std::fabs(Current + FECurrent(U-Step, 0)*R) < std::fabs(Current + FECurrent(U, 0)*R))
            U -= Step;
    }

    return U;
}

double PhysicalStuff::FilmsCurrent(double U)
{
    U = U - Vp;
    if (U >= 0)
        return Ions(U) - Electrons(U) - Gun(U) - ExtraElectrons(U) + SIE(U) + SEE(U);
    else
        return FilmsCurrent(Vp);
}

//============================================
//========    Maxwellian electrons    ========
//============================================
double PhysicalStuff::Electrons(double U)
{
    //return S*(e*n/4)*sqrt((8*Te*e)/(M_PI*m_e))*exp(-U/Te);
    return S*e*e*Integration(&PhysicalStuff::eDist, U, Ee);
}

double PhysicalStuff::eDist(double E)
{
    return n*exp(-E/Te)/sqrt(2*M_PI*m_e*Te*e);
}

//==============================================
//========    Ion saturation current    ========
//==============================================
double PhysicalStuff::Ions(double U)
{
    return (sqrt(2)/exp(1))*S*e*n*sqrt(Te*e/m_i);
}

//============================================
//========    Electron gun current    ========
//============================================
double PhysicalStuff::Gun(double U)
{
    return I*Integration(&PhysicalStuff::gunDist, U, Ee);
}

double PhysicalStuff::gunDist(double E)
{
    if ((E > 0) && (E < Ee))
        return 1/(2*sqrt(E*Ee));
    else
        return 0;
    /*if (E < Ee)
        return 1/Ee;
    else
        return 0;*/
}

//=====================================================
//========    Additional electron component    ========
//=====================================================
double PhysicalStuff::ExtraElectrons(double U)
{
    return ExtraI*Integration(&PhysicalStuff::extraDist, U, ExtraEe);
}

double PhysicalStuff::extraDist(double E)
{
    if (E < ExtraEe)
        return 1/ExtraEe;
    else
        return 0;
}

//===================================
//========    SEE current    ========
//===================================
double PhysicalStuff::SEE(double U)
{
    actualU = U;
    return Integration(&PhysicalStuff::seeIntegrand, U, inf);
}

double PhysicalStuff::seeIntegrand(double U)
{
    return Sigma(U - actualU)*(I*gunDist(U) + ExtraI*extraDist(U));// + S*e*e*eDist(U));
}

double PhysicalStuff::Sigma(double Energy)
{
    double SigmaMax = 6.4, EqPower = 1.2, MaxEnergy = 600;
    return SigmaMax*sqrt(Energy/MaxEnergy)*exp(-(pow(Energy, EqPower) - pow(MaxEnergy, EqPower))/(2*EqPower*pow(MaxEnergy, EqPower)));
}

//===================================
//========    SIE current    ========
//===================================
double PhysicalStuff::SIE(double U)
{
    return Ions(U)*Gamma(U);
}

double PhysicalStuff::Gamma(double Energy)
{
    return 0.3*pow(0.001*Energy, 0.53);
}

//==================================
//========    FE current    ========
//==================================
double PhysicalStuff::FECurrent(double U, double Phy)
{
    double E = std::fabs(U - Phy)/(filmThickness*1e-7);
    double u = (4.39/pow(WorkFunction_Al2O3, 0.5)) - (2.82e7*pow(WorkFunction_Al2O3, 1.5)/E);
    double Result = (1e4*S)*1.4e-6*pow(E,2)*pow(10, u)/WorkFunction_Al2O3;

    if (Phy <= U)
        return -Result;
    else
        return Result;
}

void PhysicalStuff::SetFECurrentLimit()
{
    double Phy = 0;
    while (FECurrent(0, Phy) < 10000.0)
        Phy += 0.01;
    RightFELimit = Phy;
}

//===========================================================
//================    Auxiliary functions    ================
//===========================================================
double PhysicalStuff::Integration(double (PhysicalStuff::* Function)(double), double Left, double Right, double Step)
{
    double Integral = 0;
    for(double X=Left; X<Right; X+=Step)
        Integral += 0.5*Step*((this->*Function)(X) + (this->*Function)(X+Step));
    return Integral;
}
