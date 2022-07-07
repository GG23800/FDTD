#include"heat_fdtd_1d.hpp"

int main(int argc, char* argv[])
{
    unsigned int nnode = 32;
    FDTD_H1D Sim0;
    Sim0.SetNumberOfNode(nnode);
    Sim0.SetNumberOfCalculationLoop(100000);
    Sim0.HMList.print();
    Sim0.SetBoundaryCondition(0,BoundaryCondition::Dirichlet,-12.2f);
    Sim0.SetBoundaryCondition(nnode-1, BoundaryCondition::Dirichlet,23.8f);

    
    Vector<float> lv(nnode);
    lv.SetRandom(-500.f, 200.f);
    Sim0.SetInitialCondition(lv);

    Sim0.SetDeltat(0.1f);
    std::cout << "Sources: ";
    Sim0.GetSources().print();
    std::cout << "Materials: ";
    Sim0.GetMaterials().print();
    std::cout << "Initial conditions: ";
    Sim0.GetInitialCondition().print();

    Sim0.SetAutoConvergenceTestOn();
    Sim0.SetConvergenceTestStep(100);
    Sim0.SetConvergenceTestCondition(0.00001f);
    Sim0.Run();
    std::cout << "result: ";
    Sim0.GetLastSimulationOutput().print();
    lv.linspace(-12.2f,23.8f);
    lv.print();

    return 0;
}
