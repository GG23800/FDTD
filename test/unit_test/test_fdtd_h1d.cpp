#include<unistd.h>
#include"heat_fdtd_1d.hpp"

int main(int argc, char* argv[])
{
    unsigned int nnode = 32;
    FDTD_H1D Sim0;
    Sim0.SetNumberOfNode(nnode);
    Sim0.SetNumberOfCalculationLoop(50000);
    Sim0.HMList.print();
    Sim0.SetBoundaryCondition(0,BoundaryCondition::Dirichlet,-12.2f);
    Sim0.SetBoundaryCondition(nnode-1, BoundaryCondition::Dirichlet,23.8f);

    
    Vector<float> lv(nnode);
    lv.SetRandom(-500.f, 200.f);
    Sim0.SetInitialCondition(lv);

    Sim0.SetDeltat(0.51f);
    std::cout << "Sources: ";
    Sim0.GetSources().print();
    std::cout << "Materials: ";
    Sim0.GetMaterials().print();
    std::cout << "Initial conditions: ";
    Sim0.GetInitialCondition().print();

    // Convergence
    Sim0.SetAutoConvergenceTestOff();
    Sim0.SetConvergenceTestStep(100);
    Sim0.SetConvergenceTestCondition(0.001f);
    // Periodic save
    Sim0.SetPeriodicSaveState(true);
    Sim0.SetPeriodicSaveStep(10);
    Sim0.SetPeriodicSaveFileName("TestPSave");

    Sim0.Run();
    int outref = 0;
    while ( Sim0.IsPreparing() || Sim0.IsRunning() )
    {
        if (Sim0.IsPreparing())
        {
            std::cout << "main is preparing loop: " << outref++ << std::endl;
        }
        //if (Sim0.IsRunning())
        //{
        //    std::cout << "main is running loop: " << outref++ << std::endl;
        //}
        Sim0.GetLastSimulationOutput();//.print();
        usleep(100);
    }/*
    outref = 0;
    while ( Sim0.IsRunning() )
    {
        std::cout << "main is running loop: " << outref++ << std::endl;
        Sim0.GetLastSimulationOutput().print();
        usleep(100);
    }*/
    std::cout << "result: ";
    Sim0.GetLastSimulationOutput().print();
    lv.linspace(-12.2f,23.8f);
    std::cout << "target output: ";
    lv.print();

    return 0;
}
